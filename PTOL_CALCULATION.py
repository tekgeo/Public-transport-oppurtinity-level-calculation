# -*- coding: utf-8 -*-
"""
PTOL for Milan: H3 grid + PT stops + walking isochrones.
Outputs a per-hex score = number of accessible PT modes (0..4).
"""
#@author: Juanita Sanchez

import os
import time

import h3
import networkx as nx
import numpy as np
import osmnx as ox
import pandas as pd
import geopandas as gpd

from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from shapely.geometry import Polygon, LineString, MultiPolygon
from shapely.ops import substring


#paths and CRS 

CRS_WGS84  = "EPSG:4326"
CRS_METRIC = "EPSG:3035"  # LAEA Europe, equal-area, meters

SHAPES_DIR = r"C:\Users\juani\Desktop\2nd sem\Space Econ\Data\shapes"
DATA_DIR   = r"C:\Users\juani\Desktop\2nd sem\Space Econ\PTOL\Data"
OUTPUT_DIR = r"C:\Users\juani\Desktop\2nd sem\Space Econ\PTOL\Outputs"

os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)

MILAN_BOUNDARY    = os.path.join(SHAPES_DIR, "Citta Metropolitana di Milano.shp")
WALK_NETWORK_FILE = os.path.join(DATA_DIR, "citta_metropolitana_milano_walk_network.graphml")

# One point layer per PT mode. Comment out any line for a mode you don't have.
STOPS_FILES = {
    "Bus":     os.path.join(DATA_DIR, "bus_stops.shp"),
    "Tram":    os.path.join(DATA_DIR, "tram_stops.shp"),
    "Metro":   os.path.join(DATA_DIR, "metro_stations.shp"),
    "Railway": os.path.join(DATA_DIR, "railway_stations.shp"),
}


#analysis parameters 

H3_RESOLUTION = 9                # ~0.1 km² per cell
WALK_SPEED_M_PER_MIN = 80.0      # 4.8 km/h, standard PTAL value

# Walking-time thresholds per PT mode (PTAL convention).
MODE_CONFIG = {
    "Bus":     {"max_time_min": 8.0},
    "Tram":    {"max_time_min": 8.0},
    "Metro":   {"max_time_min": 12.0},
    "Railway": {"max_time_min": 12.0},
}

# Isochrone polygon: each reached edge is buffered and unioned.
EDGE_BUFFER_M          = 25      # half-width of the walkable corridor
EDGE_BUFFER_RESOLUTION = 4       # octagon buffer; faster union, no visible loss at 25 m
SMALL_HOLE_AREA_M2     = 5_000   # fill holes <= this (m²) - inner courtyards;
                                 # larger gaps (parks, water) remain holes

# Stop deduplication (clusters of points closer than tolerance get collapsed).
DEDUPLICATE_STOPS = True
DEDUP_TOLERANCE_M = 5

# Stops more than this far from any pedestrian node are dropped as data errors.
MAX_SNAP_DISTANCE_M = 100

# A hex counts as "accessible" if at least this fraction of its area falls
# inside the isochrone. 0 = any intersection, 0.5 = majority.
MIN_HEX_OVERLAP_FRACTION = 0.25

# helpers


def read_vector_checked(filepath, target_crs=None, layer_name="vector layer"):
    gdf = gpd.read_file(filepath)
    if gdf.crs is None:
        raise ValueError(
            f"The layer '{layer_name}' has no CRS defined. "
            "Assign the correct CRS before running this script."
        )
    if target_crs is not None:
        gdf = gdf.to_crs(target_crs)
    return gdf


def union_geometries(gdf_or_geoseries):
    # union_all() is the modern method; fall back to unary_union on older
    # geopandas versions.
    geoms = (
        gdf_or_geoseries.geometry
        if hasattr(gdf_or_geoseries, "geometry")
        else gdf_or_geoseries
    )
    if hasattr(geoms, "union_all"):
        return geoms.union_all()
    return geoms.unary_union


def ensure_point_geometry(gdf):
    # If the layer is polygons (stations as buildings, etc.), use centroids.
    gdf = gdf[gdf.geometry.notna() & ~gdf.geometry.is_empty].copy()
    if gdf.empty:
        return gdf
    if (gdf.geometry.geom_type == "Point").all():
        return gdf
    original_crs = gdf.crs
    gdf_metric = gdf.to_crs(CRS_METRIC)
    gdf_metric["geometry"] = gdf_metric.geometry.centroid
    return gdf_metric.to_crs(original_crs)


def build_h3_grid_from_polygon(polygon, resolution):
    geojson = polygon.__geo_interface__
    hex_ids = set(h3.geo_to_cells(geojson, resolution))
    records = []
    for hid in hex_ids:
        clat, clon = h3.cell_to_latlng(hid)
        boundary = h3.cell_to_boundary(hid)
        records.append({
            "h3_index": hid,
            "clat": clat,
            "clon": clon,
            "geometry": Polygon([(lng, lat) for lat, lng in boundary]),
        })
    return gpd.GeoDataFrame(records, crs=CRS_WGS84)


def meters_to_seconds(distance_m):
    try:
        value = float(distance_m)
    except (TypeError, ValueError):
        return np.nan
    if np.isnan(value) or np.isinf(value) or value < 0:
        return np.nan
    return value / (WALK_SPEED_M_PER_MIN / 60.0)


def load_graphml_compat(filepath):
    # OSMnx changed the location of load_graphml between v1 and v2.
    if hasattr(ox, "io") and hasattr(ox.io, "load_graphml"):
        return ox.io.load_graphml(filepath)
    return ox.load_graphml(filepath)


def project_graph_compat(G, target_crs):
    if hasattr(ox, "projection") and hasattr(ox.projection, "project_graph"):
        return ox.projection.project_graph(G, to_crs=target_crs)
    if hasattr(ox, "project_graph"):
        return ox.project_graph(G, to_crs=target_crs)
    raise AttributeError("No compatible OSMnx projection function found.")


def add_walk_time_to_edges(G):
    # Each edge gets walk_time_s = length_m / walk_speed.
    iterator = (
        G.edges(data=True, keys=True) if G.is_multigraph() else G.edges(data=True)
    )
    for tup in iterator:
        data = tup[-1]
        length_m = float(data.get("length", 0.0) or 0.0)
        data["walk_time_s"] = meters_to_seconds(length_m)
    return G


def to_undirected_safe(G):
    # OSMnx returns a MultiDiGraph; for Dijkstra we want a simple undirected
    # graph keeping the lightest valid edge per node pair.
    H = nx.Graph()
    for node, data in G.nodes(data=True):
        H.add_node(node, **data)

    def _add(u, v, data):
        wt = data.get("walk_time_s", np.nan)
        try:
            wt = float(wt)
        except (TypeError, ValueError):
            return
        if np.isnan(wt) or np.isinf(wt) or wt < 0:
            return
        new_data = data.copy()
        new_data["walk_time_s"] = wt
        if H.has_edge(u, v):
            if wt < H[u][v].get("walk_time_s", np.inf):
                H[u][v].update(new_data)
        else:
            H.add_edge(u, v, **new_data)

    if G.is_multigraph():
        for u, v, _, data in G.edges(keys=True, data=True):
            _add(u, v, data)
    else:
        for u, v, data in G.edges(data=True):
            _add(u, v, data)
    return H


def deduplicate_stops_by_mode(stops_gdf, tolerance_m=5):
  
    if stops_gdf.empty:
        return stops_gdf

    stops_metric = stops_gdf.to_crs(CRS_METRIC).reset_index(drop=True).copy()
    keep_mask = np.zeros(len(stops_metric), dtype=bool)
    before = len(stops_metric)

    for _, group in stops_metric.groupby("mode"):
        idx = group.index.values
        coords = np.column_stack([group.geometry.x.values, group.geometry.y.values])

        if len(coords) == 1:
            keep_mask[idx] = True
            continue

        tree = cKDTree(coords)
        # output_type="ndarray" needs scipy >= 1.6; older versions return a set.
        try:
            pairs = tree.query_pairs(r=tolerance_m, output_type="ndarray")
        except TypeError:
            pair_set = tree.query_pairs(r=tolerance_m)
            if pair_set:
                pairs = np.fromiter(
                    (i for pair in pair_set for i in pair),
                    dtype=np.intp, count=2 * len(pair_set),
                ).reshape(-1, 2)
            else:
                pairs = np.empty((0, 2), dtype=np.intp)

        if len(pairs) == 0:
            keep_mask[idx] = True
            continue

        # Build a sparse graph from the close-pair list, then label clusters.
        n = len(coords)
        rows = np.concatenate([pairs[:, 0], pairs[:, 1]])
        cols = np.concatenate([pairs[:, 1], pairs[:, 0]])
        data = np.ones(len(rows), dtype=np.int8)
        graph = csr_matrix((data, (rows, cols)), shape=(n, n))
        _, labels = connected_components(graph, directed=False)

        # Keep one survivor per cluster (the first occurrence).
        _, first_local = np.unique(labels, return_index=True)
        keep_mask[idx[first_local]] = True

    deduped = stops_metric[keep_mask].copy()
    print(f"  Near-duplicates removed: {before - len(deduped):,}")
    return deduped.to_crs(CRS_WGS84)


def remove_small_holes(geom, max_hole_area_m2):
    # Only fill interior rings up to a size threshold. Big legitimate gaps
    # (parks, water, non-walkable land) are kept as actual holes.
    if geom is None or geom.is_empty or max_hole_area_m2 <= 0:
        return geom
    if geom.geom_type == "Polygon":
        kept = [
            ring for ring in geom.interiors
            if Polygon(ring).area > max_hole_area_m2
        ]
        return Polygon(geom.exterior, kept)
    if geom.geom_type == "MultiPolygon":
        return MultiPolygon([
            remove_small_holes(p, max_hole_area_m2) for p in geom.geoms
        ])
    return geom


def classify_stop_count(n_stops):
    # Descriptive band for stop density (NOT the PTOL output).
    if n_stops <= 0:    return "0 - No stops"
    elif n_stops <= 2:  return "1 - Very low"
    elif n_stops <= 4:  return "2 - Low"
    elif n_stops <= 9:  return "3 - Medium"
    elif n_stops <= 19: return "4 - High"
    else:               return "5 - Very high"


def classify_ptol_modes(n_modes):
    if n_modes <= 0:   return "0 - No PT access"
    elif n_modes == 1: return "1 - Access to 1 mode"
    elif n_modes == 2: return "2 - Access to 2 modes"
    elif n_modes == 3: return "3 - Access to 3 modes"
    else:              return "4 - Access to 4 modes"


def generate_catchment_distances(G, stop_nodes, max_time_s):
    # Multi-source Dijkstra with a time cutoff. Returns {node: walk_time_s}
    # for every node reachable from at least one stop within max_time_s.
    valid_sources = {n for n in stop_nodes if n in G}
    if not valid_sources:
        return {}
    return dict(nx.multi_source_dijkstra_path_length(
        G, sources=valid_sources, cutoff=max_time_s, weight="walk_time_s"
    ))


def _partial_from_endpoint(line, frac, endpoint_xy):
    # Return the portion of `line` of normalized length `frac` starting from
    # whichever line endpoint sits closer to `endpoint_xy`. The (u,v) order
    # in NetworkX may not match the LineString's stored direction, so we
    # detect direction geometrically instead of trusting the tuple order.
    start = line.coords[0]
    end   = line.coords[-1]
    d_start = (start[0] - endpoint_xy[0])**2 + (start[1] - endpoint_xy[1])**2
    d_end   = (end[0]   - endpoint_xy[0])**2 + (end[1]   - endpoint_xy[1])**2
    if d_start <= d_end:
        return substring(line, 0.0, frac, normalized=True)
    return substring(line, 1.0 - frac, 1.0, normalized=True)


def build_isochrone_polygon(G, distances, max_time_s,
                            edge_buffer_m=25, buffer_resolution=4,
                            hole_area_m2=5_000):
   
    if not distances:
        return None

    reached = set(distances.keys())
    edge_geoms = []

    for u, v, data in G.edges(data=True):
        u_in = u in reached
        v_in = v in reached
        if not (u_in or v_in):
            continue

        # Prefer the edge's stored geometry; otherwise build a straight line.
        line = data.get("geometry")
        if line is None:
            try:
                line = LineString([
                    (G.nodes[u]["x"], G.nodes[u]["y"]),
                    (G.nodes[v]["x"], G.nodes[v]["y"]),
                ])
            except KeyError:
                continue

        if u_in and v_in:
            edge_geoms.append(line)
            continue

        edge_time = float(data.get("walk_time_s", 0.0) or 0.0)
        if edge_time <= 0:
            continue

        # Partial segment from whichever endpoint was reached.
        if u_in:
            remaining = max_time_s - distances[u]
            if remaining <= 0:
                continue
            frac = min(remaining / edge_time, 1.0)
            try:
                u_xy = (G.nodes[u]["x"], G.nodes[u]["y"])
                edge_geoms.append(_partial_from_endpoint(line, frac, u_xy))
            except Exception:
                pass
        else:
            remaining = max_time_s - distances[v]
            if remaining <= 0:
                continue
            frac = min(remaining / edge_time, 1.0)
            try:
                v_xy = (G.nodes[v]["x"], G.nodes[v]["y"])
                edge_geoms.append(_partial_from_endpoint(line, frac, v_xy))
            except Exception:
                pass

    if not edge_geoms:
        return None

    edge_series = gpd.GeoSeries(edge_geoms, crs=CRS_METRIC)
    buffered = edge_series.buffer(edge_buffer_m, resolution=buffer_resolution)
    iso = buffered.union_all() if hasattr(buffered, "union_all") else buffered.unary_union
    iso = iso.buffer(0)  # repair tiny topology artifacts

    if hole_area_m2 > 0:
        iso = remove_small_holes(iso, hole_area_m2)

    return iso


# main pipeline


print("=" * 70)
print("PTOL MILAN")
print("=" * 70)

t0 = time.time()
mode_keys = list(MODE_CONFIG.keys())


# step 1: study-area boundary 

print("\n[1/6] Loading Milan boundary...")

if not os.path.isfile(MILAN_BOUNDARY):
    raise FileNotFoundError(f"Boundary file not found: {MILAN_BOUNDARY}")

milan_3035 = read_vector_checked(MILAN_BOUNDARY, target_crs=CRS_METRIC, layer_name="Milan boundary")
milan_3035 = milan_3035[milan_3035.geometry.notna() & ~milan_3035.geometry.is_empty].copy()
if milan_3035.empty:
    raise ValueError("Milan boundary contains no valid geometries.")

milan_wgs84 = milan_3035.to_crs(CRS_WGS84)
milan_union_wgs84 = union_geometries(milan_wgs84)
print(f"  Area: {milan_3035.area.sum() / 1e6:.1f} km²")


#step 2: H3 grid covering the boundary 
print(f"\n[2/6] Building H3 grid, resolution {H3_RESOLUTION}...")

hex_gdf = build_h3_grid_from_polygon(milan_union_wgs84, H3_RESOLUTION)
if hex_gdf.empty:
    raise ValueError("No H3 cells generated.")

# Clip hexes at the boundary so border cells are not counted whole.
hex_gdf = gpd.clip(hex_gdf, milan_wgs84, keep_geom_type=True)
if hex_gdf.empty:
    raise ValueError("No H3 cells after clipping.")

hex_gdf = hex_gdf.reset_index(drop=True)
hex_gdf["area_km2"] = (hex_gdf.to_crs(CRS_METRIC).area / 1e6).round(4)
print(f"  H3 cells: {len(hex_gdf):,}")


# step 3: PT stops, one mode at a time 

print("\n[3/6] Loading public transport stops...")

all_stops = []
for mode, filepath in STOPS_FILES.items():
    if not os.path.isfile(filepath):
        print(f"  {mode}: file not found, skipped.")
        continue
    gdf = read_vector_checked(filepath, target_crs=CRS_WGS84, layer_name=f"{mode} stops")
    gdf = ensure_point_geometry(gdf)
    if gdf.empty:
        print(f"  {mode}: empty, skipped.")
        continue
    gdf = gdf[gdf.geometry.within(milan_union_wgs84)].copy()
    if gdf.empty:
        print(f"  {mode}: no points inside Milan, skipped.")
        continue
    gdf["mode"] = mode
    all_stops.append(gdf[["geometry", "mode"]].copy())
    print(f"  {mode}: {len(gdf):,} points")

if not all_stops:
    raise FileNotFoundError("No valid stop layers loaded.")

stops_gdf = gpd.GeoDataFrame(pd.concat(all_stops, ignore_index=True), crs=CRS_WGS84)

# Single dedup pass: density and isochrones use the same stop universe.
if DEDUPLICATE_STOPS:
    stops_gdf = deduplicate_stops_by_mode(stops_gdf, DEDUP_TOLERANCE_M)

loaded_modes = sorted(stops_gdf["mode"].unique())
print(f"  Total stops after dedup: {len(stops_gdf):,}")
print(f"  Modes loaded: {', '.join(loaded_modes)}")

missing_modes = sorted(set(mode_keys) - set(loaded_modes))
if missing_modes:
    print(f"  WARNING: Missing modes: {', '.join(missing_modes)}")


# step 4: stop counts and density per hex 
#
# This is a descriptive stop-density measure, NOT the PTOL output.
# Useful as context but does not measure walking accessibility.

print("\n[4/6] Stop counts and density per H3 cell...")

stops_in_hex = gpd.sjoin(
    stops_gdf[["geometry", "mode"]],
    hex_gdf[["h3_index", "geometry"]],
    how="inner", predicate="within",
)
counts_by_mode = stops_in_hex.groupby(["h3_index", "mode"]).size().unstack(fill_value=0)

count_cols, density_cols = [], []
for mode in mode_keys:
    n_col = f"n_{mode.lower()}"
    d_col = f"density_{mode.lower()}"
    count_cols.append(n_col)
    density_cols.append(d_col)

    if mode in counts_by_mode.columns:
        hex_gdf[n_col] = hex_gdf["h3_index"].map(counts_by_mode[mode]).fillna(0).astype(int)
    else:
        hex_gdf[n_col] = 0

    hex_gdf[d_col] = np.where(
        hex_gdf["area_km2"] > 0,
        hex_gdf[n_col] / hex_gdf["area_km2"], 0.0
    ).round(1)

hex_gdf["n_stops_total"] = hex_gdf[count_cols].sum(axis=1)
hex_gdf["density_total_km2"] = np.where(
    hex_gdf["area_km2"] > 0,
    hex_gdf["n_stops_total"] / hex_gdf["area_km2"], 0.0
).round(1)
hex_gdf["stop_count_band"] = hex_gdf["n_stops_total"].apply(classify_stop_count)

print(f"  Stops assigned to cells: {len(stops_in_hex):,}")
print(f"  Cells with >=1 stop: {(hex_gdf['n_stops_total'] > 0).sum():,}")
for mode, col in zip(mode_keys, count_cols):
    print(f"    {mode}: {hex_gdf[col].sum():,} stops")


# step 5: PTOL = walking accessibility to each PT mode 

print("\n[5/6] Calculating isochrone-based PTOL...")

if not os.path.isfile(WALK_NETWORK_FILE):
    raise FileNotFoundError(
        f"Walk network not found: {WALK_NETWORK_FILE}\n"
        "Generate the pedestrian network GraphML before running."
    )

# Load OSM pedestrian graph, project to meters, add walk_time per edge,
# then collapse to a simple undirected graph for Dijkstra.
G_walk = load_graphml_compat(WALK_NETWORK_FILE)
G_walk_projected = project_graph_compat(G_walk, CRS_METRIC)
G_walk_projected = add_walk_time_to_edges(G_walk_projected)
G_undirected = to_undirected_safe(G_walk_projected)
graph_node_set = set(G_undirected.nodes)

print(f"  Walk graph: {G_undirected.number_of_nodes():,} nodes, "
      f"{G_undirected.number_of_edges():,} edges")

# Snap each stop to its nearest pedestrian-network node.
stops_metric = stops_gdf.to_crs(CRS_METRIC).copy()
stop_nodes, stop_distances = ox.distance.nearest_nodes(
    G_walk_projected,
    X=stops_metric.geometry.x.values,
    Y=stops_metric.geometry.y.values,
    return_dist=True,
)
stops_metric["nn_node"] = stop_nodes
stops_metric["snap_dist_m"] = stop_distances

# Drop stops whose nearest node was lost going from multidigraph -> Graph.
before_snap = len(stops_metric)
stops_metric = stops_metric[
    stops_metric["nn_node"].isin(graph_node_set)
    & stops_metric["snap_dist_m"].notna()
].copy()
print(f"  Stops dropped (invalid node): {before_snap - len(stops_metric):,}")

# Drop stops snapped suspiciously far (likely data errors).
before_far = len(stops_metric)
stops_metric = stops_metric[stops_metric["snap_dist_m"] <= MAX_SNAP_DISTANCE_M].copy()
print(f"  Stops dropped (snap > {MAX_SNAP_DISTANCE_M} m): "
      f"{before_far - len(stops_metric):,}")

if stops_metric.empty:
    raise ValueError(
        "No valid stops after snapping/filtering. Check the pedestrian "
        "network, stops, CRS, or MAX_SNAP_DISTANCE_M."
    )

print("  Snapping diagnostics (after filtering):")
print(f"    Mean:   {stops_metric['snap_dist_m'].mean():.1f} m")
print(f"    Median: {stops_metric['snap_dist_m'].median():.1f} m")
print(f"    P95:    {stops_metric['snap_dist_m'].quantile(0.95):.1f} m")
print(f"    Max:    {stops_metric['snap_dist_m'].max():.1f} m")
print(f"  Stops kept for isochrones: {len(stops_metric):,}")

# One isochrone polygon per mode (union of all stop catchments).
iso_records = []
for mode in mode_keys:
    mode_stops = stops_metric[stops_metric["mode"] == mode]
    if mode_stops.empty:
        print(f"\n  {mode}: no valid stops after snapping.")
        continue

    max_time_s = MODE_CONFIG[mode]["max_time_min"] * 60.0
    stop_node_set = set(mode_stops["nn_node"].unique())

    print(f"\n  {mode}: {len(stop_node_set):,} unique stop nodes, "
          f"threshold {MODE_CONFIG[mode]['max_time_min']:.0f} min")

    distances = generate_catchment_distances(G_undirected, stop_node_set, max_time_s)
    if not distances:
        print(f"    Empty catchment.")
        continue
    print(f"    Reachable nodes: {len(distances):,}")

    iso_geom = build_isochrone_polygon(
        G_undirected, distances, max_time_s,
        edge_buffer_m=EDGE_BUFFER_M,
        buffer_resolution=EDGE_BUFFER_RESOLUTION,
        hole_area_m2=SMALL_HOLE_AREA_M2,
    )
    if iso_geom is None or iso_geom.is_empty:
        print(f"    Could not generate isochrone polygon.")
        continue

    iso_records.append({
        "mode": mode,
        "threshold_min": MODE_CONFIG[mode]["max_time_min"],
        "geometry": iso_geom,
    })
    print(f"    Isochrone polygon generated.")

if not iso_records:
    raise ValueError("No isochrone polygons were generated.")

iso_gdf = gpd.GeoDataFrame(iso_records, crs=CRS_METRIC)
iso_modes = sorted(iso_gdf["mode"].unique())
missing_iso = sorted(set(mode_keys) - set(iso_modes))
if missing_iso:
    print(f"\n  WARNING: No isochrone for: {', '.join(missing_iso)}")

# Cross each isochrone with the H3 grid. A hex is "accessible to mode X"
# if a configurable fraction of its area sits inside mode X's isochrone.
print(f"\n  Intersecting isochrones with H3 grid "
      f"(min overlap = {MIN_HEX_OVERLAP_FRACTION:.0%} of hex area)...")

hex_metric_geom = hex_gdf.to_crs(CRS_METRIC).geometry
hex_area = hex_metric_geom.area

iso_cols = []
for mode in mode_keys:
    col = f"iso_access_{mode.lower()}"
    iso_cols.append(col)
    hex_gdf[col] = 0

for _, row in iso_gdf.iterrows():
    col = f"iso_access_{row['mode'].lower()}"
    if MIN_HEX_OVERLAP_FRACTION <= 0:
        accessible = hex_metric_geom.intersects(row["geometry"])
    else:
        inter_area = hex_metric_geom.intersection(row["geometry"]).area
        accessible = (inter_area / hex_area) >= MIN_HEX_OVERLAP_FRACTION
    hex_gdf[col] = accessible.astype(int).values

# Final PTOL score: sum of accessible modes (0..4).
hex_gdf["ptol_modes_iso"] = hex_gdf[iso_cols].sum(axis=1)
hex_gdf["ptol_iso_band"] = hex_gdf["ptol_modes_iso"].apply(classify_ptol_modes)

print(f"\n  PTOL summary:")
print(f"    Mean accessible modes: {hex_gdf['ptol_modes_iso'].mean():.2f}")
print(f"    Max accessible modes:  {hex_gdf['ptol_modes_iso'].max():.0f}")
for band, count in hex_gdf["ptol_iso_band"].value_counts().sort_index().items():
    pct = count / len(hex_gdf) * 100
    print(f"    {band:24s} {count:7,} cells  {pct:5.1f}%")


# step 6: export 
print("\n[6/6] Exporting...")

export_cols = (
    ["h3_index", "clat", "clon", "area_km2"]
    + count_cols + ["n_stops_total"]
    + density_cols + ["density_total_km2"]
    + ["stop_count_band"]
    + iso_cols
    + ["ptol_modes_iso", "ptol_iso_band"]
    + ["geometry"]
)
hex_export = hex_gdf[[c for c in export_cols if c in hex_gdf.columns]].copy()

gpkg_path = os.path.join(OUTPUT_DIR, "PTOL_MILAN_H3R9_final.gpkg")
hex_export.to_crs(CRS_METRIC).to_file(gpkg_path, driver="GPKG")
print(f"  GeoPackage: {gpkg_path}")

csv_path = os.path.join(OUTPUT_DIR, "PTOL_MILAN_H3R9_FINAL.csv")
hex_export.drop(columns="geometry").to_csv(csv_path, index=False, na_rep="")
print(f"  CSV: {csv_path}")

iso_path = os.path.join(OUTPUT_DIR, "PTOL_MILAN_H3R9_isochrones.gpkg")
iso_gdf.to_crs(CRS_METRIC).to_file(iso_path, driver="GPKG")
print(f"  Isochrones: {iso_path}")

elapsed = time.time() - t0
print("\n" + "=" * 70)
print("PTOL COMPLETE")
print(f"  Runtime: {elapsed / 60:.1f} min")
print(f"  H3 cells: {len(hex_gdf):,}")
print(f"  Stops counted: {int(hex_gdf['n_stops_total'].sum()):,}")
print(f"  Modes: {', '.join(loaded_modes)}")
print(f"  Output: {OUTPUT_DIR}")
print("=" * 70)
