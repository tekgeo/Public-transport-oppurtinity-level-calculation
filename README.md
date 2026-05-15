# Public-transport-oppurtinity-level-calculation
Public Transport Opportunity Level (PTOL) for Milan on an H3 hex grid, computed from walking isochrones over the OSM pedestrian network.

`PTOL_CALCULATION.py` computes a Public Transport Opportunity Level
(PTOL) for the Città Metropolitana di Milano on an H3 hexagonal grid.


## What PTOL measures

For each hex cell, PTOL counts **how many public transport modes are
accessible on foot within a mode-specific walking time**:

| Mode    | Walking-time threshold |
|---------|------------------------|
| Bus     | 8 min                  |
| Tram    | 8 min                  |
| Metro   | 12 min                 |
| Railway | 12 min                 |

A walking speed of 80 m/min (4.8 km/h) is used throughout, matching the
PTAL convention. A hex is "accessible to mode X" when a configurable
fraction of its area falls inside the walking isochrone generated from
mode X's stops on the OSM pedestrian network.

The PTOL value of a hex is the **sum of accessible modes**, so it takes
integer values from 0 (no PT accessible on foot) to 4 (access to all
four modes within their thresholds).

PTOL is intentionally simpler than PTAL: it does not use GTFS, service
frequencies, or waiting times. It measures *opportunity to reach a stop
on foot*, not service quality.

## Pipeline

```
[1] Boundary       → load Città Metropolitana polygon
[2] H3 grid        → generate cells (resolution 9 ≈ 0.1 km²), clip to boundary
[3] PT stops       → load one point layer per mode, deduplicate per mode
[4] Stop density   → describe offer per hex (counts, density/km²); NOT PTOL
[5] PTOL           → walking isochrones per mode → intersect with hex grid
[6] Export         → GeoPackage, CSV, and the isochrone polygons themselves
```

### Step-by-step

**Boundary.** A shapefile of the study area in EPSG:4326 (or any defined
CRS — it gets reprojected). The H3 grid is generated to cover the union
of its polygons and then clipped at the border.

**H3 grid.** Resolution 9 produces ~0.1 km² hexes, roughly 360 m across.
Rome used resolution 10 (≈ 0.015 km²) for finer detail; res 9 trades
detail for runtime. Border hexes are clipped to the boundary so their
area reflects the in-study-area portion.

**Stops.** One shapefile per mode in EPSG:4326. Polygon inputs are
collapsed to centroids. Stops outside the boundary are dropped.
Near-duplicate stops within a configurable tolerance (default 5 m) are
clustered using a KD-tree + connected-components and reduced to one
survivor per cluster. The same deduplicated stop set feeds both step 4
(density) and step 5 (isochrones).

**Stop density (step 4).** Spatial join of stops into hexes, with counts
and density per km² per mode. This is a descriptive metric of *offer*,
not of *opportunity*. The PTOL output lives in step 5.

**Isochrones (step 5).** The OSM pedestrian graph (downloaded once and
saved as `.graphml`) is loaded, projected to EPSG:3035, and converted to
a simple undirected graph keeping the lightest edge per node pair.
Walking time per edge = `length / 80 m·min⁻¹`.

Stops are snapped to the nearest pedestrian-network node; those snapped
more than `MAX_SNAP_DISTANCE_M` (default 100 m) are dropped as data
errors.

For each mode, a multi-source Dijkstra with a time cutoff returns every
network node reachable from any stop within the mode's threshold.

The isochrone polygon is built from those reachable nodes. Edges with
both endpoints reached are added in full; edges with one endpoint
reached contribute a partial segment up to the remaining time budget
(`shapely.ops.substring`). This is the standard correction for
edge-cropping bias at the isochrone boundary. Segments are buffered by
`EDGE_BUFFER_M` and unioned. Small interior holes (default ≤ 5 000 m² —
roughly inner courtyards) are filled; large holes such as parks, water,
or the Cimitero Monumentale remain as actual holes in the polygon.

**Hex × isochrone.** A hex is flagged as accessible to a mode if at
least `MIN_HEX_OVERLAP_FRACTION` (default 25 %) of its area falls
inside the mode's isochrone. PTOL is the count of accessible modes
across all four PT layers.

## Inputs

```
shapes/
  Citta Metropolitana di Milano.shp      # study-area boundary

Data/
  bus_stops.shp                          # one point layer per PT mode
  tram_stops.shp
  metro_stations.shp
  railway_stations.shp
  citta_metropolitana_milano_walk_network.graphml
                                         # OSMnx pedestrian graph
                                         # of the study area
```

All vector inputs must have a defined CRS. PT-stop layers should
contain only stops within the study area (or a slight buffer of it);
they are filtered to the boundary at load time.

The walking-network graph must cover the boundary *with a buffer*
(≥ 1 km recommended). Generating it with only the boundary itself
prevents paths from briefly leaving and re-entering Milan, which
truncates isochrones at the border.

To download the pedestrian graph once:

```python
import osmnx as ox
import geopandas as gpd

boundary = gpd.read_file("shapes/Citta Metropolitana di Milano.shp").to_crs(4326)
poly = boundary.geometry.union_all().buffer(0.02)   # ~2 km buffer in degrees

G = ox.graph_from_polygon(poly, network_type="walk", simplify=True)
ox.save_graphml(G, "Data/citta_metropolitana_milano_walk_network.graphml")
```

## Outputs

| File                                       | Content |
|--------------------------------------------|---------|
| `PTOL_MILAN_H3R9_final.gpkg`               | H3 grid in EPSG:3035 with all attributes |
| `PTOL_MILAN_H3R9_FINAL.csv`                | Same attributes, no geometry |
| `PTOL_MILAN_H3R9_isochrones.gpkg`          | One polygon per mode (the input to step 5's intersection) |

### Key columns in the hex output

| Column                  | Meaning |
|-------------------------|---------|
| `h3_index`              | H3 cell ID |
| `clat`, `clon`          | Cell centroid (lat/lon) |
| `area_km2`              | Cell area after clipping to boundary |
| `n_<mode>`              | Count of stops of `<mode>` falling inside the cell |
| `density_<mode>`        | `n_<mode>` per km² |
| `n_stops_total`         | All modes combined |
| `density_total_km2`     | All modes combined per km² |
| `stop_count_band`       | Descriptive band on total stop count |
| `iso_access_<mode>`     | 1 if the cell's overlap with mode's isochrone meets the threshold, else 0 |
| **`ptol_modes_iso`**    | **PTOL score: sum of `iso_access_*`. Values 0–4.** |
| `ptol_iso_band`         | Text label for `ptol_modes_iso` |

## Configuration

All knobs are at the top of the script. The ones most likely to be
adjusted:

| Variable                   | Default | What it does |
|----------------------------|---------|--------------|
| `H3_RESOLUTION`            | 9       | Cell size. Use 10 for ~6× more cells and finer detail. |
| `MODE_CONFIG[*]["max_time_min"]` | 8 / 12 | Walking thresholds per mode. |
| `EDGE_BUFFER_M`            | 25      | Half-width of the walkable corridor around each reached edge. |
| `EDGE_BUFFER_RESOLUTION`   | 4       | Buffer segments. 4 = octagon, much faster than the default circle. |
| `SMALL_HOLE_AREA_M2`       | 5 000   | Maximum hole size that gets filled. Larger gaps (parks, water) stay open. |
| `DEDUP_TOLERANCE_M`        | 5       | Stops within this distance get clustered into one. Raise to ~50–100 if your stop layer represents station entrances/platforms rather than logical stops. |
| `MAX_SNAP_DISTANCE_M`      | 100     | Stops farther than this from any pedestrian node are dropped. |
| `MIN_HEX_OVERLAP_FRACTION` | 0.25    | 0 = any intersection makes a hex accessible; 0.5 = majority. |

## Methodological notes and caveats

- **No frequency weighting.** A hex with metro access (3-min headways) and
  one with railway access (60-min headways) both score 1 on that mode.
  PTOL is *opportunity to reach*, not *quality of service*. Adding GTFS
  frequency weighting is a natural next step but moves the metric toward
  PTAL.

- **Long-edge approximation.** When both endpoints of an edge are
  reached, the entire edge is included even though a midpoint of a very
  long edge might in fact be > threshold from any stop. With typical
  urban edges of ~50 m and thresholds of 8–12 min this approximation is
  negligible.

- **Stop-layer hygiene.** If a stop layer represents station *entrances*
  rather than logical stations, a single metro station can appear as 8+
  points within tens of meters. The 5-m deduplication does not cluster
  these. Either filter the source layer (`public_transport=station`
  rather than `subway_entrance` in OSM) or raise `DEDUP_TOLERANCE_M` to
  ~80–120 m for metro/railway.

- **Border effects.** The walk graph must extend past the study-area
  boundary or isochrones get clipped artificially at the edge.

- **Hex × isochrone threshold.** `MIN_HEX_OVERLAP_FRACTION = 0.25` is a
  pragmatic middle ground. Any-intersection is too lax (a hex touching
  the isochrone by 1 m² counts); centroid-only is too strict in border
  cells; 25–50 % area overlap is a reasonable default for H3 res 9.

## Dependencies

```bash
pip install "osmnx>=2" "networkx>=3" "geopandas>=0.14" "h3>=4" \
            shapely numpy pandas pyogrio scipy
```

`osmnx` v1 also works (the compatibility shims `load_graphml_compat`
and `project_graph_compat` cover both APIs). `scipy >= 1.6` is
recommended for the KD-tree `output_type="ndarray"` shortcut, but the
script falls back to the set-based API automatically.

## Running

Edit the paths in the configuration block to point at your data, then:

```bash
python PTOL_CALCULATION.py
```

Typical runtime for the Città Metropolitana di Milano at H3 resolution
9 with ~10–15 k stops is in the order of minutes on a laptop. Most of
the time is in the per-mode buffer-and-union step; lower
`EDGE_BUFFER_RESOLUTION` if it becomes a bottleneck.
