# project_fem2geo

Small CLI tool to project geographic data (lon/lat/depth) and simple VTK meshes into a FEM-style ENU coordinate system.

It supports:
- **points**: CSV tables (lon/lat + depth) → projected **x/y/z** and optional VTP export
- **mesh**: VTK meshes (VTP/VTU/VTK) in ENU cartesian → projected FEM coordinates
- optional **slicing** by lon/lat/depth (deg, deg, km down)
- optional **anchoring** to match a real-world location to a FEM coordinate origin
- optional **rotation** around the output anchor


## Installation (venv + pip)

```bash
python -m venv venv
source venv/bin/activate
pip install .
```

If you are developing locally:

```bash
pip install -e .
```


## Running the examples

Examples live in the `examples/` folder. Each example contains `config.yml` files showing a workflow.

Run them from inside each example folder:

```bash
cd examples/a_points2crs
project_fem2geo points config_interface.yml
project_fem2geo points config_intraslab.yml

cd ../b_mesh2crs
project_fem2geo mesh config_top.yml
project_fem2geo mesh config_bot.yml

cd ../c_points2fem
project_fem2geo points config.yml

cd ../d_mesh2fem
project_fem2geo mesh config.yml
```


## What each example means

### `a_points2crs`
**Goal:** Convert a CSV catalog of points from lon/lat/depth into a projected CRS (e.g., UTM), without anchoring.

- Input: CSV with columns like `[longitude, latitude, depth]`
- Output: usually a `*.vtp` point cloud (for ParaView) or `*.csv`
- Slicing: `lon_range`, `lat_range`, `depth_range_km` reduce the dataset before writing

Use this when you just want the points in the same projected CRS as your other data (slab, model grids, etc).


### `b_mesh2crs`
**Goal:** Transform an existing ENU mesh into a target CRS/units, no anchor.

- Input: mesh in ENU cartesian coordinates
- You specify `input.mesh_units` (`m` or `km`)
- The mesh is written to the output path

Use this when the mesh is already in the same global CRS but you need consistent units or a clean subset via slicing.

### `c_points2fem`
**Goal:** Convert lon/lat/depth points into FEM coordinates using an **anchor**.

You provide:
- `anchor.src_geo = [lon, lat, depth_km]` (depth positive down)
- `anchor.dst = [x, y, z]` in FEM coordinates

The tool computes a translation (and optional rotation) so that the anchor maps exactly to your FEM reference point.

Use this when the FEM model uses its own local origin (0,0,0) and you want geographic data to land in the FEM frame.



### `d_mesh2fem`
**Goal:** Transform an ENU mesh into FEM coordinates using an anchor.

Same anchoring idea as for points, but applied to mesh vertices.
Slicing (lon/lat/depth) is applied before output.

Use this when you have a surface mesh (e.g., slab interface) in a projected CRS and you need it aligned into your FEM local coordinate system.


## Output formats

### Points
- `output.format: csv` writes the input table + extra `x,y,z` columns
- `output.format: vtp` writes a ParaView-friendly point cloud (`.vtp`)

### Mesh
- Output format is inferred from the output filename extension.
- If slicing causes the mesh type to change (e.g., polydata → unstructured), the tool may switch the extension and print a warning.


## Notes on slicing

Slicing is always described in:
- `lon_range`: degrees (EPSG:4326)
- `lat_range`: degrees (EPSG:4326)
- `depth_range_km`: km, positive down

This is consistent across both points and meshes.

For meshes, depth is derived from ENU z assuming:
- ENU z is “up”
- depth_km = `-z_m / 1000`

