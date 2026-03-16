from pathlib import Path

import numpy as np
import pandas as pd
import pyvista as pv
import yaml
from pyproj import CRS, Transformer

UNIT_TO_M = {"m": 1.0, "km": 1000.0}


def unit_to_m(unit):
    """
    Return meters-per-unit for 'm' or 'km'.
    """
    u = str(unit).strip().lower()
    if u not in UNIT_TO_M:
        raise ValueError("unit must be 'm' or 'km'.")
    return UNIT_TO_M[u]


def check_updown(value, label):
    """
    Validate 'up' / 'down'.
    """
    v = str(value).strip().lower()
    if v not in ("up", "down"):
        raise ValueError(label + " must be 'up' or 'down'.")
    return v


def read_config(path):
    """
    Read YAML config and resolve input/output file paths relative to the config file.
    """
    path = Path(path)
    cfg = yaml.safe_load(path.read_text())
    base = path.parent

    if "crs" not in cfg or "src" not in cfg["crs"] or "dst" not in cfg["crs"]:
        raise ValueError("crs must have 'src' and 'dst'.")

    if "input" not in cfg or "file" not in cfg["input"]:
        raise ValueError("input.file is required.")
    if "output" not in cfg or "file" not in cfg["output"]:
        raise ValueError("output.file is required.")

    inp = cfg["input"]
    out = cfg["output"]

    inp_file = Path(inp["file"])
    out_file = Path(out["file"])
    inp["file"] = str(inp_file if inp_file.is_absolute() else (base / inp_file))
    out["file"] = str(out_file if out_file.is_absolute() else (base / out_file))

    if "xy_units" not in cfg:
        cfg["xy_units"] = {"src": "deg", "dst": "m"}
    else:
        if "src" not in cfg["xy_units"] or "dst" not in cfg["xy_units"]:
            raise ValueError("xy_units must have 'src' and 'dst'.")

    return cfg


def to_lonlat(cfg, x, y):
    """
    Convert arrays x,y from crs.src to EPSG:4326 lon/lat (degrees).
    """
    crs_src = CRS.from_user_input(cfg["crs"]["src"])
    crs_ll = CRS.from_epsg(4326)
    tfm = Transformer.from_crs(crs_src, crs_ll, always_xy=True)
    lon, lat = tfm.transform(np.asarray(x, dtype=float), np.asarray(y, dtype=float))
    return np.asarray(lon, dtype=float), np.asarray(lat, dtype=float)


def slice_mask_lonlatdepth(lon, lat, depth_km, lon_range, lat_range, depth_range_km):
    """
    Build a boolean mask from optional lon/lat/depth ranges.
    """
    m = np.ones(len(lon), dtype=bool)

    if lon_range is not None:
        lo0, lo1 = float(lon_range[0]), float(lon_range[1])
        m &= (lon >= lo0) & (lon <= lo1)

    if lat_range is not None:
        la0, la1 = float(lat_range[0]), float(lat_range[1])
        m &= (lat >= la0) & (lat <= la1)

    if depth_range_km is not None:
        d0, d1 = float(depth_range_km[0]), float(depth_range_km[1])
        m &= (depth_km >= d0) & (depth_km <= d1)

    return m


def depth_km_from_point_z(z_vals, z_units, z_positive):
    """
    Convert point Z values to depth (km, positive down).
    """
    z = np.asarray(z_vals, dtype=float) * unit_to_m(z_units)
    zp = check_updown(z_positive, "depth.src_positive")
    depth_m = z if zp == "down" else -z
    return depth_m / 1000.0


def convert_point_z(z_vals, src_units, src_positive, dst_units, dst_positive):
    """
    Convert point Z values from (units, sign) to (units, sign).
    """
    z = np.asarray(z_vals, dtype=float) * unit_to_m(src_units)
    sp = check_updown(src_positive, "depth.src_positive")
    dp = check_updown(dst_positive, "depth.dst_positive")

    if sp != dp:
        z = -z

    return z / unit_to_m(dst_units)


def compute_anchor_shift(cfg):
    """
    Compute (dx, dy, dz) in output units, and (x0, y0, az_deg) for optional rotation.

    Anchor input is lon/lat/depth_km (depth positive down).
    Anchor output is x/y/z in output XY units and output z convention.
    """
    anchor = cfg.get("anchor", None)
    rotation = cfg.get("rotation", {})

    if anchor is None:
        return 0.0, 0.0, 0.0, 0.0, 0.0, None

    if "src_geo" not in anchor or "dst" not in anchor:
        raise ValueError("anchor must have 'src_geo' and 'dst'.")

    src_geo = anchor["src_geo"]
    dst = anchor["dst"]
    if (not isinstance(src_geo, (list, tuple))) or len(src_geo) != 3:
        raise ValueError("anchor.src_geo must be [lon, lat, depth_km].")
    if (not isinstance(dst, (list, tuple))) or len(dst) != 3:
        raise ValueError("anchor.dst must be [x, y, z].")

    lon0 = float(src_geo[0])
    lat0 = float(src_geo[1])
    depth0_km = float(src_geo[2])

    x0 = float(dst[0])
    y0 = float(dst[1])
    z0 = float(dst[2])

    crs_dst = CRS.from_user_input(cfg["crs"]["dst"])
    tfm = Transformer.from_crs(CRS.from_epsg(4326), crs_dst, always_xy=True)
    Xm, Ym = tfm.transform(lon0, lat0)

    dst_xy_units = str(cfg["xy_units"]["dst"]).strip().lower()
    if dst_xy_units not in ("m", "km"):
        raise ValueError("xy_units.dst must be 'm' or 'km'.")

    ax = float(Xm) / unit_to_m(dst_xy_units)
    ay = float(Ym) / unit_to_m(dst_xy_units)

    dp = cfg.get("depth", {})
    dst_pos = check_updown(dp.get("dst_positive", "up"), "depth.dst_positive")

    depth_m = depth0_km * 1000.0
    depth_dst = depth_m / unit_to_m(dst_xy_units)
    z_anchor = -depth_dst if dst_pos == "up" else depth_dst

    dx = x0 - ax
    dy = y0 - ay
    dz = z0 - z_anchor

    az = rotation.get("azimuth_deg", None)
    az_deg = None if az is None else float(az)

    return dx, dy, dz, x0, y0, az_deg


def rotate_xy(x, y, x0, y0, az_deg):
    """
    Rotate x,y around x0,y0 by azimuth degrees (CCW).
    """
    t = np.deg2rad(float(az_deg))
    c, s = np.cos(t), np.sin(t)
    xy = np.vstack([x - x0, y - y0])
    xy2 = np.array([[c, -s], [s, c]], dtype=float) @ xy
    return xy2[0] + x0, xy2[1] + y0


def transform_xy_points(cfg, x_vals, y_vals):
    """
    Transform point XY from crs.src to crs.dst and scale to xy_units.dst.
    """
    crs_src = CRS.from_user_input(cfg["crs"]["src"])
    crs_dst = CRS.from_user_input(cfg["crs"]["dst"])
    tfm = Transformer.from_crs(crs_src, crs_dst, always_xy=True)

    src_xy_units = str(cfg["xy_units"]["src"]).strip().lower()
    dst_xy_units = str(cfg["xy_units"]["dst"]).strip().lower()
    if dst_xy_units not in ("m", "km"):
        raise ValueError("xy_units.dst must be 'm' or 'km'.")

    x = np.asarray(x_vals, dtype=float)
    y = np.asarray(y_vals, dtype=float)

    if crs_src.is_geographic:
        if src_xy_units != "deg":
            raise ValueError("xy_units.src must be 'deg' when crs.src is geographic.")
        X, Y = tfm.transform(x, y)
    else:
        if src_xy_units not in ("m", "km"):
            raise ValueError("xy_units.src must be 'm' or 'km' when crs.src is projected.")
        X, Y = tfm.transform(x * unit_to_m(src_xy_units), y * unit_to_m(src_xy_units))

    X = np.asarray(X, dtype=float) / unit_to_m(dst_xy_units)
    Y = np.asarray(Y, dtype=float) / unit_to_m(dst_xy_units)
    return X, Y


def project_points(cfg):
    """
    Transform a CSV point table. Adds x,y,z columns.
    """
    inp = cfg["input"]

    cols = inp.get("columns", None)
    if (not isinstance(cols, (list, tuple))) or len(cols) != 3:
        raise ValueError("input.columns must be [x_col, y_col, z_col].")
    x_col, y_col, z_col = cols

    df = pd.read_csv(inp["file"])

    lon_range = inp.get("lon_range", None)
    lat_range = inp.get("lat_range", None)
    depth_range_km = inp.get("depth_range_km", None)

    dp = cfg.get("depth", {})
    z_units = dp.get("src_units", "km")
    z_pos = dp.get("src_positive", "down")

    if lon_range is not None or lat_range is not None or depth_range_km is not None:
        x_raw = df[x_col].to_numpy(float)
        y_raw = df[y_col].to_numpy(float)

        if CRS.from_user_input(cfg["crs"]["src"]).is_geographic:
            lon = x_raw
            lat = y_raw
        else:
            src_xy_units = str(cfg["xy_units"]["src"]).strip().lower()
            if src_xy_units not in ("m", "km"):
                raise ValueError("xy_units.src must be 'm' or 'km' when slicing projected points.")
            lon, lat = to_lonlat(cfg, x_raw * unit_to_m(src_xy_units), y_raw * unit_to_m(src_xy_units))

        depth_km = depth_km_from_point_z(df[z_col].to_numpy(float), z_units, z_pos)
        m = slice_mask_lonlatdepth(lon, lat, depth_km, lon_range, lat_range, depth_range_km)
        df = df[m]
        if df.empty:
            raise ValueError("No points left after slicing.")

    X, Y = transform_xy_points(cfg, df[x_col].to_numpy(float), df[y_col].to_numpy(float))

    dst_units = str(cfg["xy_units"]["dst"]).strip().lower()
    dst_pos = dp.get("dst_positive", "up")
    z_out_units = dp.get("dst_units", dst_units)

    z = convert_point_z(df[z_col].to_numpy(float), z_units, z_pos, z_out_units, dst_pos)

    dx, dy, dz, x0, y0, az_deg = compute_anchor_shift(cfg)
    X = X + dx
    Y = Y + dy
    z = z + dz

    if az_deg is not None:
        X, Y = rotate_xy(X, Y, x0, y0, az_deg)

    out_df = df.copy()
    out_df["x"] = X
    out_df["y"] = Y
    out_df["z"] = z
    return out_df


def mesh_depth_km_from_enu_z(z_vals, mesh_units):
    """
    Mesh Z is ENU 'up'. Convert to depth km (positive down).
    """
    z_m = np.asarray(z_vals, dtype=float) * unit_to_m(mesh_units)
    return (-z_m) / 1000.0


def project_mesh(cfg):
    """
    Transform a mesh (vtp/vtu/vtk). Mesh is ENU; only mesh_units matters.
    Slicing is applied in lon/lat/depth (deg, deg, km down).
    """
    inp = cfg["input"]
    out = cfg["output"]

    mesh = pv.read(inp["file"])
    pts = np.asarray(mesh.points, dtype=float)
    if pts.size == 0:
        raise ValueError("Empty mesh.")

    mesh_units = inp.get("mesh_units", "m")
    mesh_units = str(mesh_units).strip().lower()
    if mesh_units not in ("m", "km"):
        raise ValueError("input.mesh_units must be 'm' or 'km'.")

    lon_range = inp.get("lon_range", None)
    lat_range = inp.get("lat_range", None)
    depth_range_km = inp.get("depth_range_km", None)

    if lon_range is not None or lat_range is not None or depth_range_km is not None:
        crs_src = CRS.from_user_input(cfg["crs"]["src"])
        to_ll = Transformer.from_crs(crs_src, CRS.from_epsg(4326), always_xy=True)

        Xm = pts[:, 0] * unit_to_m(mesh_units)
        Ym = pts[:, 1] * unit_to_m(mesh_units)
        lon, lat = to_ll.transform(Xm, Ym)
        lon = np.asarray(lon, dtype=float)
        lat = np.asarray(lat, dtype=float)

        depth_km = mesh_depth_km_from_enu_z(pts[:, 2], mesh_units)
        m = slice_mask_lonlatdepth(lon, lat, depth_km, lon_range, lat_range, depth_range_km)

        mesh = mesh.extract_points(m, adjacent_cells=True, include_cells=True)
        pts = np.asarray(mesh.points, dtype=float)
        if pts.size == 0:
            raise ValueError("No mesh left after slicing.")

    X, Y = transform_xy_points(cfg, pts[:, 0], pts[:, 1])

    dst_xy_units = str(cfg["xy_units"]["dst"]).strip().lower()
    z = (pts[:, 2].astype(float) * unit_to_m(mesh_units)) / unit_to_m(dst_xy_units)

    dx, dy, dz, x0, y0, az_deg = compute_anchor_shift(cfg)
    X = X + dx
    Y = Y + dy
    z = z + dz

    if az_deg is not None:
        X, Y = rotate_xy(X, Y, x0, y0, az_deg)

    mesh.points = np.c_[X, Y, z]

    out_path = Path(out["file"])
    ext = out_path.suffix.lower()
    is_poly = isinstance(mesh, pv.PolyData)

    if is_poly and ext == ".vtu":
        new_path = out_path.with_suffix(".vtp")
        print(f"warning: output is PolyData; writing {new_path.name} instead of {out_path.name}")
        out_path = new_path
    if (not is_poly) and ext == ".vtp":
        new_path = out_path.with_suffix(".vtu")
        print(f"warning: output is UnstructuredGrid; writing {new_path.name} instead of {out_path.name}")
        out_path = new_path

    mesh.save(str(out_path))
    return str(out_path)


def write_vtp_from_df(df, path, xyz=("x", "y", "z")):
    """
    Write a VTP point cloud from a dataframe.
    """
    x, y, z = xyz
    pts = np.c_[df[x].to_numpy(float), df[y].to_numpy(float), df[z].to_numpy(float)]
    poly = pv.PolyData(pts)

    for col in df.columns:
        if col in xyz:
            continue
        arr = df[col].to_numpy()
        if np.issubdtype(arr.dtype, np.number):
            poly.point_data[col] = arr

    poly.save(str(path))