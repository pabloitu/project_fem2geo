
import argparse
from pathlib import Path

from .core import read_config, project_mesh, project_points, write_vtp_from_df


def main(argv=None):
    """
    Command line interface for project_fem2geo.
    """
    p = argparse.ArgumentParser(prog="project_fem2geo")
    sub = p.add_subparsers(dest="cmd", required=True)

    p_points = sub.add_parser("points", help="Transform a CSV point table")
    p_points.add_argument("config", help="Path to config.yml")
    p_points.add_argument("--out", default=None, help="Override output path")

    p_mesh = sub.add_parser("mesh", help="Transform a mesh (vtp/vtu/vtk)")
    p_mesh.add_argument("config", help="Path to config.yml")
    p_mesh.add_argument("--out", default=None, help="Override output path")

    args = p.parse_args(argv)
    cfg = read_config(args.config)

    base = Path(Path(args.config).resolve().parent)

    if args.cmd == "points":
        df = project_points(cfg)

        out_cfg = cfg.get("output", {})
        fmt = str(out_cfg.get("format", "vtp")).strip().lower()

        out_file = args.out if args.out else out_cfg.get("file", "projected_points.vtp")
        out_path = Path(out_file)
        if not out_path.is_absolute():
            out_path = base / out_path
        out_path.parent.mkdir(parents=True, exist_ok=True)

        if fmt == "csv":
            df.to_csv(out_path, index=False)
        elif fmt == "vtp":
            write_vtp_from_df(df, out_path, xyz=("x", "y", "z"))
        else:
            raise ValueError("output.format must be 'csv' or 'vtp'.")

        return 0

    if args.cmd == "mesh":
        out_cfg = cfg.get("output", {})
        out_file = args.out if args.out else out_cfg.get("file", "out.vtp")
        out_path = Path(out_file)
        if not out_path.is_absolute():
            out_path = base / out_path
        out_path.parent.mkdir(parents=True, exist_ok=True)

        cfg["output"]["file"] = str(out_path)
        project_mesh(cfg)
        return 0

    return 1