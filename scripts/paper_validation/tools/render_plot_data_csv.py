from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


class PlotDataRenderError(ValueError):
    pass


def _default_output_path(csv_path: Path) -> Path:
    stem = csv_path.stem
    if stem.endswith("_plot_data") and csv_path.parent.name == "data":
        return csv_path.parent.parent / f"{stem.removesuffix('_plot_data')}.png"
    return csv_path.with_suffix(".png")


def _float_cell(row: dict[str, str], field: str, *, required: bool = True) -> float:
    value = row.get(field, "")
    if value == "":
        if required:
            raise PlotDataRenderError(f"Missing numeric field {field!r} in row for {row.get('figure_file', '<unknown>')}.")
        return math.nan
    try:
        numeric = float(value)
    except ValueError as exc:
        raise PlotDataRenderError(f"Non-numeric field {field!r}={value!r} in row for {row.get('figure_file', '<unknown>')}.") from exc
    if not math.isfinite(numeric):
        raise PlotDataRenderError(f"Non-finite field {field!r}={value!r} in row for {row.get('figure_file', '<unknown>')}.")
    return numeric


def _read_rows(csv_path: Path) -> list[dict[str, str]]:
    with csv_path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _axes_key(row: dict[str, str]) -> int:
    raw = row.get("axes_index", "0") or "0"
    try:
        return int(float(raw))
    except ValueError as exc:
        raise PlotDataRenderError(f"Invalid axes_index {raw!r}.") from exc


def _series_key(row: dict[str, str]) -> tuple[int, str, str]:
    return (_axes_key(row), row.get("artist_type", ""), row.get("series_index", "0") or "0")


def _color(row: dict[str, str]) -> str | None:
    color = (row.get("color") or "").strip()
    return color or None


def _label(row: dict[str, str]) -> str | None:
    label = (row.get("artist_label") or "").strip()
    return label or None


def _style(row: dict[str, str]) -> str:
    style = (row.get("linestyle") or "").strip()
    return style or "-"


def _marker(row: dict[str, str]) -> str | None:
    marker = (row.get("marker") or "").strip()
    return marker or None


def render_csv_to_figure(csv_path: Path):
    rows = _read_rows(csv_path)
    if not rows:
        raise PlotDataRenderError(f"{csv_path} has no plot rows.")
    artist_types = {row.get("artist_type", "") for row in rows}
    if artist_types == {"no_numeric_artists"}:
        raise PlotDataRenderError(f"{csv_path} contains no numeric artists to render.")
    unsupported = sorted(artist_types - {"line", "scatter", "bar"})
    if unsupported:
        raise PlotDataRenderError(f"{csv_path} contains unsupported artist_type values: {', '.join(unsupported)}.")

    axes_indexes = sorted({_axes_key(row) for row in rows})
    fig, axes = plt.subplots(len(axes_indexes), 1, figsize=(7.2, max(4.2, 3.2 * len(axes_indexes))), squeeze=False)
    axes_by_index = {axis_index: axes[row_index][0] for row_index, axis_index in enumerate(axes_indexes)}
    grouped: dict[tuple[int, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[_series_key(row)].append(row)

    for (axes_index, artist_type, _series_index), series_rows in sorted(grouped.items()):
        ax = axes_by_index[axes_index]
        series_rows = sorted(series_rows, key=lambda row: float(row.get("point_index", "0") or 0.0))
        first = series_rows[0]
        label = _label(first)
        color = _color(first)
        if artist_type == "line":
            x = [_float_cell(row, "x") for row in series_rows]
            y = [_float_cell(row, "y") for row in series_rows]
            kwargs = {"label": label, "linestyle": _style(first)}
            if color:
                kwargs["color"] = color
            marker = _marker(first)
            if marker:
                kwargs["marker"] = marker
            ax.plot(x, y, **kwargs)
        elif artist_type == "scatter":
            x = [_float_cell(row, "x") for row in series_rows]
            y = [_float_cell(row, "y") for row in series_rows]
            kwargs = {"label": label}
            if color:
                kwargs["color"] = color
            ax.scatter(x, y, **kwargs)
        elif artist_type == "bar":
            x = [_float_cell(row, "x") for row in series_rows]
            y = [_float_cell(row, "y", required=False) for row in series_rows]
            width = [_float_cell(row, "width", required=False) for row in series_rows]
            height = [_float_cell(row, "height") for row in series_rows]
            bottoms = [0.0 if math.isnan(value) else value for value in y]
            widths = [0.8 if math.isnan(value) else value for value in width]
            kwargs = {"label": label}
            if color:
                kwargs["color"] = color
            ax.bar(x, height, width=widths, bottom=bottoms, **kwargs)

    for axis_index, ax in axes_by_index.items():
        axis_rows = [row for row in rows if _axes_key(row) == axis_index]
        first = axis_rows[0]
        if first.get("axes_title"):
            ax.set_title(first["axes_title"])
        if first.get("x_label"):
            ax.set_xlabel(first["x_label"])
        if first.get("y_label"):
            ax.set_ylabel(first["y_label"])
        if any(_label(row) for row in axis_rows):
            ax.legend(loc="best", frameon=False)
        ax.grid(color="0.9", linewidth=0.7)
    fig.tight_layout(pad=1.25)
    return fig


def render_csv_to_static_assets(csv_path: Path, output_path: Path | None = None, *, dpi: int = 140) -> Path:
    csv_path = Path(csv_path)
    output_path = Path(output_path) if output_path is not None else _default_output_path(csv_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig = render_csv_to_figure(csv_path)
    try:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        fig.savefig(output_path.with_suffix(".svg"), format="svg", bbox_inches="tight")
    finally:
        plt.close(fig)
    return output_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Render a canonical plot-data CSV file to static PNG and SVG assets.")
    parser.add_argument("csv_paths", nargs="+", type=Path, help="CSV files using scripts.plot_outputs.export_plot_data schema.")
    parser.add_argument("--output", type=Path, help="Output PNG path. Only valid with one CSV input.")
    parser.add_argument("--dpi", type=int, default=140, help="PNG output DPI.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.output is not None and len(args.csv_paths) != 1:
        raise SystemExit("--output can only be used with one CSV input.")
    for csv_path in args.csv_paths:
        output = render_csv_to_static_assets(csv_path, args.output, dpi=args.dpi)
        print(f"Rendered {csv_path} -> {output} and {output.with_suffix('.svg')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
