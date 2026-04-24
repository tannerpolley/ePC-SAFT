from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any, Iterable


REPO_ROOT = Path(__file__).resolve().parents[1]
PLOTS_ROOT = REPO_ROOT / "docs" / "plots"
PAPER_VALIDATION_SOURCE_ROOT = REPO_ROOT / "scripts" / "paper_validation"
PAPER_VALIDATION_PLOTS_ROOT = PLOTS_ROOT / "paper_validation"
FITS_PLOTS_ROOT = PLOTS_ROOT / "fits"


def _clean_analysis_name(name: str) -> str:
    return name.removesuffix("_analysis")


def paper_validation_path(source_path: str | Path, filename: str | None = None) -> Path:
    source = Path(source_path).resolve()
    source_dir = source if source.is_dir() else source.parent
    rel_dir = source_dir.relative_to(PAPER_VALIDATION_SOURCE_ROOT)
    rel_parts = list(rel_dir.parts)
    if rel_parts:
        rel_parts[0] = _clean_analysis_name(rel_parts[0])
    target_dir = PAPER_VALIDATION_PLOTS_ROOT.joinpath(*rel_parts)
    target = target_dir / (filename if filename is not None else source.name)
    target.parent.mkdir(parents=True, exist_ok=True)
    return target


def paper_validation_dir(source_path: str | Path) -> Path:
    source = Path(source_path).resolve()
    source_dir = source if source.is_dir() else source.parent
    rel_dir = source_dir.relative_to(PAPER_VALIDATION_SOURCE_ROOT)
    rel_parts = list(rel_dir.parts)
    if rel_parts:
        rel_parts[0] = _clean_analysis_name(rel_parts[0])
    target = PAPER_VALIDATION_PLOTS_ROOT.joinpath(*rel_parts)
    target.mkdir(parents=True, exist_ok=True)
    return target


def paper_validation_output_path(path: str | Path) -> Path:
    source = Path(path).resolve()
    if source.is_relative_to(PAPER_VALIDATION_SOURCE_ROOT):
        return paper_validation_path(source.parent, source.name)
    source.parent.mkdir(parents=True, exist_ok=True)
    return source


def fits_plot_path(*parts: str | Path) -> Path:
    target = FITS_PLOTS_ROOT.joinpath(*(str(part) for part in parts))
    target.parent.mkdir(parents=True, exist_ok=True)
    return target


def plot_data_path(image_path: str | Path) -> Path:
    image = Path(image_path)
    return image.parent / "data" / f"{image.stem}_plot_data.csv"


def _format_cell(value: Any) -> str:
    if value is None:
        return ""
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return str(value)
    if math.isnan(numeric):
        return "nan"
    return format(numeric, ".17g")


def _label_for_artist(artist: Any) -> str:
    label = ""
    if hasattr(artist, "get_label"):
        label = str(artist.get_label() or "")
    return "" if label.startswith("_") else label


def _line_rows(fig: Any, image_path: Path) -> Iterable[dict[str, Any]]:
    for axes_index, ax in enumerate(fig.axes):
        axes_title = ax.get_title() if hasattr(ax, "get_title") else ""
        for series_index, line in enumerate(ax.get_lines()):
            x_data = list(line.get_xdata(orig=False))
            y_data = list(line.get_ydata(orig=False))
            for point_index, (x_value, y_value) in enumerate(zip(x_data, y_data, strict=False)):
                yield {
                    "figure_file": image_path.name,
                    "axes_index": axes_index,
                    "axes_title": axes_title,
                    "artist_type": "line",
                    "artist_label": _label_for_artist(line),
                    "series_index": series_index,
                    "point_index": point_index,
                    "x": x_value,
                    "y": y_value,
                    "width": "",
                    "height": "",
                }


def _scatter_rows(fig: Any, image_path: Path) -> Iterable[dict[str, Any]]:
    for axes_index, ax in enumerate(fig.axes):
        axes_title = ax.get_title() if hasattr(ax, "get_title") else ""
        for series_index, collection in enumerate(ax.collections):
            if not hasattr(collection, "get_offsets"):
                continue
            offsets = collection.get_offsets()
            if len(offsets) == 0:
                continue
            for point_index, offset in enumerate(offsets):
                yield {
                    "figure_file": image_path.name,
                    "axes_index": axes_index,
                    "axes_title": axes_title,
                    "artist_type": "scatter",
                    "artist_label": _label_for_artist(collection),
                    "series_index": series_index,
                    "point_index": point_index,
                    "x": offset[0],
                    "y": offset[1],
                    "width": "",
                    "height": "",
                }


def _bar_rows(fig: Any, image_path: Path) -> Iterable[dict[str, Any]]:
    for axes_index, ax in enumerate(fig.axes):
        axes_title = ax.get_title() if hasattr(ax, "get_title") else ""
        for series_index, patch in enumerate(ax.patches):
            required = ("get_x", "get_y", "get_width", "get_height")
            if not all(hasattr(patch, name) for name in required):
                continue
            width = patch.get_width()
            height = patch.get_height()
            if abs(float(width)) <= 1.0e-12 and abs(float(height)) <= 1.0e-12:
                continue
            yield {
                "figure_file": image_path.name,
                "axes_index": axes_index,
                "axes_title": axes_title,
                "artist_type": "bar",
                "artist_label": _label_for_artist(patch),
                "series_index": series_index,
                "point_index": 0,
                "x": patch.get_x() + 0.5 * width,
                "y": patch.get_y(),
                "width": width,
                "height": height,
            }


def export_plot_data(fig: Any, image_path: str | Path) -> Path:
    image = Path(image_path)
    rows = [*_line_rows(fig, image), *_scatter_rows(fig, image), *_bar_rows(fig, image)]
    if not rows:
        rows = [
            {
                "figure_file": image.name,
                "axes_index": "",
                "axes_title": "",
                "artist_type": "no_numeric_artists",
                "artist_label": "",
                "series_index": "",
                "point_index": "",
                "x": "",
                "y": "",
                "width": "",
                "height": "",
            }
        ]
    path = plot_data_path(image)
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "figure_file",
        "axes_index",
        "axes_title",
        "artist_type",
        "artist_label",
        "series_index",
        "point_index",
        "x",
        "y",
        "width",
        "height",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: _format_cell(row.get(field, "")) for field in fieldnames})
    return path


def save_plot_figure(fig: Any, path: str | Path, *, dpi: int = 300, bbox_inches: str | None = "tight", **savefig_kwargs: Any) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    kwargs = dict(savefig_kwargs)
    if bbox_inches is not None:
        kwargs["bbox_inches"] = bbox_inches
    fig.savefig(output_path, dpi=dpi, **kwargs)
    export_plot_data(fig, output_path)
    return output_path
