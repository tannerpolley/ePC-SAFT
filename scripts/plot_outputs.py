from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
from matplotlib import colors as mcolors


REPO_ROOT = Path(__file__).resolve().parents[1]
PLOTS_ROOT = REPO_ROOT / "docs" / "plots"
PAPER_VALIDATION_SOURCE_ROOT = REPO_ROOT / "scripts" / "paper_validation"
PAPER_VALIDATION_PLOTS_ROOT = PAPER_VALIDATION_SOURCE_ROOT
FITS_PLOTS_ROOT = REPO_ROOT / "scripts" / "fits" / "out"
TEST_PLOTS_ROOT = REPO_ROOT / "tests" / "plots" / "out"
OUTPUT_DIR_NAME = "out"


def _clean_analysis_name(name: str) -> str:
    return name.removesuffix("_analysis")


def paper_validation_path(source_path: str | Path, filename: str | None = None) -> Path:
    source = Path(source_path).resolve()
    source_dir = source if source.is_dir() else source.parent
    source_dir.relative_to(PAPER_VALIDATION_SOURCE_ROOT)
    target_dir = source_dir / OUTPUT_DIR_NAME
    target = target_dir / (filename if filename is not None else source.name)
    target.parent.mkdir(parents=True, exist_ok=True)
    return target


def paper_validation_dir(source_path: str | Path) -> Path:
    source = Path(source_path).resolve()
    source_dir = source if source.is_dir() else source.parent
    source_dir.relative_to(PAPER_VALIDATION_SOURCE_ROOT)
    target = source_dir / OUTPUT_DIR_NAME
    target.mkdir(parents=True, exist_ok=True)
    return target


def paper_validation_output_path(path: str | Path) -> Path:
    source = Path(path).resolve()
    if source.is_relative_to(PAPER_VALIDATION_SOURCE_ROOT):
        return paper_validation_path(source.parent, source.name)
    if source.parent.name == OUTPUT_DIR_NAME:
        target = source
    else:
        target = source.parent / OUTPUT_DIR_NAME / source.name
    target.parent.mkdir(parents=True, exist_ok=True)
    return target


def fits_plot_path(*parts: str | Path) -> Path:
    target = FITS_PLOTS_ROOT.joinpath(*(str(part) for part in parts))
    target.parent.mkdir(parents=True, exist_ok=True)
    return target


def _test_plot_category_parts(category: str | Path | Iterable[str | Path] | None) -> list[str] | None:
    if category is None:
        return None
    raw_parts: list[str | Path]
    if isinstance(category, (str, Path)):
        raw_parts = [category]
    else:
        raw_parts = list(category)

    parts: list[str] = []
    for raw_part in raw_parts:
        path_part = Path(raw_part)
        if path_part.is_absolute():
            raise ValueError(f"test plot category must be relative: {raw_part}")
        for part in path_part.parts:
            if part in ("", "."):
                continue
            if part == "..":
                raise ValueError(f"test plot category cannot contain '..': {raw_part}")
            parts.append(part)
    return parts


def test_plot_path(
    source_path: str | Path,
    filename: str | Path,
    *,
    category: str | Path | Iterable[str | Path] | None = None,
) -> Path:
    category_parts = _test_plot_category_parts(category)
    if category_parts is not None:
        target = TEST_PLOTS_ROOT.joinpath(*category_parts) / Path(filename)
        target.parent.mkdir(parents=True, exist_ok=True)
        return target

    source = Path(source_path)
    parts = list(source.parts)
    try:
        tests_index = max(index for index, part in enumerate(parts) if part == "tests")
        rel_parts = parts[tests_index + 1 :]
    except ValueError:
        rel_parts = [source.name]
    if rel_parts and Path(rel_parts[-1]).suffix:
        module_name = Path(rel_parts[-1]).stem
        if module_name.startswith("test_"):
            module_name = module_name.removeprefix("test_")
        rel_parts = [*rel_parts[:-1], module_name]
    target = TEST_PLOTS_ROOT.joinpath(*rel_parts) / Path(filename)
    target.parent.mkdir(parents=True, exist_ok=True)
    return target


def plot_data_path(image_path: str | Path) -> Path:
    image = Path(image_path)
    return image.parent / f"{image.stem}_plot_data.csv"


def plot_svg_path(image_path: str | Path) -> Path:
    image = Path(image_path)
    return image.with_suffix(".svg")


def _strip_trailing_whitespace(path: Path) -> None:
    text = path.read_text(encoding="utf-8")
    normalized = "\n".join(line.rstrip() for line in text.splitlines())
    if text.endswith("\n"):
        normalized += "\n"
    path.write_text(normalized, encoding="utf-8", newline="\n")


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


def _color_cell(value: Any) -> str:
    try:
        return mcolors.to_hex(value, keep_alpha=False)
    except (TypeError, ValueError):
        return ""


def _line_style_cell(line: Any) -> str:
    if not hasattr(line, "get_linestyle"):
        return ""
    style = str(line.get_linestyle() or "")
    return "" if style in {"None", "none"} else style


def _marker_cell(artist: Any) -> str:
    if not hasattr(artist, "get_marker"):
        return ""
    marker = str(artist.get_marker() or "")
    return "" if marker in {"None", "none"} else marker


def _linewidth_cell(artist: Any) -> Any:
    if hasattr(artist, "get_linewidth"):
        return artist.get_linewidth()
    if hasattr(artist, "get_linewidths"):
        widths = artist.get_linewidths()
        if len(widths):
            return widths[0]
    return ""


def _collection_color(collection: Any) -> str:
    if hasattr(collection, "get_facecolors"):
        colors = collection.get_facecolors()
        if len(colors):
            return _color_cell(colors[0])
    if hasattr(collection, "get_edgecolors"):
        colors = collection.get_edgecolors()
        if len(colors):
            return _color_cell(colors[0])
    return ""


def _line_rows(fig: Any, image_path: Path) -> Iterable[dict[str, Any]]:
    for axes_index, ax in enumerate(fig.axes):
        axes_title = ax.get_title() if hasattr(ax, "get_title") else ""
        x_label = ax.get_xlabel() if hasattr(ax, "get_xlabel") else ""
        y_label = ax.get_ylabel() if hasattr(ax, "get_ylabel") else ""
        for series_index, line in enumerate(ax.get_lines()):
            x_data = list(line.get_xdata(orig=False))
            y_data = list(line.get_ydata(orig=False))
            for point_index, (x_value, y_value) in enumerate(zip(x_data, y_data, strict=False)):
                yield {
                    "figure_file": image_path.name,
                    "axes_index": axes_index,
                    "axes_title": axes_title,
                    "x_label": x_label,
                    "y_label": y_label,
                    "artist_type": "line",
                    "artist_label": _label_for_artist(line),
                    "color": _color_cell(line.get_color() if hasattr(line, "get_color") else ""),
                    "linestyle": _line_style_cell(line),
                    "marker": _marker_cell(line),
                    "linewidth": _linewidth_cell(line),
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
        x_label = ax.get_xlabel() if hasattr(ax, "get_xlabel") else ""
        y_label = ax.get_ylabel() if hasattr(ax, "get_ylabel") else ""
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
                    "x_label": x_label,
                    "y_label": y_label,
                    "artist_type": "scatter",
                    "artist_label": _label_for_artist(collection),
                    "color": _collection_color(collection),
                    "linestyle": "",
                    "marker": "",
                    "linewidth": _linewidth_cell(collection),
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
        x_label = ax.get_xlabel() if hasattr(ax, "get_xlabel") else ""
        y_label = ax.get_ylabel() if hasattr(ax, "get_ylabel") else ""
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
                "x_label": x_label,
                "y_label": y_label,
                "artist_type": "bar",
                "artist_label": _label_for_artist(patch),
                "color": _color_cell(patch.get_facecolor() if hasattr(patch, "get_facecolor") else ""),
                "linestyle": "",
                "marker": "",
                "linewidth": _linewidth_cell(patch),
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
                "x_label": "",
                "y_label": "",
                "artist_type": "no_numeric_artists",
                "artist_label": "",
                "color": "",
                "linestyle": "",
                "marker": "",
                "linewidth": "",
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
        "x_label",
        "y_label",
        "artist_type",
        "artist_label",
        "color",
        "linestyle",
        "marker",
        "linewidth",
        "series_index",
        "point_index",
        "x",
        "y",
        "width",
        "height",
    ]
    frame = pd.DataFrame([{field: _format_cell(row.get(field, "")) for field in fieldnames} for row in rows], columns=fieldnames)
    frame.to_csv(path, index=False)
    return path


def save_plot_figure(
    fig: Any,
    path: str | Path,
    *,
    dpi: int = 300,
    bbox_inches: str | None = "tight",
    svg_companion: bool = False,
    **savefig_kwargs: Any,
) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    kwargs = dict(savefig_kwargs)
    if bbox_inches is not None:
        kwargs["bbox_inches"] = bbox_inches
    fig.savefig(output_path, dpi=dpi, **kwargs)
    if svg_companion:
        svg_path = plot_svg_path(output_path)
        svg_kwargs = dict(kwargs)
        fig.savefig(svg_path, format="svg", **svg_kwargs)
        _strip_trailing_whitespace(svg_path)
    export_plot_data(fig, output_path)
    return output_path
