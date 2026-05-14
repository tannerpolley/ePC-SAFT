from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
from matplotlib import colors as mcolors

REPO_ROOT = Path(__file__).resolve().parents[1]
ANALYSES_ROOT = REPO_ROOT / "analyses"
PAPER_VALIDATION_SOURCE_ROOT = ANALYSES_ROOT / "paper_validation"
PAPER_VALIDATION_NATIVE_ROOT = PAPER_VALIDATION_SOURCE_ROOT / "native"
PAPER_VALIDATION_APPLICATION_ROOT = PAPER_VALIDATION_SOURCE_ROOT / "application"
DATA_VALIDATION_ROOT = ANALYSES_ROOT / "data_validation"
PACKAGE_VALIDATION_ROOT = ANALYSES_ROOT / "package_validation"
FITS_ANALYSIS_ROOT = DATA_VALIDATION_ROOT / "miac_fits"
FITS_CATEGORY_ROOTS = {
    "dielectric": DATA_VALIDATION_ROOT / "dielectric_fits",
    "osmotic": DATA_VALIDATION_ROOT / "osmotic_validation",
}
TEST_PLOTS_ANALYSIS_ROOT = PACKAGE_VALIDATION_ROOT / "package_plot_smokes"
TEST_PLOT_CATEGORY_ROOTS = {
    ("api", "parity"): "api_parity",
    ("contributions", "neutral"): "contribution_outputs",
    ("contributions", "ionic"): "contribution_outputs",
    ("equilibrium", "vle"): "equilibrium_outputs",
    ("equilibrium", "lle"): "equilibrium_outputs",
    ("equilibrium", "electrolyte_lle"): "equilibrium_outputs",
    ("equilibrium", "stability"): "equilibrium_outputs",
    ("native", "branches"): "native_outputs",
    ("native", "derivatives"): "native_outputs",
    ("properties", "residual_energy"): "property_outputs",
    ("properties", "activity_fugacity"): "property_outputs",
    ("regression", "hydrocarbon"): "regression_outputs",
    ("regression", "gradients"): "regression_outputs",
    ("regression", "binary_ethanol_water"): "regression_outputs",
}
RESULTS_DIR_NAME = "results"
RUNS_DIR_NAME = "runs"
FIGURES_DIR_NAME = "figures"
FIGURE_INPUT_DIR_NAME = "input"
FIGURE_OUTPUT_DIR_NAME = "output"
FIGURE_SCRIPTS_DIR_NAME = "scripts"
ANALYSIS_METADATA_NAME = "analysis.yaml"


def _clean_analysis_name(name: str) -> str:
    return name.removesuffix("_analysis")


def _analysis_root_for(source_path: str | Path) -> Path:
    source = Path(source_path).resolve()
    source_dir = source if source.is_dir() else source.parent
    try:
        source_dir.relative_to(ANALYSES_ROOT)
    except ValueError as exc:
        raise ValueError(f"path is not inside analyses/: {source}") from exc

    for candidate in (source_dir, *source_dir.parents):
        if candidate == ANALYSES_ROOT.parent:
            break
        if (candidate / ANALYSIS_METADATA_NAME).is_file():
            return candidate

    raise ValueError(f"path does not include an analysis metadata root: {source}")


def analysis_root(source_path: str | Path) -> Path:
    return _analysis_root_for(source_path)


def repo_root_for(source_path: str | Path) -> Path:
    source = Path(source_path).resolve()
    try:
        source.relative_to(REPO_ROOT)
    except ValueError as exc:
        raise ValueError(f"path is not inside repo root {REPO_ROOT}: {source}") from exc
    return REPO_ROOT


def _figure_root_for(source_path: str | Path) -> Path | None:
    source = Path(source_path).resolve()
    source_dir = source if source.is_dir() else source.parent
    analysis_root = _analysis_root_for(source_dir)
    try:
        relative = source_dir.relative_to(analysis_root)
    except ValueError:
        return None
    parts = relative.parts
    if len(parts) >= 2 and parts[0] == FIGURES_DIR_NAME:
        return analysis_root / FIGURES_DIR_NAME / parts[1]
    return None


def figure_root_dir(source_path: str | Path) -> Path:
    figure_root = _figure_root_for_any(source_path)
    if figure_root is None:
        raise ValueError(f"path is not inside analyses/*/figures/*: {Path(source_path).resolve()}")
    return figure_root


def analysis_scripts_dir(source_path: str | Path) -> Path:
    return analysis_root(source_path) / FIGURE_SCRIPTS_DIR_NAME


def _figure_root_for_script(source_path: str | Path) -> Path | None:
    source = Path(source_path).resolve()
    source_dir = source if source.is_dir() else source.parent
    analysis_root = _analysis_root_for(source_dir)
    scripts_root = analysis_root / "scripts"
    try:
        relative = source_dir.relative_to(scripts_root)
    except ValueError:
        return None
    parts = [part for part in relative.parts if part not in ("", ".")]
    if not parts:
        return None
    return analysis_root / FIGURES_DIR_NAME / Path(*parts)


def _figure_root_for_any(source_path: str | Path) -> Path | None:
    return _figure_root_for(source_path) or _figure_root_for_script(source_path)


def _relative_script_parts(source_path: str | Path) -> list[str]:
    source = Path(source_path).resolve()
    source_dir = source if source.is_dir() else source.parent
    analysis_root = _analysis_root_for(source_dir)
    try:
        relative = source_dir.relative_to(analysis_root / "scripts")
    except ValueError:
        relative = source_dir.relative_to(analysis_root)
    return [part for part in relative.parts if part not in ("", ".")]


def _is_placeholder_filename(filename: str | Path) -> bool:
    return Path(filename).stem.startswith("_placeholder")


def _plot_set_dir(root: Path, parts: Iterable[str | Path], filename: str | Path | None = None) -> Path:
    clean_parts: list[str] = []
    for raw_part in parts:
        part = Path(raw_part)
        if part.is_absolute():
            raise ValueError(f"plot-set path part must be relative: {raw_part}")
        for path_part in part.parts:
            if path_part in ("", "."):
                continue
            if path_part == "..":
                raise ValueError(f"plot-set path part cannot contain '..': {raw_part}")
            clean_parts.append(path_part)

    if filename is not None and not _is_placeholder_filename(filename):
        stem = Path(filename).stem
        if stem and (not clean_parts or clean_parts[-1] != stem):
            clean_parts.append(stem)

    target = root / RESULTS_DIR_NAME
    if clean_parts:
        target = target.joinpath(*clean_parts)
    target.mkdir(parents=True, exist_ok=True)
    return target


def analysis_plot_set_dir(
    source_path: str | Path,
    filename: str | Path | None = None,
    *,
    category: str | Path | Iterable[str | Path] | None = None,
) -> Path:
    figure_root = _figure_root_for_any(source_path)
    if figure_root is not None:
        target = figure_root / FIGURE_OUTPUT_DIR_NAME
        if category is not None:
            if isinstance(category, (str, Path)):
                category_parts = [category]
            else:
                category_parts = list(category)
            if category_parts:
                target = target.joinpath(*[part for part in category_parts if str(part) not in ("", ".")])
        target.mkdir(parents=True, exist_ok=True)
        return target

    analysis_root = _analysis_root_for(source_path)
    if category is None:
        plot_set_parts: list[str | Path] = _relative_script_parts(source_path)
    elif isinstance(category, (str, Path)):
        plot_set_parts = [category]
    else:
        plot_set_parts = list(category)
    return _plot_set_dir(analysis_root, plot_set_parts, filename)


def _data_kind_parts(kind: str | Path | Iterable[str | Path]) -> list[str]:
    if isinstance(kind, (str, Path)):
        raw_parts = [kind]
    else:
        raw_parts = list(kind)

    parts: list[str] = []
    for raw_part in raw_parts:
        path_part = Path(raw_part)
        if path_part.is_absolute():
            raise ValueError(f"data-path kind must be relative: {raw_part}")
        for part in path_part.parts:
            if part in ("", "."):
                continue
            if part == "..":
                raise ValueError(f"data-path kind cannot contain '..': {raw_part}")
            parts.append(part)
    return parts


def analysis_data_dir(
    source_path: str | Path,
    *,
    kind: str | Path | Iterable[str | Path] = "input",
    category: str | Path | Iterable[str | Path] | None = None,
) -> Path:
    kind_parts = _data_kind_parts(kind)
    figure_root = _figure_root_for_any(source_path)
    explicit_category = category is not None
    if category is None:
        category_parts: list[str | Path] = []
    elif isinstance(category, (str, Path)):
        category_parts = [category]
    else:
        category_parts = list(category)

    if figure_root is not None:
        if kind_parts and kind_parts[0] == "input":
            target = figure_root / FIGURE_INPUT_DIR_NAME
            kind_tail = kind_parts[1:]
        else:
            target = figure_root / FIGURE_OUTPUT_DIR_NAME
            kind_tail = kind_parts
        if kind_tail:
            target = target.joinpath(*kind_tail)
        if explicit_category and category_parts:
            target = target.joinpath(*[part for part in category_parts if str(part) not in ("", ".")])
        target.mkdir(parents=True, exist_ok=True)
        return target

    analysis_root = _analysis_root_for(source_path)
    target = analysis_root / "data"
    if kind_parts:
        target = target.joinpath(*kind_parts)
    if category_parts or not explicit_category:
        if not explicit_category:
            category_parts = _relative_script_parts(source_path)
    if category_parts:
        target = target.joinpath(*[part for part in category_parts if str(part) not in ("", ".")])
    target.mkdir(parents=True, exist_ok=True)
    return target


def analysis_data_path(
    source_path: str | Path,
    filename: str | Path,
    *,
    kind: str | Path | Iterable[str | Path] = "input",
    category: str | Path | Iterable[str | Path] | None = None,
) -> Path:
    target = analysis_data_dir(source_path, kind=kind, category=category)
    path = target / Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def analysis_runs_dir(
    source_path: str | Path,
    *,
    category: str | Path | Iterable[str | Path] | None = None,
) -> Path:
    explicit_category = category is not None
    if category is None:
        category_parts: list[str | Path] = []
    elif isinstance(category, (str, Path)):
        category_parts = [category]
    else:
        category_parts = list(category)

    figure_root = _figure_root_for_any(source_path)
    if figure_root is not None:
        target = figure_root / FIGURE_OUTPUT_DIR_NAME / RUNS_DIR_NAME
        if explicit_category and category_parts:
            target = target.joinpath(*[part for part in category_parts if str(part) not in ("", ".")])
        target.mkdir(parents=True, exist_ok=True)
        return target

    analysis_root = _analysis_root_for(source_path)
    target = analysis_root / RESULTS_DIR_NAME / RUNS_DIR_NAME
    if category_parts or not explicit_category:
        if not explicit_category:
            category_parts = _relative_script_parts(source_path)
    if category_parts:
        target = target.joinpath(*[part for part in category_parts if str(part) not in ("", ".")])
    target.mkdir(parents=True, exist_ok=True)
    return target


def analysis_runs_path(
    source_path: str | Path,
    filename: str | Path,
    *,
    category: str | Path | Iterable[str | Path] | None = None,
) -> Path:
    target = analysis_runs_dir(source_path, category=category)
    path = target / Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def figure_input_dir(source_path: str | Path) -> Path:
    figure_root = _figure_root_for_any(source_path)
    if figure_root is None:
        raise ValueError(f"path is not inside analyses/*/figures/*: {Path(source_path).resolve()}")
    target = figure_root / FIGURE_INPUT_DIR_NAME
    target.mkdir(parents=True, exist_ok=True)
    return target


def figure_input_path(source_path: str | Path, filename: str | Path) -> Path:
    target = figure_input_dir(source_path) / Path(filename)
    target.parent.mkdir(parents=True, exist_ok=True)
    return target


def figure_output_dir(source_path: str | Path) -> Path:
    figure_root = _figure_root_for_any(source_path)
    if figure_root is None:
        raise ValueError(f"path is not inside analyses/*/figures/*: {Path(source_path).resolve()}")
    target = figure_root / FIGURE_OUTPUT_DIR_NAME
    target.mkdir(parents=True, exist_ok=True)
    return target


def figure_output_path(source_path: str | Path, filename: str | Path | None = None) -> Path:
    source = Path(source_path).resolve()
    target = figure_output_dir(source)
    name = Path(filename) if filename is not None else source.name
    path = target / name
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def analysis_final_dir(source_path: str | Path, category: str = "figures") -> Path:
    """Compatibility alias for analysis-owned curated result directories.

    New analysis outputs should use figure-owned ``figures/<figure_id>/output``
    folders. Historical analyses that still use analysis-level curated result
    directories map this helper to ``results/<category>/`` instead of the
    older ``results/final/<category>/`` layout.
    """

    return analysis_plot_set_dir(source_path, category=category)


def analysis_final_path(source_path: str | Path, filename: str | Path, category: str = "figures") -> Path:
    return analysis_plot_set_dir(source_path, filename=filename, category=category) / Path(filename)


def paper_validation_path(source_path: str | Path, filename: str | None = None) -> Path:
    source = Path(source_path).resolve()
    figure_root = _figure_root_for_any(source)
    if figure_root is not None:
        return figure_output_path(source, filename if filename is not None else source.name)
    target = analysis_plot_set_dir(source, filename if filename is not None else source.name) / Path(
        filename if filename is not None else source.name
    )
    target.parent.mkdir(parents=True, exist_ok=True)
    return target


def paper_validation_dir(source_path: str | Path) -> Path:
    figure_root = _figure_root_for_any(source_path)
    if figure_root is not None:
        return figure_output_dir(source_path)
    return analysis_plot_set_dir(source_path)


def paper_validation_output_path(path: str | Path) -> Path:
    source = Path(path).resolve()
    figure_root = _figure_root_for_any(source)
    if figure_root is not None:
        if FIGURE_OUTPUT_DIR_NAME in source.parts:
            target = source
        else:
            target = figure_output_path(source)
        target.parent.mkdir(parents=True, exist_ok=True)
        return target
    if source.is_relative_to(ANALYSES_ROOT):
        if RESULTS_DIR_NAME in source.parts:
            target = source
        else:
            target = analysis_plot_set_dir(source.parent, source.name) / source.name
        target.parent.mkdir(parents=True, exist_ok=True)
        return target
    if RESULTS_DIR_NAME in source.parts:
        target = source
    else:
        target = _plot_set_dir(source.parent, [source.stem], source.name) / source.name
    target.parent.mkdir(parents=True, exist_ok=True)
    return target


def fits_plot_path(*parts: str | Path) -> Path:
    raw_parts = [str(part) for part in parts]
    analysis_root = FITS_CATEGORY_ROOTS.get(raw_parts[0], FITS_ANALYSIS_ROOT) if raw_parts else FITS_ANALYSIS_ROOT
    figure_root = analysis_root / FIGURES_DIR_NAME
    if raw_parts and Path(raw_parts[-1]).suffix:
        filename = Path(raw_parts[-1])
        figure_parts = raw_parts[:-1]
        target = figure_root
        if figure_parts:
            target = target.joinpath(*figure_parts)
        target = target / FIGURE_OUTPUT_DIR_NAME / filename
    else:
        target = figure_root
        if raw_parts:
            target = target.joinpath(*raw_parts)
        target = target / FIGURE_OUTPUT_DIR_NAME
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
        target = TEST_PLOTS_ANALYSIS_ROOT / FIGURES_DIR_NAME
        mapped_figure_id = TEST_PLOT_CATEGORY_ROOTS.get(tuple(category_parts))
        if mapped_figure_id is not None:
            target = target / mapped_figure_id / FIGURE_OUTPUT_DIR_NAME / Path(filename)
            target.parent.mkdir(parents=True, exist_ok=True)
            return target
        if category_parts:
            target = target.joinpath(*category_parts)
        target = target / FIGURE_OUTPUT_DIR_NAME / Path(filename)
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
    target = TEST_PLOTS_ANALYSIS_ROOT / FIGURES_DIR_NAME
    if rel_parts:
        target = target.joinpath(*rel_parts)
    target = target / FIGURE_OUTPUT_DIR_NAME / Path(filename)
    target.parent.mkdir(parents=True, exist_ok=True)
    return target


def plot_data_path(image_path: str | Path) -> Path:
    image = Path(image_path)
    return image.parent / f"{image.stem}.csv"


def plot_svg_path(image_path: str | Path) -> Path:
    image = Path(image_path)
    return image.with_suffix(".svg")


def plot_style_path(image_path: str | Path) -> Path:
    image = Path(image_path)
    return image.parent / f"{image.stem}.mpl.yaml"


def _strip_trailing_whitespace(path: Path) -> None:
    text = path.read_text(encoding="utf-8")
    normalized = "\n".join(line.rstrip() for line in text.splitlines())
    if text.endswith("\n"):
        normalized += "\n"
    path.write_text(normalized, encoding="utf-8", newline="\n")


def _yaml_scalar(value: Any) -> str:
    text = "" if value is None else str(value)
    return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'


def write_mpl_style_contract(fig: Any, image_path: str | Path) -> Path:
    image = Path(image_path)
    path = plot_style_path(image)
    axes_blocks: list[str] = []
    for index, ax in enumerate(getattr(fig, "axes", [])):
        axes_blocks.extend(
            [
                f"  - index: {index}",
                f"    title: {_yaml_scalar(ax.get_title() if hasattr(ax, 'get_title') else '')}",
                f"    xlabel: {_yaml_scalar(ax.get_xlabel() if hasattr(ax, 'get_xlabel') else '')}",
                f"    ylabel: {_yaml_scalar(ax.get_ylabel() if hasattr(ax, 'get_ylabel') else '')}",
            ]
        )
    if not axes_blocks:
        axes_blocks.append("  []")

    content = [
        "# Matplotlib plot-set style contract.",
        "# Edit this sidecar, then rerun the owning render script.",
        "figure:",
        f"  file: {_yaml_scalar(image.name)}",
        "  format: " + _yaml_scalar(image.suffix.lstrip(".")),
        "axes:",
        *axes_blocks,
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        path.write_text("\n".join(content), encoding="utf-8", newline="\n")
    return path


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
    frame = pd.DataFrame(
        [{field: _format_cell(row.get(field, "")) for field in fieldnames} for row in rows], columns=fieldnames
    )
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
    write_mpl_style_contract(fig, output_path)
    return output_path
