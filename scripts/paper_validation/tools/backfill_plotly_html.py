from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import plotly.graph_objects as go
from plotly.subplots import make_subplots


REPO_ROOT = Path(__file__).resolve().parents[3]
PLOTS_ROOT = REPO_ROOT / "docs" / "plots"
DEFAULT_ROOTS = ("fits", "paper_validation")
BACKFILL_MARKER = 'epcsaft-interactive-source="csv_backfill"'
MATPLOTLIB_COLORWAY = (
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
)
LINE_DASH_MAP = {
    "--": "dash",
    "-.": "dashdot",
    ":": "dot",
    "dashed": "dash",
    "dashdot": "dashdot",
    "dotted": "dot",
    "solid": "solid",
    "-": "solid",
}


@dataclass(frozen=True)
class PlotlyBackfillResult:
    candidates: int = 0
    created: int = 0
    skipped: dict[str, int] = field(default_factory=dict)


def _plot_data_path(png_path: Path) -> Path:
    return png_path.parent / "data" / f"{png_path.stem}_plot_data.csv"


def _html_path(png_path: Path) -> Path:
    return png_path.with_suffix(".html")


def _as_float(value: object) -> float | None:
    if value is None:
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    return numeric if math.isfinite(numeric) else None


def _as_int(value: object, default: int = 0) -> int:
    numeric = _as_float(value)
    return default if numeric is None else int(numeric)


def _read_rows(csv_path: Path) -> list[dict[str, str]]:
    with csv_path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def _is_supported_row(row: dict[str, str]) -> bool:
    artist_type = row.get("artist_type", "")
    if artist_type in {"line", "scatter"}:
        return _as_float(row.get("x")) is not None and _as_float(row.get("y")) is not None
    if artist_type == "bar":
        return _as_float(row.get("x")) is not None and _as_float(row.get("height")) is not None
    return False


def _skip_reason(rows: list[dict[str, str]]) -> str | None:
    if not rows:
        return "empty_csv"
    artist_types = {row.get("artist_type", "") for row in rows}
    if artist_types <= {"existing_png_backfill", "no_numeric_artists", ""}:
        return "no_numeric_artists"
    if not any(_is_supported_row(row) for row in rows):
        return "no_supported_rows"
    return None


def _trace_name(row: dict[str, str]) -> str:
    label = (row.get("artist_label") or "").strip()
    if label:
        return label
    artist_type = row.get("artist_type", "trace")
    series_number = _as_int(row.get("series_index")) + 1
    if artist_type == "line":
        return f"Curve {series_number}"
    if artist_type == "scatter":
        return f"Data points {series_number}"
    if artist_type == "bar":
        return f"Bar series {series_number}"
    return f"Trace {series_number}"


def _style_color(row: dict[str, str]) -> str | None:
    color = (row.get("color") or "").strip()
    return color or None


def _line_width(row: dict[str, str]) -> float | None:
    width = _as_float(row.get("linewidth"))
    if width is None or width <= 0.0:
        return None
    return width


def _line_dash(row: dict[str, str]) -> str | None:
    style = (row.get("linestyle") or "").strip()
    return LINE_DASH_MAP.get(style)


def _line_mode(row: dict[str, str], artist_type: str) -> str:
    if artist_type != "line":
        return "markers"
    marker = (row.get("marker") or "").strip()
    return "lines+markers" if marker else "lines"


def _title_from_png(png_path: Path) -> str:
    return png_path.stem.replace("_", " ").replace("-", " ")


def _axis_title(rows: list[dict[str, str]], axes_index: int) -> str:
    for row in rows:
        if _as_int(row.get("axes_index")) == axes_index and (row.get("axes_title") or "").strip():
            return str(row["axes_title"])
    return f"Axis {axes_index + 1}"


def _path_context(png_path: Path, rows: list[dict[str, str]], axes_index: int) -> str:
    title = _axis_title(rows, axes_index)
    return f"{png_path.as_posix()} {png_path.stem} {title}".lower()


def _inferred_axis_label(
    png_path: Path,
    rows: list[dict[str, str]],
    axes_index: int,
    field: str,
    fallback: str,
) -> str:
    context = _path_context(png_path, rows, axes_index)
    if field == "x_label":
        if any(token in context for token in ("miac", "activity coefficient", "osmotic", "molality")):
            return "Molality, m / mol kg^-1"
        if any(token in context for token in ("composition", "mole fraction", "x_")):
            return "Mole fraction, x_i"
        if "density" in context:
            return "Density, rho"
        if any(token in context for token in ("figure_3", "figure 3", "contribution", "bar")):
            return "Reference case index"
        return "Independent variable, x"
    if any(token in context for token in ("miac_m", "maic_m", "mean ionic activity")):
        return "Mean ionic activity coefficient, γ±"
    if any(token in context for token in ("activity coefficient", "miac")):
        return "Activity coefficient, γ"
    if "osmotic" in context:
        return "Osmotic coefficient, φ_osm"
    if "solvation" in context:
        return "Solvation free energy, ΔG^solv"
    if any(token in context for token in ("fugacity", "lnphi", "ln phi")):
        return "Fugacity coefficient, ln φ"
    if any(token in context for token in ("contribution", "born", "debye", "hard-chain", "figure_3", "figure 3")):
        return "Contribution value"
    if "density" in context:
        return "Density or pressure response"
    return fallback


def _axis_label(rows: list[dict[str, str]], axes_index: int, field: str, fallback: str, *, png_path: Path) -> str:
    for row in rows:
        if _as_int(row.get("axes_index")) == axes_index and (row.get(field) or "").strip():
            return str(row[field])
    if fallback in {"x", "value"}:
        return _inferred_axis_label(png_path, rows, axes_index, field, fallback)
    return fallback


def _add_trace(fig: go.Figure, trace: go.BaseTraceType, *, row: int | None, col: int | None) -> None:
    if row is None or col is None:
        fig.add_trace(trace)
    else:
        fig.add_trace(trace, row=row, col=col)


def _update_axes(
    fig: go.Figure,
    *,
    x_title: str,
    y_title: str,
    row: int | None,
    col: int | None,
) -> None:
    if row is None or col is None:
        fig.update_xaxes(title_text=x_title, automargin=True)
        fig.update_yaxes(title_text=y_title, automargin=True)
    else:
        fig.update_xaxes(title_text=x_title, automargin=True, row=row, col=col)
        fig.update_yaxes(title_text=y_title, automargin=True, row=row, col=col)


def _add_group_trace(
    fig: go.Figure,
    rows: list[dict[str, str]],
    *,
    row: int | None,
    col: int | None,
) -> None:
    first = rows[0]
    artist_type = first.get("artist_type", "")
    name = _trace_name(first)
    sorted_rows = sorted(rows, key=lambda item: _as_int(item.get("point_index")))

    if artist_type in {"line", "scatter"}:
        x_values = [_as_float(item.get("x")) for item in sorted_rows]
        y_values = [_as_float(item.get("y")) for item in sorted_rows]
        valid = [(x, y) for x, y in zip(x_values, y_values, strict=False) if x is not None and y is not None]
        if not valid:
            return
        x, y = zip(*valid, strict=True)
        line_style: dict[str, object] = {}
        marker_style: dict[str, object] = {}
        color = _style_color(first)
        if color:
            line_style["color"] = color
            marker_style["color"] = color
        width = _line_width(first)
        if width is not None:
            line_style["width"] = width
        dash = _line_dash(first)
        if dash is not None:
            line_style["dash"] = dash
        trace = go.Scatter(
            x=list(x),
            y=list(y),
            mode=_line_mode(first, artist_type),
            name=name,
            line=line_style or None,
            marker=marker_style or None,
            hovertemplate=f"{name}<br>x=%{{x:.6g}}<br>y=%{{y:.6g}}<extra></extra>",
        )
        _add_trace(fig, trace, row=row, col=col)
        return

    if artist_type == "bar":
        x_values = [_as_float(item.get("x")) for item in sorted_rows]
        heights = [_as_float(item.get("height")) for item in sorted_rows]
        valid_bar = [(x, height) for x, height in zip(x_values, heights, strict=False) if x is not None and height is not None]
        if not valid_bar:
            return
        x, height = zip(*valid_bar, strict=True)
        trace = go.Bar(
            x=list(x),
            y=list(height),
            name=name,
            marker_color=_style_color(first),
            hovertemplate=f"{name}<br>x=%{{x:.6g}}<br>height=%{{y:.6g}}<extra></extra>",
        )
        _add_trace(fig, trace, row=row, col=col)


def figure_from_plot_csv(csv_path: Path, png_path: Path) -> go.Figure | None:
    rows = _read_rows(csv_path)
    if _skip_reason(rows) is not None:
        return None

    axes_indexes = sorted({_as_int(row.get("axes_index")) for row in rows if _is_supported_row(row)})
    if not axes_indexes:
        return None

    if len(axes_indexes) == 1:
        fig = go.Figure()
        axis_positions = {axes_indexes[0]: (None, None)}
    else:
        fig = make_subplots(
            rows=len(axes_indexes),
            cols=1,
            subplot_titles=[_axis_title(rows, axes_index) for axes_index in axes_indexes],
            vertical_spacing=min(0.12, 0.08 + 0.02 * len(axes_indexes)),
        )
        axis_positions = {axes_index: (index + 1, 1) for index, axes_index in enumerate(axes_indexes)}

    grouped: dict[tuple[int, str, int, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if not _is_supported_row(row):
            continue
        axes_index = _as_int(row.get("axes_index"))
        artist_type = row.get("artist_type", "")
        series_index = _as_int(row.get("series_index"))
        label = row.get("artist_label", "")
        grouped[(axes_index, artist_type, series_index, label)].append(row)

    for (axes_index, _artist_type, _series_index, _label), group_rows in sorted(grouped.items()):
        subplot_row, subplot_col = axis_positions[axes_index]
        _add_group_trace(fig, group_rows, row=subplot_row, col=subplot_col)

    if not fig.data:
        return None

    fig.update_layout(
        title=_title_from_png(png_path),
        template="plotly_white",
        colorway=MATPLOTLIB_COLORWAY,
        margin={"l": 76, "r": 34, "t": 88, "b": 104},
        hovermode="closest",
        title_font={"size": 16},
        legend={
            "title": {"text": "Trace"},
            "orientation": "h",
            "yanchor": "top",
            "y": -0.18,
            "xanchor": "center",
            "x": 0.5,
            "font": {"size": 11},
        },
    )
    for axes_index, (subplot_row, subplot_col) in axis_positions.items():
        _update_axes(
            fig,
            x_title=_axis_label(rows, axes_index, "x_label", "x", png_path=png_path),
            y_title=_axis_label(rows, axes_index, "y_label", "value", png_path=png_path),
            row=subplot_row,
            col=subplot_col,
        )
    return fig


def write_backfill_html(fig: go.Figure, html_path: Path) -> Path:
    html_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(
        html_path,
        include_plotlyjs="cdn",
        include_mathjax="cdn",
        full_html=True,
        config={"displaylogo": False, "responsive": True},
    )
    text = html_path.read_text(encoding="utf-8")
    if BACKFILL_MARKER not in text:
        text = text.replace("<head>", f"<head>\n  <meta name=\"epcsaft-interactive-source\" content=\"csv_backfill\">", 1)
        if BACKFILL_MARKER not in text:
            text = f"<!-- {BACKFILL_MARKER} -->\n{text}"
        html_path.write_text(text, encoding="utf-8", newline="\n")
    return html_path


def _iter_candidate_pngs(plots_root: Path, roots: Iterable[str]) -> Iterable[Path]:
    for root_name in roots:
        root = plots_root / root_name
        if not root.exists():
            continue
        yield from sorted(root.rglob("*.png"), key=lambda path: path.relative_to(plots_root).as_posix().lower())


def backfill_plotly_html(
    plots_root: str | Path = PLOTS_ROOT,
    *,
    roots: Iterable[str] = DEFAULT_ROOTS,
    force: bool = False,
    dry_run: bool = False,
) -> PlotlyBackfillResult:
    root = Path(plots_root)
    counters: Counter[str] = Counter()
    created = 0
    candidates = 0

    for png_path in _iter_candidate_pngs(root, roots):
        data_path = _plot_data_path(png_path)
        html_path = _html_path(png_path)
        if not data_path.exists():
            counters["missing_csv"] += 1
            continue
        candidates += 1
        if html_path.exists() and not force:
            counters["existing_html"] += 1
            continue

        rows = _read_rows(data_path)
        reason = _skip_reason(rows)
        if reason is not None:
            counters[reason] += 1
            continue

        fig = figure_from_plot_csv(data_path, png_path)
        if fig is None:
            counters["no_supported_rows"] += 1
            continue
        created += 1
        if not dry_run:
            write_backfill_html(fig, html_path)

    return PlotlyBackfillResult(candidates=candidates, created=created, skipped=dict(sorted(counters.items())))


def _format_result(result: PlotlyBackfillResult, *, dry_run: bool) -> str:
    action = "Would create" if dry_run else "Created"
    skipped = ", ".join(f"{reason}={count}" for reason, count in result.skipped.items()) or "none"
    return f"{action} {result.created} Plotly HTML file(s) from {result.candidates} CSV-backed candidate(s). Skipped: {skipped}."


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Backfill Plotly HTML companions from plot CSV data.")
    parser.add_argument("--root", type=Path, default=PLOTS_ROOT, help="Plot gallery root. Defaults to docs/plots.")
    parser.add_argument(
        "--scan-root",
        action="append",
        dest="scan_roots",
        help="Relative plot subfolder to scan. Can be provided more than once. Defaults to fits and paper_validation.",
    )
    parser.add_argument("--force", action="store_true", help="Overwrite existing HTML companions.")
    parser.add_argument("--dry-run", action="store_true", help="Report work without writing HTML files.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    result = backfill_plotly_html(
        args.root,
        roots=tuple(args.scan_roots or DEFAULT_ROOTS),
        force=args.force,
        dry_run=args.dry_run,
    )
    print(_format_result(result, dry_run=args.dry_run))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
