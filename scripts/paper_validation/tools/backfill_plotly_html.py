from __future__ import annotations

import argparse
import csv
import math
import re
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


def _series_key(row: dict[str, str]) -> int:
    return _as_int(row.get("series_index"))


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


def _latex_to_plotly_text(label: str) -> str:
    replacements = {
        r"$\gamma_{\pm}^{m,*}$": "γ±",
        r"$\gamma_{\pm}$": "γ±",
        r"$\phi_m$": "φ_m",
        r"$\Delta G_i^{solv,\infty,x}$": "ΔG_i^solv,∞,x",
        r"$\Delta G_{\mathrm{hyd},i}^{\infty}$": "ΔG_hyd,i^∞",
        r"$m$": "m",
        r"mol kg$^{-1}$": "mol kg⁻¹",
        r"mol$^{-1}$": "mol⁻¹",
    }
    out = label
    for source, replacement in replacements.items():
        out = out.replace(source, replacement)
    out = re.sub(r"\$([A-Za-z0-9_{}\\,^+\-\s]+)\$", lambda match: match.group(1), out)
    out = out.replace(r"\mathrm{", "").replace("}", "").replace("{", "").replace("\\", "")
    return out


def _context_text(png_path: Path, rows: list[dict[str, str]], axes_index: int = 0) -> str:
    return f"{png_path.as_posix()} {png_path.stem} {_axis_title(rows, axes_index)}".lower()


def _presentation_series_labels(context: str) -> list[str] | None:
    rules: tuple[tuple[tuple[str, ...], list[str]], ...] = (
        (("sodium salts in ethanol", "ethanol_sodium_salts"), ["NaCl", "NaBr", "NaI"]),
        (("chlorides in methanol", "methanol_chlorides"), ["LiCl", "NaCl", "KCl"]),
        (("iodides in methanol", "methanol_iodides"), ["LiI", "NaI", "KI"]),
        (("potassium halides in water", "water_potassium_halides"), ["KCl", "KBr", "KI"]),
        (("nacl in water, methanol, and ethanol", "nacl_solvents"), ["Water", "Methanol", "Ethanol"]),
        (("nabr in water, methanol, and ethanol", "nabr_solvents"), ["Water", "Methanol", "Ethanol"]),
        (("lii in water, methanol, and ethanol", "lii_solvents"), ["Water", "Methanol", "Ethanol"]),
        (("libr in methanol and ethanol", "libr_nonaqueous"), ["Methanol", "Ethanol"]),
    )
    for tokens, labels in rules:
        if any(token in context for token in tokens):
            return labels
    return None


def _formula_from_context(context: str) -> str | None:
    for formula in ("LiCl", "LiBr", "LiI", "NaCl", "NaBr", "NaI", "KCl", "KBr", "KI"):
        if formula.lower() in context:
            return formula
    words = {
        "sodium chloride": "NaCl",
        "sodium bromide": "NaBr",
        "sodium iodide": "NaI",
        "lithium chloride": "LiCl",
        "lithium bromide": "LiBr",
        "lithium iodide": "LiI",
        "potassium chloride": "KCl",
        "potassium bromide": "KBr",
        "potassium iodide": "KI",
    }
    for phrase, formula in words.items():
        if phrase in context:
            return formula
    return None


def _inferred_series_base(row: dict[str, str], rows: list[dict[str, str]], png_path: Path) -> str | None:
    label = (row.get("artist_label") or "").strip()
    if label:
        return _latex_to_plotly_text(label)
    axes_index = _as_int(row.get("axes_index"))
    context = _context_text(png_path, rows, axes_index)
    labels = _presentation_series_labels(context)
    series_index = _series_key(row)
    if labels is not None and 0 <= series_index < len(labels):
        return labels[series_index]
    formula = _formula_from_context(context)
    if formula is not None:
        return formula
    return None


def _trace_name(row: dict[str, str], rows: list[dict[str, str]], png_path: Path, *, bar_group_index: int | None = None) -> str:
    artist_type = row.get("artist_type", "trace")
    if artist_type == "bar" and bar_group_index is not None:
        label = _bar_group_label(rows, png_path, _as_int(row.get("axes_index")), bar_group_index)
        if label is not None:
            return label
    explicit_label = (row.get("artist_label") or "").strip()
    base = _inferred_series_base(row, rows, png_path)
    if explicit_label:
        return base or _latex_to_plotly_text(explicit_label)
    series_number = _series_key(row) + 1
    if base:
        if artist_type == "line":
            return f"{base} fit"
        if artist_type == "scatter":
            return f"{base} data"
        return base
    if artist_type == "line":
        return f"Fit {series_number}"
    if artist_type == "scatter":
        return f"Data {series_number}"
    if artist_type == "bar":
        return f"Set {series_number}"
    return f"Series {series_number}"


def _style_color(row: dict[str, str]) -> str | None:
    color = (row.get("color") or "").strip()
    return color or None


def _fallback_series_color(row: dict[str, str]) -> str:
    return MATPLOTLIB_COLORWAY[_series_key(row) % len(MATPLOTLIB_COLORWAY)]


def _trace_color(row: dict[str, str]) -> str:
    return _style_color(row) or _fallback_series_color(row)


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


def _inferred_axis_label(
    png_path: Path,
    rows: list[dict[str, str]],
    axes_index: int,
    field: str,
    fallback: str,
) -> str:
    context = _context_text(png_path, rows, axes_index)
    if field == "x_label":
        if any(token in context for token in ("miac", "activity coefficient", "osmotic", "molality")):
            return "m / mol kg⁻¹"
        if any(token in context for token in ("composition", "mole fraction", "x_")):
            return "x_i"
        if "density" in context:
            return "ρ"
        if any(token in context for token in ("figure_3", "figure 3", "contribution", "bar")):
            return "Reference case index"
        return "Independent variable, x"
    if any(token in context for token in ("miac_m", "maic_m", "mean ionic activity")):
        return "γ±"
    if any(token in context for token in ("activity coefficient", "miac")):
        return "γ"
    if "osmotic" in context:
        return "φ_m"
    if "solvation" in context:
        return "ΔG_i^solv,∞"
    if any(token in context for token in ("fugacity", "lnphi", "ln phi")):
        return "ln φ_i"
    if any(token in context for token in ("contribution", "born", "debye", "hard-chain", "figure_3", "figure 3")):
        return "Contribution value"
    if "density" in context:
        return "ρ or P"
    return fallback


def _axis_label(rows: list[dict[str, str]], axes_index: int, field: str, fallback: str, *, png_path: Path) -> str:
    for row in rows:
        if _as_int(row.get("axes_index")) == axes_index and (row.get(field) or "").strip():
            return _latex_to_plotly_text(str(row[field]))
    if fallback in {"x", "value"}:
        return _inferred_axis_label(png_path, rows, axes_index, field, fallback)
    return fallback


def _supported_rows_for_axis(rows: list[dict[str, str]], axes_index: int, artist_type: str | None = None) -> list[dict[str, str]]:
    return [
        row
        for row in rows
        if _is_supported_row(row)
        and _as_int(row.get("axes_index")) == axes_index
        and (artist_type is None or row.get("artist_type") == artist_type)
    ]


def _bar_offset(row: dict[str, str]) -> float:
    x = _as_float(row.get("x")) or 0.0
    return round(x - round(x), 6)


def _bar_offsets(rows: list[dict[str, str]], axes_index: int) -> list[float]:
    return sorted({_bar_offset(row) for row in _supported_rows_for_axis(rows, axes_index, "bar")})


def _bar_group_index(row: dict[str, str], rows: list[dict[str, str]]) -> int | None:
    if row.get("artist_type") != "bar" or (row.get("artist_label") or "").strip():
        return None
    offsets = _bar_offsets(rows, _as_int(row.get("axes_index")))
    if len(offsets) <= 1:
        return None
    try:
        return offsets.index(_bar_offset(row))
    except ValueError:
        return None


def _bar_group_label(rows: list[dict[str, str]], png_path: Path, axes_index: int, group_index: int) -> str | None:
    context = _context_text(png_path, rows, axes_index)
    if "figure_4" in context and "solvation" in context:
        labels = ["Literature", "ePC-SAFT 2025", "ePC-SAFT 2020"]
        return labels[group_index] if group_index < len(labels) else None
    if "figure_3_total" in context or "total check" in context:
        labels = [
            "Figure 2 total (paper)",
            "Figure 2 total (ePC-SAFT)",
            "Figure 3 sum (paper)",
            "Figure 3 sum (ePC-SAFT)",
        ]
        return labels[group_index] if group_index < len(labels) else None
    return f"Set {group_index + 1}"


def _bar_category_labels(rows: list[dict[str, str]], png_path: Path, axes_index: int) -> list[str] | None:
    bar_rows = _supported_rows_for_axis(rows, axes_index, "bar")
    if not bar_rows:
        return None
    category_count = len({_bar_category_index(row) for row in bar_rows})
    context = _context_text(png_path, rows, axes_index)
    if "solvation" in context:
        labels = ["Li+", "Na+", "K+", "Cl-", "Br-", "I-"]
        return labels[:category_count] if category_count <= len(labels) else None
    return None


def _bar_category_index(row: dict[str, str]) -> int:
    return int(round((_as_float(row.get("x")) or 0.0) - _bar_offset(row)))


def _bar_x_value(row: dict[str, str], rows: list[dict[str, str]], png_path: Path) -> float | str:
    axes_index = _as_int(row.get("axes_index"))
    category_index = _bar_category_index(row)
    labels = _bar_category_labels(rows, png_path, axes_index)
    if labels is not None and 0 <= category_index < len(labels):
        return labels[category_index]
    return category_index


def _hide_unlabeled_geometry(row: dict[str, str], rows: list[dict[str, str]], png_path: Path) -> bool:
    if row.get("artist_type") != "line" or (row.get("artist_label") or "").strip():
        return False
    context = _context_text(png_path, rows, _as_int(row.get("axes_index")))
    if "2026_khudaida" not in context:
        return False
    axis_rows = _supported_rows_for_axis(rows, _as_int(row.get("axes_index")), "line")
    return len(axis_rows) > 20


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
    all_rows: list[dict[str, str]],
    png_path: Path,
    row: int | None,
    col: int | None,
) -> None:
    first = rows[0]
    artist_type = first.get("artist_type", "")
    bar_group = _bar_group_index(first, all_rows)
    name = _trace_name(first, all_rows, png_path, bar_group_index=bar_group)
    sorted_rows = sorted(rows, key=lambda item: _as_int(item.get("point_index")))
    axes_index = _as_int(first.get("axes_index"))
    x_label = _axis_label(all_rows, axes_index, "x_label", "x", png_path=png_path)
    y_label = _axis_label(all_rows, axes_index, "y_label", "value", png_path=png_path)
    showlegend = not _hide_unlabeled_geometry(first, all_rows, png_path)

    if artist_type in {"line", "scatter"}:
        x_values = [_as_float(item.get("x")) for item in sorted_rows]
        y_values = [_as_float(item.get("y")) for item in sorted_rows]
        valid = [(x, y) for x, y in zip(x_values, y_values, strict=False) if x is not None and y is not None]
        if not valid:
            return
        x, y = zip(*valid, strict=True)
        line_style: dict[str, object] = {}
        marker_style: dict[str, object] = {}
        color = _trace_color(first)
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
            showlegend=showlegend,
            hovertemplate=f"{name}<br>{x_label}=%{{x:.6g}}<br>{y_label}=%{{y:.6g}}<extra></extra>",
        )
        _add_trace(fig, trace, row=row, col=col)
        return

    if artist_type == "bar":
        if bar_group is not None:
            sorted_rows = sorted(rows, key=_bar_category_index)
            x_values = [_bar_x_value(item, all_rows, png_path) for item in sorted_rows]
        else:
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
            marker_color=_trace_color(first),
            hovertemplate=f"{name}<br>{x_label}=%{{x}}<br>{y_label}=%{{y:.6g}}<extra></extra>",
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

    grouped: dict[tuple[int, str, int, str | int | None], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if not _is_supported_row(row):
            continue
        axes_index = _as_int(row.get("axes_index"))
        artist_type = row.get("artist_type", "")
        bar_group = _bar_group_index(row, rows)
        series_index = bar_group if bar_group is not None else _as_int(row.get("series_index"))
        label = row.get("artist_label", "")
        group_label = label if bar_group is None else f"bar:{bar_group}"
        grouped[(axes_index, artist_type, series_index, group_label)].append(row)

    for (axes_index, _artist_type, _series_index, _label), group_rows in sorted(grouped.items()):
        subplot_row, subplot_col = axis_positions[axes_index]
        _add_group_trace(fig, group_rows, all_rows=rows, png_path=png_path, row=subplot_row, col=subplot_col)

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
