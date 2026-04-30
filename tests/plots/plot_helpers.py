from __future__ import annotations

from pathlib import Path
import re
import textwrap
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["svg.hashsalt"] = "epcsaft-test-plots"

from epcsaft import ePCSAFTMixture
from epcsaft.regression import _debug_native_pure_neutral_objective
from scripts import plot_outputs
from tests.helpers.regression_cases import _methane_like_records
from tests.helpers.regression_cases import _minimal_neutral_metadata

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


def hydrocarbon_basis_mixture() -> ePCSAFTMixture:
    params = {
        "m": np.asarray([1.0, 1.6069, 2.0020]),
        "s": np.asarray([3.7039, 3.5206, 3.6184]),
        "e": np.asarray([150.03, 191.42, 208.11]),
        "k_ij": np.asarray(
            [
                [0.0, 3.0e-4, 1.15e-2],
                [3.0e-4, 0.0, 5.10e-3],
                [1.15e-2, 5.10e-3, 0.0],
            ]
        ),
    }
    return ePCSAFTMixture.from_params(params, species=["Methane", "Ethane", "Propane"])


def methanol_cyclohexane_mixture() -> ePCSAFTMixture:
    params = {
        "MW": np.asarray([32.042e-3, 84.147e-3]),
        "m": np.asarray([1.5255, 2.5303]),
        "s": np.asarray([3.2300, 3.8499]),
        "e": np.asarray([188.90, 278.11]),
        "e_assoc": np.asarray([2899.5, 0.0]),
        "vol_a": np.asarray([0.035176, 0.0]),
        "assoc_scheme": ["2B", None],
        "k_ij": np.asarray([[0.0, 0.051], [0.051, 0.0]]),
        "z": np.asarray([0.0, 0.0]),
        "dielc": np.asarray([33.05, 2.02]),
    }
    return ePCSAFTMixture.from_params(params, species=["Methanol", "Cyclohexane"])


def assert_plot_with_data(path: Path) -> None:
    csv_path = path.parent / "data" / f"{path.stem}_plot_data.csv"
    svg_path = path.with_suffix(".svg")
    html_path = path.with_suffix(".html")
    assert path.exists()
    assert svg_path.exists()
    assert html_path.exists()
    assert csv_path.exists()
    assert csv_path.stat().st_size > 0


def save_plotly_html(fig: go.Figure, image_path: Path) -> Path:
    html_path = plot_outputs.plot_html_path(image_path)
    html_path.parent.mkdir(parents=True, exist_ok=True)
    _translate_plotly_axis_titles(fig)
    fig.update_layout(
        template="plotly_white",
        colorway=MATPLOTLIB_COLORWAY,
        margin={"l": 72, "r": 30, "t": 74, "b": 88},
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
        height=390,
    )
    fig.update_xaxes(automargin=True, title_standoff=12)
    fig.update_yaxes(automargin=True, title_standoff=12)
    fig.write_html(
        html_path,
        include_plotlyjs="cdn",
        include_mathjax="cdn",
        full_html=True,
        config={"displaylogo": False, "responsive": True},
    )
    return html_path


_SPECIES_LABELS = {
    "A": r"A",
    "B": r"B",
    "C": r"C",
    "water": r"\mathrm{H_2O}",
    "Na": r"\mathrm{Na^+}",
    "Cl": r"\mathrm{Cl^-}",
    "methane": r"\mathrm{CH_4}",
    "ethane": r"\mathrm{C_2H_6}",
    "propane": r"\mathrm{C_3H_8}",
    "methanol": r"\mathrm{MeOH}",
    "cyclohexane": r"\mathrm{C_6H_{12}}",
}

_EXACT_MATH_LABELS = {
    "p": r"$P$",
    "P": r"$P$",
    "pressure": r"$P$",
    "diag pressure": r"diagnostics $P$",
    "rho": r"$\rho$",
    "rho_molar": r"$\rho$ (molar)",
    "rho_mass": r"$\rho_{\mathrm{mass}}$",
    "rho mass": r"$\rho_{\mathrm{mass}}$",
    "density": r"$\rho$",
    "z": r"$Z$",
    "Z": r"$Z$",
    "ares": r"$A^{res}$",
    "dadt": r"$\partial{}A^{res}/\partial{}T$",
    "hres": r"$h^{res}$",
    "sres": r"$s^{res}$",
    "gres": r"$g^{res}$",
    "epsr": r"$\epsilon_r$",
    "epsr mixture": r"$\epsilon_r$ mixture",
    "osmotic": r"$\phi_{\mathrm{osm}}$",
    "osmotic_coef": r"$\phi_{\mathrm{osm}}$",
    "material balance": "material balance",
    "material tolerance": "material balance tolerance",
    "fugacity residual": r"$\ln f_i$ residual",
    "fugacity tolerance": r"$\ln f_i$ tolerance",
    "mean gamma x": r"$\gamma_{\pm}^{(x)}$",
    "mean gamma m": r"$\gamma_{\pm}^{(m)}$",
    "ionic mean gamma mole": r"ionic $\gamma_{\pm}^{(x)}$",
    "ionic mean gamma molality": r"ionic $\gamma_{\pm}^{(m)}$",
    "stable min TPD": r"stable $\min(\mathrm{TPD})$",
    "unstable min TPD": r"unstable $\min(\mathrm{TPD})$",
    "sigma": r"$\sigma$",
    "epsilon": r"$\epsilon/k_B$",
    "m": r"$m$",
}

_PLOTLY_EXACT_LABELS = {
    "p": "P",
    "P": "P",
    "pressure": "P",
    "diag pressure": "diagnostics P",
    "rho": "rho (ρ)",
    "rho_molar": "rho (ρ, molar)",
    "rho_mass": "rho_mass (ρ_mass)",
    "rho mass": "rho_mass (ρ_mass)",
    "density": "rho (ρ)",
    "z": "Z",
    "Z": "Z",
    "ares": "A^res",
    "dadt": "dA^res/dT",
    "hres": "h^res",
    "sres": "s^res",
    "gres": "g^res",
    "epsr": "epsilon_r (ε_r)",
    "epsr mixture": "epsilon_r (ε_r) mixture",
    "osmotic": "phi_osm (φ_osm)",
    "osmotic_coef": "phi_osm (φ_osm)",
    "material balance": "material balance",
    "material tolerance": "material balance tolerance",
    "fugacity residual": "ln f_i residual",
    "fugacity tolerance": "ln f_i tolerance",
    "mean gamma x": "gamma_pm^(x)",
    "mean gamma m": "gamma_pm^(m)",
    "ionic mean gamma mole": "ionic gamma_pm^(x)",
    "ionic mean gamma molality": "ionic gamma_pm^(m)",
    "stable min TPD": "stable min(TPD)",
    "unstable min TPD": "unstable min(TPD)",
    "sigma": "sigma",
    "epsilon": "epsilon/kB",
    "m": "m",
}

_PLOTLY_SPECIES_LABELS = {
    "A": "A",
    "B": "B",
    "C": "C",
    "water": "H2O",
    "Na": "Na+",
    "Cl": "Cl-",
    "methane": "CH4",
    "ethane": "C2H6",
    "propane": "C3H8",
    "methanol": "MeOH",
    "cyclohexane": "C6H12",
}


def math_label(label: object) -> str:
    text = str(label)
    if text in _EXACT_MATH_LABELS:
        return _EXACT_MATH_LABELS[text]
    if text in _SPECIES_LABELS:
        return f"${_SPECIES_LABELS[text]}$"
    if text.lower() in _SPECIES_LABELS:
        return f"${_SPECIES_LABELS[text.lower()]}$"

    for prefix, rendered_prefix in (("mures", r"$\mu^{res}"), ("lnphi", r"$\ln \phi")):
        if text.startswith(f"{prefix}[") and text.endswith("]"):
            index = text[len(prefix) + 1 : -1]
            return f"{rendered_prefix}_{{{index}}}$"

    for prefix in ("neutral", "ionic", "vap", "liq", "stable", "unstable"):
        marker = f"{prefix} "
        if text.startswith(marker):
            return f"{prefix} {math_label(text[len(marker):])}"

    for prefix, rendered_prefix in (("beta", r"$\beta"), ("gamma", r"$\gamma"), ("mean gamma", r"$\gamma_{\pm}"), ("mures", r"$\mu^{res}"), ("lnphi", r"$\ln \phi"), ("phi", r"$\phi"), ("fugcoef", r"$\ln \phi"), ("gsolv", r"$\Delta G^{solv}")):
        marker = f"{prefix} "
        if text.startswith(marker):
            species = text[len(marker):]
            species_label = _SPECIES_LABELS.get(species, rf"\mathrm{{{species}}}")
            return f"{rendered_prefix}_{{{species_label}}}$"

    for suffix in ("total", "hc", "disp", "branch", "extreme"):
        marker = f" {suffix}"
        if text.endswith(marker):
            return f"{math_label(text[: -len(marker)])} {suffix}"

    if text.startswith("d") and "-d" in text:
        return rf"$\Delta(\partial{{}}A^{{res}}/\partial{{}}x)$ {text}"

    words = text.split()
    if len(words) > 1:
        return " ".join(
            math_label(word)
            if word in _EXACT_MATH_LABELS or word in _SPECIES_LABELS or word.lower() in _SPECIES_LABELS
            else word
            for word in words
        )
    return text


def math_labels(labels: Iterable[object]) -> list[str]:
    return [math_label(label) for label in labels]


def plotly_label(label: object) -> str:
    text = str(label)
    if text in _PLOTLY_EXACT_LABELS:
        return _PLOTLY_EXACT_LABELS[text]
    if text in _PLOTLY_SPECIES_LABELS:
        return _PLOTLY_SPECIES_LABELS[text]
    if text.lower() in _PLOTLY_SPECIES_LABELS:
        return _PLOTLY_SPECIES_LABELS[text.lower()]

    for prefix, rendered_prefix in (("mures", "mu^res"), ("lnphi", "ln phi")):
        if text.startswith(f"{prefix}[") and text.endswith("]"):
            index = text[len(prefix) + 1 : -1]
            return f"{rendered_prefix}_{index}"

    for prefix in ("neutral", "ionic", "vap", "liq", "stable", "unstable"):
        marker = f"{prefix} "
        if text.startswith(marker):
            return f"{prefix} {plotly_label(text[len(marker):])}"

    for prefix, rendered_prefix in (
        ("beta", "beta (β)"),
        ("gamma", "gamma (γ)"),
        ("mean gamma", "gamma_pm (γ_pm)"),
        ("mures", "mu^res (μ^res)"),
        ("lnphi", "ln phi (ln φ)"),
        ("phi", "phi (φ)"),
        ("fugcoef", "ln phi"),
        ("gsolv", "Delta G^solv (ΔG^solv)"),
    ):
        marker = f"{prefix} "
        if text.startswith(marker):
            species = text[len(marker):]
            species_label = _PLOTLY_SPECIES_LABELS.get(species, species)
            return f"{rendered_prefix}_{species_label}"

    for suffix in ("total", "hc", "disp", "branch", "extreme"):
        marker = f" {suffix}"
        if text.endswith(marker):
            return f"{plotly_label(text[: -len(marker)])} {suffix}"

    if text.startswith("d") and "-d" in text:
        return f"Delta(dA^res/dx) {text}"

    words = text.split()
    if len(words) > 1:
        return " ".join(
            plotly_label(word)
            if word in _PLOTLY_EXACT_LABELS
            or word in _PLOTLY_SPECIES_LABELS
            or word.lower() in _PLOTLY_SPECIES_LABELS
            else word
            for word in words
        )
    return text.replace("$", "")


def plotly_labels(labels: Iterable[object]) -> list[str]:
    return [plotly_label(label) for label in labels]


def plotly_axis_label(label: object) -> str:
    text = str(label)
    replacements = {
        r"$\rho$": "rho (ρ)",
        r"$A^{res}$": "A^res",
        r"$h^{res}$": "h^res",
        r"$s^{res}$": "s^res",
        r"$g^{res}$": "g^res",
        r"$Z$": "Z",
        r"$\ln \phi$": "ln φ",
        r"$\phi$": "φ",
        r"$\beta$": "β",
        r"$\gamma$": "γ",
        r"$x_i$": "x_i",
        r"$y_i$": "y_i",
        r"$\gamma_{\pm}$": "γ±",
        r"$\phi_{\mathrm{osm}}$": "φ_osm",
        r"$\Delta G^{solv}$": "ΔG^solv",
    }
    for source, replacement in replacements.items():
        text = text.replace(source, replacement)
    text = re.sub(r"\$x_\{\\mathrm\{([^}]+)\}\}\$", r"x_\1", text)
    text = re.sub(r"\$y_\{\\mathrm\{([^}]+)\}\}\$", r"y_\1", text)
    text = re.sub(r"\$\\ln\s+\\phi_\{\\mathrm\{([^}]+)\}\}\$", r"ln φ_\1", text)
    text = text.replace(r"\mathrm{", "").replace("}", "")
    return text.replace("$", "")


def _translate_plotly_axis_titles(fig: go.Figure) -> None:
    for axis_name in fig.layout:
        if not (str(axis_name).startswith("xaxis") or str(axis_name).startswith("yaxis")):
            continue
        axis = fig.layout[axis_name]
        if axis.title and axis.title.text:
            axis.title.text = plotly_axis_label(axis.title.text)


def _wrap_label(label: str, width: int = 18) -> str:
    return "\n".join(textwrap.wrap(str(label), width=width, break_long_words=False, break_on_hyphens=False)) or str(label)


def _wrapped_labels(labels: Iterable[str], *, width: int = 18) -> list[str]:
    return [_wrap_label(label, width=width) for label in labels]


def _max_wrapped_lines(labels: Iterable[str], *, width: int = 18) -> int:
    wrapped = _wrapped_labels(labels, width=width)
    return max((label.count("\n") + 1 for label in wrapped), default=1)


def _comparison_size(labels: list[str], *, relative_error: bool) -> tuple[float, float]:
    label_count = max(len(labels), 1)
    longest = max((len(label) for label in labels), default=1)
    width = max(7.5, min(18.0, label_count * 0.92 + longest * 0.08))
    height = (5.8 if relative_error else 4.5) + min(2.2, 0.22 * _max_wrapped_lines(labels))
    return width, height


def _maybe_use_symmetric_log_scale(ax, values: np.ndarray) -> None:
    finite = np.abs(values[np.isfinite(values)])
    nonzero = finite[finite > 0.0]
    if nonzero.size == 0:
        return
    if float(np.max(nonzero) / np.min(nonzero)) <= 1.0e4:
        return
    linthresh = max(float(np.min(nonzero)) * 0.5, 1.0e-12)
    ax.set_yscale("symlog", linthresh=linthresh)
    ax.text(
        0.995,
        0.04,
        "symlog scale",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize="x-small",
        color="0.35",
    )


def _maybe_use_symmetric_log_x_scale(ax, values: np.ndarray) -> None:
    finite = np.abs(values[np.isfinite(values)])
    nonzero = finite[finite > 0.0]
    if nonzero.size == 0:
        return
    if float(np.max(nonzero) / np.min(nonzero)) <= 1.0e4:
        return
    linthresh = max(float(np.min(nonzero)) * 0.5, 1.0e-12)
    ax.set_xscale("symlog", linthresh=linthresh)
    ax.text(
        0.98,
        0.04,
        "symlog scale",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize="x-small",
        color="0.35",
    )


def _finish_figure(fig) -> None:
    fig.align_labels()
    fig.tight_layout(pad=1.35)


def assert_figure_text_is_inside_canvas(fig, *, margin_px: float = 1.0) -> None:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    canvas = fig.bbox
    for text in fig.findobj(match=lambda artist: hasattr(artist, "get_window_extent") and hasattr(artist, "get_text")):
        if not text.get_visible() or not text.get_text():
            continue
        bbox = text.get_window_extent(renderer=renderer)
        assert bbox.x0 >= canvas.x0 - margin_px
        assert bbox.y0 >= canvas.y0 - margin_px
        assert bbox.x1 <= canvas.x1 + margin_px
        assert bbox.y1 <= canvas.y1 + margin_px


def _plotly_axis_type(values: np.ndarray) -> str:
    finite = np.abs(values[np.isfinite(values)])
    nonzero = finite[finite > 0.0]
    if nonzero.size == 0:
        return "linear"
    return "log" if float(np.max(nonzero) / np.min(nonzero)) > 1.0e4 and np.all(values > 0.0) else "linear"


def _save_interactive_comparison_plot(
    image_path: Path,
    title: str,
    labels: list[str],
    actual: np.ndarray,
    expected: np.ndarray,
    *,
    ylabel: str,
    relative_error: bool,
) -> Path:
    display_labels = plotly_labels(labels)
    if relative_error:
        fig = make_subplots(
            rows=1,
            cols=2,
            subplot_titles=(title, "Relative error"),
            column_widths=[0.72, 0.28],
        )
    else:
        fig = go.Figure()

    values = np.concatenate([actual, expected])
    scale = np.maximum(np.abs(expected), 1.0e-30)
    rel_error = (actual - expected) / scale
    orientation = "h" if len(display_labels) > 8 or max((len(label) for label in display_labels), default=0) > 18 else "v"

    if orientation == "h":
        fig.add_trace(
            go.Bar(y=display_labels, x=actual, name="Actual model output", orientation="h"),
            row=1 if relative_error else None,
            col=1 if relative_error else None,
        )
        fig.add_trace(
            go.Bar(y=display_labels, x=expected, name="Expected/reference", orientation="h"),
            row=1 if relative_error else None,
            col=1 if relative_error else None,
        )
        if relative_error:
            fig.add_trace(go.Bar(y=display_labels, x=rel_error, name="Relative error", orientation="h"), row=1, col=2)
            fig.update_yaxes(autorange="reversed", row=1, col=1)
            fig.update_yaxes(autorange="reversed", showticklabels=False, row=1, col=2)
            fig.update_xaxes(title_text=ylabel, type=_plotly_axis_type(values[values > 0.0]), row=1, col=1)
            fig.update_xaxes(title_text="Relative error", row=1, col=2)
        else:
            fig.update_yaxes(autorange="reversed")
            fig.update_xaxes(title_text=ylabel, type=_plotly_axis_type(values[values > 0.0]))
    else:
        fig.add_trace(
            go.Bar(x=display_labels, y=actual, name="Actual model output"),
            row=1 if relative_error else None,
            col=1 if relative_error else None,
        )
        fig.add_trace(
            go.Bar(x=display_labels, y=expected, name="Expected/reference"),
            row=1 if relative_error else None,
            col=1 if relative_error else None,
        )
        if relative_error:
            fig.add_trace(go.Bar(x=display_labels, y=rel_error, name="Relative error"), row=1, col=2)
            fig.update_yaxes(title_text=ylabel, type=_plotly_axis_type(values[values > 0.0]), row=1, col=1)
            fig.update_yaxes(title_text="Relative error", row=1, col=2)
        else:
            fig.update_yaxes(title_text=ylabel, type=_plotly_axis_type(values[values > 0.0]))

    fig.update_layout(
        title=title if not relative_error else None,
        barmode="group",
    )
    return save_plotly_html(fig, image_path)


def _save_interactive_parity_plot(
    image_path: Path,
    title: str,
    labels: list[str],
    actual: np.ndarray,
    expected: np.ndarray,
    *,
    xlabel: str,
    ylabel: str,
) -> Path:
    display_labels = plotly_labels(labels)
    finite = np.isfinite(actual) & np.isfinite(expected)
    scale = np.maximum(np.abs(expected), 1.0e-30)
    rel_error = (actual - expected) / scale
    fig = make_subplots(rows=1, cols=2, subplot_titles=(title, "Relative error"), column_widths=[0.62, 0.38])
    fig.add_trace(
        go.Scatter(
            x=expected[finite],
            y=actual[finite],
            mode="markers",
            text=np.asarray(display_labels, dtype=object)[finite],
            name="Actual vs expected points",
            marker={"color": MATPLOTLIB_COLORWAY[0]},
            hovertemplate="%{text}<br>expected=%{x:.6g}<br>actual=%{y:.6g}<extra></extra>",
        ),
        row=1,
        col=1,
    )
    if np.any(finite):
        lo = float(min(np.min(expected[finite]), np.min(actual[finite])))
        hi = float(max(np.max(expected[finite]), np.max(actual[finite])))
        pad = max((hi - lo) * 0.08, 1.0e-12)
        fig.add_trace(
            go.Scatter(
                x=[lo - pad, hi + pad],
                y=[lo - pad, hi + pad],
                mode="lines",
                name="Parity line (actual = expected)",
                line={"color": "#404040", "width": 1.2},
            ),
            row=1,
            col=1,
        )
    fig.add_trace(go.Bar(y=display_labels, x=rel_error, orientation="h", name="Relative error"), row=1, col=2)
    fig.update_xaxes(title_text=xlabel, row=1, col=1)
    fig.update_yaxes(title_text=ylabel, row=1, col=1)
    fig.update_xaxes(title_text="Relative error", row=1, col=2)
    fig.update_yaxes(autorange="reversed", row=1, col=2)
    return save_plotly_html(fig, image_path)


def _save_interactive_contribution_plot(
    image_path: Path,
    title: str,
    rows: list[dict[str, object]],
    *,
    breakdown: bool,
) -> Path:
    labels = plotly_labels(str(row["label"]) for row in rows)
    term_names = sorted({name for row in rows for name in row["terms"]})
    totals = np.asarray([float(row["total"]) for row in rows], dtype=float)
    fig = go.Figure()
    for term_name in term_names:
        values = [float(row["terms"].get(term_name, 0.0)) for row in rows]
        fig.add_trace(go.Bar(y=labels, x=values, orientation="h", name=f"Term: {term_name}"))
    fig.add_trace(
        go.Scatter(
            y=labels,
            x=totals,
            mode="markers",
            marker_symbol="x",
            marker_size=10,
            marker={"color": "#111111"},
            name="Reported total",
        )
    )
    fig.update_layout(
        title=title,
        barmode="relative" if not breakdown else "group",
        xaxis_title="Contribution value" if not breakdown else "Term contribution",
        yaxis={"autorange": "reversed"},
    )
    return save_plotly_html(fig, image_path)


def save_comparison_plot(
    filename: str,
    title: str,
    labels: list[str],
    actual: np.ndarray,
    expected: np.ndarray,
    *,
    category: Iterable[str],
    ylabel: str = "Value",
    relative_error: bool = True,
) -> Path:
    actual = np.asarray(actual, dtype=float)
    expected = np.asarray(expected, dtype=float)
    display_labels = math_labels(labels)
    is_dense = len(display_labels) > 12 or max((len(label) for label in display_labels), default=0) > 22
    if is_dense:
        wrapped_labels = _wrapped_labels(display_labels, width=24)
        y = np.arange(len(labels), dtype=float)
        fig_height = max(6.8, min(13.5, 1.8 + 0.38 * len(labels)))
        if relative_error:
            fig, axes = plt.subplots(1, 2, figsize=(14.8, fig_height), width_ratios=[2.3, 1.05])
            ax, err_ax = axes
        else:
            fig, ax = plt.subplots(figsize=(10.5, fig_height))
            err_ax = None

        ax.barh(y - 0.18, actual, height=0.36, label="Actual")
        ax.barh(y + 0.18, expected, height=0.36, label="Expected")
        ax.set_title(title)
        ax.set_xlabel(ylabel)
        ax.set_yticks(y)
        ax.set_yticklabels(wrapped_labels, fontsize="small")
        ax.invert_yaxis()
        ax.grid(axis="x", color="0.88", linewidth=0.7)
        ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2, frameon=False)
        ax.axvline(0.0, color="0.82", linewidth=0.8)
        _maybe_use_symmetric_log_x_scale(ax, np.concatenate([actual, expected]))

        if err_ax is not None:
            scale = np.maximum(np.abs(expected), 1.0e-30)
            err_ax.barh(y, (actual - expected) / scale, height=0.55, label="Relative error")
            err_ax.axvline(0.0, color="0.25", linewidth=0.8)
            err_ax.set_xlabel("Rel. err.")
            err_ax.set_yticks(y)
            err_ax.set_yticklabels([])
            err_ax.invert_yaxis()
            err_ax.grid(axis="x", color="0.9", linewidth=0.7)
    else:
        x = np.arange(len(display_labels), dtype=float)
        wrapped_labels = _wrapped_labels(display_labels)
        figsize = _comparison_size(display_labels, relative_error=relative_error)
        if relative_error:
            fig, axes = plt.subplots(2, 1, figsize=figsize, height_ratios=[3, 1.45])
            ax, err_ax = axes
        else:
            fig, ax = plt.subplots(figsize=figsize)
            err_ax = None

        ax.bar(x - 0.18, actual, width=0.36, label="Actual")
        ax.bar(x + 0.18, expected, width=0.36, label="Expected")
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.set_xticks(x)
        ax.set_xticklabels([] if err_ax is not None else wrapped_labels, rotation=0, ha="center")
        ax.tick_params(axis="x", length=0 if err_ax is not None else 3)
        ax.grid(axis="y", color="0.88", linewidth=0.7)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.14), ncol=2, frameon=False)
        ax.axhline(0.0, color="0.82", linewidth=0.8)
        _maybe_use_symmetric_log_scale(ax, np.concatenate([actual, expected]))

        if err_ax is not None:
            scale = np.maximum(np.abs(expected), 1.0e-30)
            err_ax.bar(x, (actual - expected) / scale, width=0.5, label="Relative error")
            err_ax.axhline(0.0, color="0.25", linewidth=0.8)
            err_ax.set_ylabel("Rel. err.")
            err_ax.grid(axis="y", color="0.9", linewidth=0.7)
            err_ax.set_xticks(x)
            err_ax.set_xticklabels(wrapped_labels, rotation=0, ha="center", fontsize="small")

    output_path = plot_outputs.test_plot_path(__file__, filename, category=category)
    try:
        _finish_figure(fig)
        plot_outputs.save_plot_figure(fig, output_path, dpi=120, svg_companion=True)
        _save_interactive_comparison_plot(
            output_path,
            title,
            labels,
            actual,
            expected,
            ylabel=ylabel,
            relative_error=relative_error,
        )
    finally:
        plt.close(fig)
    assert_plot_with_data(output_path)
    return output_path


def save_parity_plot(
    filename: str,
    title: str,
    labels: list[str],
    actual: np.ndarray,
    expected: np.ndarray,
    *,
    category: Iterable[str],
    xlabel: str = "Expected",
    ylabel: str = "Actual",
) -> Path:
    actual = np.asarray(actual, dtype=float)
    expected = np.asarray(expected, dtype=float)
    finite = np.isfinite(actual) & np.isfinite(expected)
    plot_actual = actual[finite]
    plot_expected = expected[finite]

    display_labels = math_labels(labels)
    wrapped_labels = _wrapped_labels(display_labels, width=22)
    fig_height = max(5.2, min(11.0, 1.6 + 0.28 * len(labels)))
    fig, axes = plt.subplots(1, 2, figsize=(13.2, fig_height), width_ratios=[2.3, 1.9])
    ax, err_ax = axes
    ax.scatter(plot_expected, plot_actual, label="Comparison points")
    if plot_expected.size:
        lo = float(min(np.min(plot_expected), np.min(plot_actual)))
        hi = float(max(np.max(plot_expected), np.max(plot_actual)))
        pad = max((hi - lo) * 0.08, 1.0e-12)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color="0.25", linewidth=1.0, label="Parity")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(color="0.9", linewidth=0.7)
    ax.legend(loc="best", frameon=False)

    scale = np.maximum(np.abs(expected), 1.0e-30)
    y = np.arange(len(display_labels), dtype=float)
    err_ax.barh(y, (actual - expected) / scale, height=0.55)
    err_ax.axvline(0.0, color="0.25", linewidth=0.8)
    err_ax.set_title("Relative error")
    err_ax.set_xlabel("Relative error")
    err_ax.set_yticks(y)
    err_ax.set_yticklabels(wrapped_labels, fontsize="small")
    err_ax.invert_yaxis()
    err_ax.grid(axis="x", color="0.9", linewidth=0.7)

    output_path = plot_outputs.test_plot_path(__file__, filename, category=category)
    try:
        _finish_figure(fig)
        plot_outputs.save_plot_figure(fig, output_path, dpi=120, svg_companion=True)
        _save_interactive_parity_plot(output_path, title, labels, actual, expected, xlabel=xlabel, ylabel=ylabel)
    finally:
        plt.close(fig)
    assert_plot_with_data(output_path)
    return output_path


def append_payload_rows(rows: list[dict[str, object]], prefix: str, payload: dict[str, object]) -> None:
    total = np.asarray(payload["total"], dtype=float)
    term_arrays = {key: np.asarray(value, dtype=float) for key, value in payload["terms"].items()}
    if total.ndim == 0:
        rows.append(
            {
                "label": prefix,
                "terms": {key: float(value) for key, value in payload["terms"].items()},
                "total": float(total),
            }
        )
        return
    for index, total_value in enumerate(total.tolist()):
        rows.append(
            {
                "label": f"{prefix}[{index}]",
                "terms": {key: float(value[index]) for key, value in term_arrays.items()},
                "total": float(total_value),
            }
        )


def save_contribution_closure_plot(
    filename: str,
    title: str,
    rows: list[dict[str, object]],
    *,
    category: Iterable[str],
) -> Path:
    labels = math_labels(str(row["label"]) for row in rows)
    wrapped_labels = _wrapped_labels(labels, width=18)
    term_names = sorted({name for row in rows for name in row["terms"]})
    y = np.arange(len(rows), dtype=float)

    fig_height = max(7.0, min(12.0, 2.2 + 0.48 * len(rows)))
    fig, axes = plt.subplots(1, 2, figsize=(13.2, fig_height), width_ratios=[2.4, 1.25])
    ax, err_ax = axes
    positive_left = np.zeros(len(rows), dtype=float)
    negative_left = np.zeros(len(rows), dtype=float)
    for term_name in term_names:
        values = np.asarray([float(row["terms"].get(term_name, 0.0)) for row in rows], dtype=float)
        lefts = np.where(values >= 0.0, positive_left, negative_left)
        ax.barh(y, values, left=lefts, height=0.58, label=term_name)
        positive_left += np.where(values >= 0.0, values, 0.0)
        negative_left += np.where(values < 0.0, values, 0.0)

    totals = np.asarray([float(row["total"]) for row in rows], dtype=float)
    term_sums = np.asarray([sum(float(value) for value in row["terms"].values()) for row in rows], dtype=float)
    ax.scatter(totals, y, marker="x", color="black", label="Reported total", zorder=5)
    ax.set_title(title)
    ax.set_xlabel("Contribution value")
    ax.set_yticks(y)
    ax.set_yticklabels(wrapped_labels, fontsize="small")
    ax.invert_yaxis()
    ax.axvline(0.0, color="0.82", linewidth=0.8)
    ax.grid(axis="x", color="0.9", linewidth=0.7)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=min(4, max(1, len(term_names))), fontsize="small", frameon=False)

    scale = np.maximum(np.abs(totals), 1.0e-30)
    err_ax.barh(y, (term_sums - totals) / scale, height=0.55, label="Closure error")
    err_ax.axvline(0.0, color="0.25", linewidth=0.8)
    err_ax.set_xlabel("Rel. closure")
    err_ax.set_yticks(y)
    err_ax.set_yticklabels([])
    err_ax.invert_yaxis()
    err_ax.grid(axis="x", color="0.9", linewidth=0.7)

    output_path = plot_outputs.test_plot_path(__file__, filename, category=category)
    try:
        _finish_figure(fig)
        plot_outputs.save_plot_figure(fig, output_path, dpi=120, svg_companion=True)
        _save_interactive_contribution_plot(output_path, title, rows, breakdown=False)
    finally:
        plt.close(fig)
    assert_plot_with_data(output_path)
    return output_path


def save_contribution_term_breakdown_plot(
    filename: str,
    title: str,
    rows: list[dict[str, object]],
    *,
    category: Iterable[str],
) -> Path:
    labels = math_labels(str(row["label"]) for row in rows)
    wrapped_labels = _wrapped_labels(labels, width=18)
    term_names = sorted({name for row in rows for name in row["terms"]})
    y = np.arange(len(rows), dtype=float)
    height = min(0.75 / max(len(term_names), 1), 0.18)
    offsets = (np.arange(len(term_names), dtype=float) - (len(term_names) - 1) / 2.0) * height

    fig_height = max(6.4, min(12.0, 1.8 + 0.46 * len(rows)))
    fig, ax = plt.subplots(figsize=(12.4, fig_height))
    for offset, term_name in zip(offsets, term_names, strict=False):
        values = np.asarray([float(row["terms"].get(term_name, 0.0)) for row in rows], dtype=float)
        ax.barh(y + offset, values, height=height, label=term_name)

    ax.scatter([float(row["total"]) for row in rows], y, marker="x", color="black", label="Reported total", zorder=5)
    ax.set_title(title)
    ax.set_xlabel("Term contribution")
    ax.set_yticks(y)
    ax.set_yticklabels(wrapped_labels, fontsize="small")
    ax.invert_yaxis()
    ax.axvline(0.0, color="0.82", linewidth=0.8)
    ax.grid(axis="x", color="0.9", linewidth=0.7)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=min(4, max(1, len(term_names))), fontsize="small", frameon=False)

    output_path = plot_outputs.test_plot_path(__file__, filename, category=category)
    try:
        _finish_figure(fig)
        plot_outputs.save_plot_figure(fig, output_path, dpi=120, svg_companion=True)
        _save_interactive_contribution_plot(output_path, title, rows, breakdown=True)
    finally:
        plt.close(fig)
    assert_plot_with_data(output_path)
    return output_path


def finite_difference_gradient_values() -> tuple[np.ndarray, np.ndarray]:
    theta = {"m": 1.05, "s": 3.68, "e": 151.0}

    def objective_at(m: float, s: float, e: float) -> float:
        debug = _debug_native_pure_neutral_objective(
            _methane_like_records(),
            "Methane",
            assoc_scheme="",
            fixed_parameters=_minimal_neutral_metadata(16.043e-3),
            initial_guess=theta,
            x={"m": m, "s": s, "e": e},
        )
        return float(debug["objective"])

    debug = _debug_native_pure_neutral_objective(
        _methane_like_records(),
        "Methane",
        assoc_scheme="",
        fixed_parameters=_minimal_neutral_metadata(16.043e-3),
        initial_guess=theta,
        x=theta,
    )
    exact = np.asarray(debug["gradient"], dtype=float)
    eps = np.asarray([1.0e-6, 1.0e-6, 1.0e-5], dtype=float)
    fd = np.empty(3, dtype=float)
    base = np.asarray([theta["m"], theta["s"], theta["e"]], dtype=float)
    for i in range(3):
        forward = base.copy()
        backward = base.copy()
        forward[i] += eps[i]
        backward[i] -= eps[i]
        fd[i] = (objective_at(*forward) - objective_at(*backward)) / (2.0 * eps[i])
    return exact, fd
