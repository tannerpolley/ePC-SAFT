"""Generate Figure 6b fit-comparison plots on the paper-style ``mu`` basis.

Each contribution plot overlays the digitized paper points against the raw
residual chemical-potential route used by the paper-style decomposition.

The total plot additionally shows:
- the digitized paper total points
- a connected paper-total guide line through those points
- the EOS total from ``lnfugcoef_total``
- the sum of the ``mu`` contribution curves
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Dict, Tuple

import matplotlib
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[6]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from scripts.plot_outputs import paper_validation_path
from figure6b_libr_ethanol_contributions import (
    AXIS_LABEL_SIZE,
    AXIS_TICK_SIZE,
    _build_params,
    _calc_ln_miac_contributions,
    _load_exp_data,
    _salt_mole_fraction_from_molality,
)
from figure6b_digitized_reference_replica import (
    SERIES_STYLES,
    _load_digitized_curves,
)

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTPUT_ROOT = paper_validation_path(Path(__file__).resolve().parent, "output")
OUTPUT_DATA_DIR = OUTPUT_ROOT / "data"
OUTPUT_PLOTS_DIR = OUTPUT_ROOT / "plots"

CONTRIBUTION_ORDER = ["born", "dh", "hc", "disp", "assoc", "total"]
FILE_LABELS = {
    "born": "born",
    "dh": "dh",
    "hc": "hard_chain",
    "disp": "dispersion",
    "assoc": "association",
    "total": "total",
}
TITLE_LABELS = {
    "born": "Born",
    "dh": "Debye-Huckel",
    "hc": "Hard-chain",
    "disp": "Dispersion",
    "assoc": "Association",
    "total": "Total",
}
METHOD_STYLES = {
    "mu": {
        "label": r"ePC-SAFT contribution from $\Delta \tilde{\mu}^{(\alpha)}$",
        "color": "#1f77b4",
        "linestyle": "--",
        "linewidth": 1.8,
    },
}


def _metric_summary(x_data: np.ndarray, y_data: np.ndarray, x_model: np.ndarray, y_model: np.ndarray) -> Tuple[float, float, float]:
    y_interp = np.interp(x_data, x_model, y_model)
    delta = y_interp - y_data
    rmse = float(np.sqrt(np.mean(delta * delta)))
    mae = float(np.mean(np.abs(delta)))
    max_abs = float(np.max(np.abs(delta)))
    return rmse, mae, max_abs


def _y_limits(y_data: np.ndarray, y_model: np.ndarray) -> Tuple[float, float]:
    y_all = np.concatenate([np.asarray(y_data, dtype=float), np.asarray(y_model, dtype=float)])
    y_min = float(np.min(y_all))
    y_max = float(np.max(y_all))
    span = max(y_max - y_min, 0.2)
    pad = 0.08 * span
    return y_min - pad, y_max + pad


def _sum_contribution_curves(curves: Dict[str, np.ndarray]) -> np.ndarray:
    return (
        np.asarray(curves["born"], dtype=float)
        + np.asarray(curves["dh"], dtype=float)
        + np.asarray(curves["hc"], dtype=float)
        + np.asarray(curves["disp"], dtype=float)
        + np.asarray(curves["assoc"], dtype=float)
    )


def run_analysis(
    miac_data_path: Path,
    digitized_path: Path,
    output_dir: Path,
    metrics_path: Path,
    grid_points: int,
    x_min: float,
    x_max: float,
    d_born_mode: int | None = None,
) -> Dict[str, Dict[str, float]]:
    m_exp, _, _ = _load_exp_data(miac_data_path)
    m_upper = float(np.max(m_exp))
    m_grid = np.linspace(0.0, m_upper, int(grid_points))
    x_grid = _salt_mole_fraction_from_molality(m_grid)

    user_options = {}
    if d_born_mode is not None:
        user_options = {"elec_model": {"born_model": {"d_Born_mode": int(d_born_mode)}}}
    params = _build_params(user_options=user_options)
    model_curves_by_method: Dict[str, Dict[str, np.ndarray]] = {}
    method_errors: Dict[str, str] = {}
    for method_name in ("lnphi", "mu"):
        try:
            model_curves_by_method[method_name] = _calc_ln_miac_contributions(m_grid, params, method=method_name)
        except Exception as exc:
            method_errors[method_name] = str(exc)
    born_compare_curves: Dict[str, np.ndarray] = {}
    if d_born_mode is None or int(d_born_mode) == 1:
        try:
            compare_params = _build_params(user_options={"elec_model": {"born_model": {"d_Born_mode": 0}}})
            born_compare_curves["dBorn0_mu"] = _calc_ln_miac_contributions(m_grid, compare_params, method="mu")["born"]
        except Exception as exc:
            method_errors["dBorn0_mu"] = str(exc)
    digitized = _load_digitized_curves(digitized_path)

    output_dir.mkdir(parents=True, exist_ok=True)
    results: Dict[str, Dict[str, float]] = {}
    summary_rows = []

    for name in CONTRIBUTION_ORDER:
        x_data, y_data = digitized[name]
        if name == "total":
            method_names = tuple(
                method_name for method_name in ("paper_total_line", "model_total", "mu_sum")
                if (
                    method_name == "paper_total_line"
                ) or (
                    method_name == "model_total"
                    and "lnphi" in model_curves_by_method
                ) or (
                    method_name == "mu_sum"
                    and "mu" in model_curves_by_method
                )
            )
        else:
            method_names = tuple(
                method_name for method_name in ("mu",)
                if method_name in model_curves_by_method
            )
        results[name] = {"n_points": float(len(x_data))}
        metric_lines = []

        style = SERIES_STYLES[name]
        fig, ax = plt.subplots(figsize=(7.6, 5.2))
        fig.patch.set_facecolor("white")
        ax.set_facecolor("white")

        ax.scatter(
            x_data,
            y_data,
            color="black",
            marker="o",
            s=34,
            facecolors="black",
            linewidths=0.8,
            label="Digitized paper data",
            zorder=7,
        )

        y_candidates = [np.asarray(y_data, dtype=float)]
        for method_name in method_names:
            if method_name == "paper_total_line":
                y_model = np.interp(x_grid, x_data, y_data)
                label = "Digitized paper total guide"
                color = "black"
                linestyle = "--"
                linewidth = 1.4
            elif method_name == "model_total":
                y_model = np.asarray(model_curves_by_method["lnphi"]["total"], dtype=float)
                label = r"ePC-SAFT EOS total from $\Delta \ln\varphi$"
                color = "#111111"
                linestyle = "-"
                linewidth = 2.2
            elif method_name == "mu_sum":
                y_model = _sum_contribution_curves(model_curves_by_method["mu"])
                label = r"Sum of ePC-SAFT $\Delta \tilde{\mu}^{(\alpha)}$ contributions"
                color = "#1f77b4"
                linestyle = ":"
                linewidth = 2.1
            else:
                y_model = np.asarray(model_curves_by_method[method_name][name], dtype=float)
                method_style = METHOD_STYLES[method_name]
                color = str(style["color"])
                linestyle = str(method_style["linestyle"])
                linewidth = float(method_style["linewidth"])
                label = str(method_style["label"])
            if method_name != "paper_total_line":
                rmse, mae, max_abs = _metric_summary(x_data, y_data, x_grid, y_model)
                results[name][f"{method_name}_rmse"] = rmse
                results[name][f"{method_name}_mae"] = mae
                results[name][f"{method_name}_max_abs"] = max_abs
                summary_rows.append(
                    {
                        "contribution": name,
                        "method": method_name,
                        "rmse": rmse,
                        "mae": mae,
                        "max_abs": max_abs,
                        "n_points": len(x_data),
                    }
                )
                metric_lines.append(f"{method_name}: RMSE={rmse:.4f}, MAE={mae:.4f}, Max |Δ|={max_abs:.4f}")
            ax.plot(
                x_grid,
                y_model,
                color=color,
                linestyle=linestyle,
                linewidth=linewidth,
                label=label,
                zorder=5,
            )
            y_candidates.append(y_model)

        if name == "born" and "dBorn0_mu" in born_compare_curves:
            y_model = np.asarray(born_compare_curves["dBorn0_mu"], dtype=float)
            rmse, mae, max_abs = _metric_summary(x_data, y_data, x_grid, y_model)
            results[name]["dBorn0_mu_rmse"] = rmse
            results[name]["dBorn0_mu_mae"] = mae
            results[name]["dBorn0_mu_max_abs"] = max_abs
            summary_rows.append(
                {
                    "contribution": name,
                    "method": "dBorn0_mu",
                    "rmse": rmse,
                    "mae": mae,
                    "max_abs": max_abs,
                    "n_points": len(x_data),
                }
            )
            metric_lines.append(f"dBorn0_mu: RMSE={rmse:.4f}, MAE={mae:.4f}, Max |Δ|={max_abs:.4f}")
            ax.plot(
                x_grid,
                y_model,
                color="#cc5500",
                linestyle="-",
                linewidth=1.9,
                label=r"ePC-SAFT Born with $d_{\mathrm{Born}}$ mode 0",
                zorder=4,
            )
            y_candidates.append(y_model)

        y_lo, y_hi = _y_limits(y_data, np.concatenate(y_candidates[1:]))
        ax.set_xlim(float(x_min), float(x_max))
        ax.set_ylim(y_lo, y_hi)
        ax.set_xlabel(r"salt mole fraction, $x_{salt}$", fontsize=AXIS_LABEL_SIZE)
        if name == "total":
            ax.set_ylabel(r"Total or summed contribution to $\ln(\gamma_{\pm}^{*})$", fontsize=AXIS_LABEL_SIZE)
            ax.set_title(r"Figure 6b fit check: paper total vs EOS total vs summed $\mu$ contributions")
        else:
            ax.set_ylabel(r"Contribution to $\ln(\gamma_{\pm}^{*})$", fontsize=AXIS_LABEL_SIZE)
            ax.set_title(f"Figure 6b fit check: {TITLE_LABELS[name]} ($\\mu$ basis)")
        ax.grid(True, alpha=0.3, color="0.7")
        ax.tick_params(colors="black", labelsize=AXIS_TICK_SIZE)
        for spine in ax.spines.values():
            spine.set_color("black")
            spine.set_linewidth(1.0)

        metric_text = "\n".join(metric_lines)
        ax.text(
            0.98,
            0.97,
            metric_text,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
            bbox={"facecolor": "white", "edgecolor": "black", "alpha": 1.0, "boxstyle": "round,pad=0.25"},
        )

        legend = ax.legend(fontsize=9)
        frame = legend.get_frame()
        frame.set_facecolor("white")
        frame.set_edgecolor("black")
        frame.set_alpha(1.0)

        out_path = output_dir / f"figure6b_fit_{FILE_LABELS[name]}.png"
        fig.tight_layout()
        fig.savefig(out_path, dpi=220)
        plt.close(fig)
        ranked_methods = tuple(method_name for method_name in method_names if method_name != "paper_total_line")
        best_method = min(ranked_methods, key=lambda method_name: float(results[name][f"{method_name}_rmse"]))
        print(f"{name}: best_method={best_method}, wrote={out_path}")

    for method_name, error_text in method_errors.items():
        print(f"skipped method {method_name}: {error_text}")

    metrics_path.parent.mkdir(parents=True, exist_ok=True)
    with metrics_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["contribution", "method", "rmse", "mae", "max_abs", "n_points"],
        )
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"wrote metrics: {metrics_path}")

    return results


def _parse_args() -> argparse.Namespace:
    default_digitized = (
        REPO_ROOT
        / "scripts"
        / "paper_validation"
        / "2020_Bulow_analysis"
        / "figure_6"
        / "figure_6b"
        / "data"
        / "Figure6b_curves.csv"
    )
    parser = argparse.ArgumentParser(description="Generate per-contribution Figure 6b model-vs-digitized fit plots")
    parser.add_argument(
        "--miac-data",
        type=Path,
        default=REPO_ROOT / "data" / "MIAC" / "ethanol" / "ethanol-LiBr.csv",
    )
    parser.add_argument(
        "--digitized",
        type=Path,
        default=default_digitized,
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=OUTPUT_PLOTS_DIR / "figure6b_fit_checks",
    )
    parser.add_argument(
        "--metrics-csv",
        type=Path,
        default=OUTPUT_DATA_DIR / "figure6b_fit_checks" / "figure6b_fit_method_metrics.csv",
    )
    parser.add_argument("--grid-points", type=int, default=1201)
    parser.add_argument("--x-min", type=float, default=0.0)
    parser.add_argument("--x-max", type=float, default=0.2)
    parser.add_argument(
        "--d-born-mode",
        type=int,
        choices=[0, 1],
        default=None,
        help="Optional override for 2020_Bulow born_model.d_Born_mode.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    run_analysis(
        miac_data_path=Path(args.miac_data),
        digitized_path=Path(args.digitized),
        output_dir=Path(args.out_dir),
        metrics_path=Path(args.metrics_csv),
        grid_points=int(args.grid_points),
        x_min=float(args.x_min),
        x_max=float(args.x_max),
        d_born_mode=None if args.d_born_mode is None else int(args.d_born_mode),
    )


if __name__ == "__main__":
    main()



