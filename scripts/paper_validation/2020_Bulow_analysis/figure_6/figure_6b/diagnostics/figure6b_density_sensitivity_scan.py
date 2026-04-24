"""Scan Figure 6b contribution sensitivity to imposed liquid-density scaling."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[6]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts._env import require_epcsaft_install
from scripts.plot_outputs import paper_validation_path, save_plot_figure

require_epcsaft_install()
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from figure6b_digitized_reference_replica import _load_digitized_curves
from figure6b_libr_ethanol_contributions import P_REF, T_REF, _build_params, _molality_for_salt_mole_fraction, _molality_to_species_molefraction
from scripts._epcsaft_oop import epcsaft_density, epcsaft_fugacity_coefficient_terms

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTPUT_ROOT = paper_validation_path(Path(__file__).resolve().parent, "output")
OUTPUT_DATA_DIR = OUTPUT_ROOT / "data"
OUTPUT_PLOTS_DIR = OUTPUT_ROOT / "plots"


CONTRIBUTION_ORDER = ["born", "dh", "hc", "disp", "assoc", "total"]
TITLE_LABELS = {
    "born": "Born",
    "dh": "Debye-Huckel",
    "hc": "Hard-chain",
    "disp": "Dispersion",
    "assoc": "Association",
    "total": "Total",
}


def _calc_point_contributions(x_salt: float, params: Dict[str, object], rho_scale: float) -> Dict[str, float]:
    molality = _molality_for_salt_mole_fraction(x_salt)
    m_eval = molality if molality > 0.0 else 1e-12
    x = _molality_to_species_molefraction(m_eval)
    rho_base = float(epcsaft_density(T_REF, P_REF, x, params, phase="liq"))
    rho = float(rho_scale) * rho_base

    terms = epcsaft_fugacity_coefficient_terms(T_REF, rho, x, params)
    eps = 1e-12
    x_inf = np.full_like(x, eps)
    x_inf[2] = max(1.0 - eps * (len(x) - 1), eps)
    x_inf /= np.sum(x_inf)
    rho_inf = float(epcsaft_density(T_REF, P_REF, x_inf, params, phase="liq"))
    terms_inf = epcsaft_fugacity_coefficient_terms(T_REF, rho_inf, x_inf, params)

    def mean_ionic_delta(key: str) -> float:
        a = np.asarray(terms[key], dtype=float)
        b = np.asarray(terms_inf[key], dtype=float)
        return float(0.5 * ((a[0] - b[0]) + (a[1] - b[1])))

    out = {
        "born": mean_ionic_delta("lnfugcoef_born"),
        "dh": mean_ionic_delta("lnfugcoef_ion"),
        "hc": mean_ionic_delta("lnfugcoef_hc"),
        "disp": mean_ionic_delta("lnfugcoef_disp"),
        "assoc": mean_ionic_delta("lnfugcoef_assoc"),
        "total": mean_ionic_delta("lnfugcoef_total"),
        "rho_base": rho_base,
        "rho_scaled": rho,
        "rho_inf": float(rho_inf),
    }
    return out


def _metric_summary(y_data: np.ndarray, y_model: np.ndarray) -> Tuple[float, float, float]:
    y_model = np.asarray(y_model, dtype=float)
    y_data = np.asarray(y_data, dtype=float)
    if (not np.all(np.isfinite(y_model))) or (not np.all(np.isfinite(y_data))):
        return math.nan, math.nan, math.nan
    delta = y_model - y_data
    rmse = float(np.sqrt(np.mean(delta * delta)))
    mae = float(np.mean(np.abs(delta)))
    max_abs = float(np.max(np.abs(delta)))
    return rmse, mae, max_abs


def _unique_digitized_x(series: Dict[str, Tuple[np.ndarray, np.ndarray]]) -> np.ndarray:
    x_all: List[float] = []
    for name in CONTRIBUTION_ORDER:
        x_all.extend(np.asarray(series[name][0], dtype=float).tolist())
    x_unique = np.asarray(sorted({round(float(x), 10) for x in x_all}), dtype=float)
    return x_unique


def run_scan(
    digitized_path: Path,
    output_csv: Path,
    output_plot: Path,
    output_points_csv: Path,
    rho_min: float,
    rho_max: float,
    rho_step: float,
) -> Dict[str, float]:
    series = _load_digitized_curves(digitized_path)
    params = _build_params(user_options={})
    x_union = _unique_digitized_x(series)

    scales = np.arange(float(rho_min), float(rho_max) + 0.5 * float(rho_step), float(rho_step))
    summary_rows: List[Dict[str, float]] = []
    point_rows: List[Dict[str, float]] = []

    for rho_scale in scales:
        point_cache = {float(x): _calc_point_contributions(float(x), params, float(rho_scale)) for x in x_union}
        summary: Dict[str, float] = {"rho_scale": float(rho_scale)}

        for name in CONTRIBUTION_ORDER:
            x_data, y_data = series[name]
            y_model = np.asarray([point_cache[float(round(float(x), 10))][name] for x in x_data], dtype=float)
            rmse, mae, max_abs = _metric_summary(np.asarray(y_data, dtype=float), y_model)
            summary[f"{name}_rmse"] = rmse
            summary[f"{name}_mae"] = mae
            summary[f"{name}_max_abs"] = max_abs

        summary_rows.append(summary)

        if abs(float(rho_scale) - 1.0) < 1e-12:
            for x in x_union:
                row = {"rho_scale": float(rho_scale), "x_salt": float(x)}
                cache = point_cache[float(x)]
                row["molality"] = float(_molality_for_salt_mole_fraction(float(x)))
                row["rho_base"] = cache["rho_base"]
                row["rho_scaled"] = cache["rho_scaled"]
                row["rho_inf"] = cache["rho_inf"]
                for name in CONTRIBUTION_ORDER:
                    row[name] = cache[name]
                point_rows.append(row)

    if not summary_rows:
        raise ValueError("No density-scan rows were generated.")

    finite_total_rows = [row for row in summary_rows if math.isfinite(row["total_rmse"])]
    best_total = min(finite_total_rows, key=lambda row: row["total_rmse"]) if finite_total_rows else None

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = ["rho_scale"]
        for name in CONTRIBUTION_ORDER:
            fieldnames.extend([f"{name}_rmse", f"{name}_mae", f"{name}_max_abs"])
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)

    with output_points_csv.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = ["rho_scale", "x_salt", "molality", "rho_base", "rho_scaled", "rho_inf"] + CONTRIBUTION_ORDER
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(point_rows)

    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    for name in CONTRIBUTION_ORDER:
        ax.plot(
            [row["rho_scale"] for row in summary_rows],
            [row[f"{name}_rmse"] for row in summary_rows],
            linewidth=2.0 if name in {"dh", "total"} else 1.6,
            linestyle="-" if name in {"born", "dh", "total"} else "--",
            label=TITLE_LABELS[name],
        )

    ax.axvline(1.0, color="black", linewidth=1.0, linestyle=":")
    if best_total is not None:
        ax.axvline(best_total["rho_scale"], color="green", linewidth=1.0, linestyle=":")
    ax.set_xlabel(r"Imposed density scale, $\rho / \rho_{\mathrm{EOS}}$")
    ax.set_ylabel("RMSE vs digitized Figure 6b series")
    ax.set_title("Figure 6b density sensitivity")
    ax.grid(True, alpha=0.3, color="0.7")
    ax.legend(fontsize=8)
    fig.tight_layout()
    save_plot_figure(fig, output_plot, dpi=220, bbox_inches=None)
    plt.close(fig)

    print(f"Wrote density scan summary: {output_csv}")
    print(f"Wrote baseline point table: {output_points_csv}")
    print(f"Wrote density scan plot: {output_plot}")
    invalid_total = sum(0 if math.isfinite(row["total_rmse"]) else 1 for row in summary_rows)
    if best_total is None:
        print("No finite total-RMSE cases were found in the scanned density range.")
        return {"rho_scale": math.nan, "total_rmse": math.nan, "invalid_total_scales": float(invalid_total)}

    print(f"Best total RMSE at rho_scale={best_total['rho_scale']:.4f}: {best_total['total_rmse']:.6f}")
    print(f"Density scales with non-finite total curve values: {invalid_total}")
    best_total = dict(best_total)
    best_total["invalid_total_scales"] = float(invalid_total)
    return best_total


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Scan Figure 6b contribution sensitivity to density scaling")
    parser.add_argument(
        "--digitized",
        type=Path,
        default=REPO_ROOT / "scripts" / "paper_validation" / "2020_Bulow_analysis" / "figure_6" / "figure_6b" / "data" / "Figure6b_curves.csv",
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=OUTPUT_DATA_DIR / "figure6b_density_scan_metrics.csv",
    )
    parser.add_argument(
        "--out-plot",
        type=Path,
        default=OUTPUT_PLOTS_DIR / "figure6b_density_scan_rmse.png",
    )
    parser.add_argument(
        "--out-points-csv",
        type=Path,
        default=OUTPUT_DATA_DIR / "figure6b_density_scan_baseline_points.csv",
    )
    parser.add_argument("--rho-min", type=float, default=0.94)
    parser.add_argument("--rho-max", type=float, default=1.06)
    parser.add_argument("--rho-step", type=float, default=0.005)
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    run_scan(
        digitized_path=Path(args.digitized),
        output_csv=Path(args.out_csv),
        output_plot=Path(args.out_plot),
        output_points_csv=Path(args.out_points_csv),
        rho_min=float(args.rho_min),
        rho_max=float(args.rho_max),
        rho_step=float(args.rho_step),
    )


if __name__ == "__main__":
    main()

