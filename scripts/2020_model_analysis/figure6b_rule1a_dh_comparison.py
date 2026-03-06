"""Compare ethanol-LiBr Figure 6b DH curves for dielectric rules 1, 1a, and 4.

Rule 1: species-mole-fraction linear mixing.
Rule 1a / rule 7: neutral-salt-mole-fraction linear mixing.
Rule 4: empirical solvent-fraction rule with numerical dielectric derivatives.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import matplotlib
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from figure6b_libr_ethanol_contributions import (
    _build_params,
    _calc_ln_miac_contributions,
    _load_exp_data,
    _salt_mole_fraction_from_molality,
)

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _load_digitized_dh(path: Path) -> tuple[np.ndarray, np.ndarray]:
    x_vals: list[float] = []
    y_vals: list[float] = []
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"No header found in {path}.")
        x_key = next((k for k in reader.fieldnames if k and k.strip().lower().startswith("x")), None)
        y_key = next((k for k in reader.fieldnames if k and k.strip().lower().startswith("ln")), None)
        if x_key is None or y_key is None:
            raise ValueError(f"Expected x and ln-gamma columns in {path}.")
        for row in reader:
            try:
                x = float(row[x_key])
                y = float(row[y_key])
            except (TypeError, ValueError, KeyError):
                continue
            if math.isfinite(x) and math.isfinite(y):
                x_vals.append(x)
                y_vals.append(y)
    if not x_vals:
        raise ValueError(f"No usable digitized DH points found in {path}.")
    x_arr = np.asarray(x_vals, dtype=float)
    y_arr = np.asarray(y_vals, dtype=float)
    order = np.argsort(x_arr)
    return x_arr[order], y_arr[order]


def _curve_for_user_options(user_options: dict, m_grid: np.ndarray) -> np.ndarray:
    params = _build_params(user_options=user_options)
    curves = _calc_ln_miac_contributions(m_grid, params)
    return np.asarray(curves["dh"], dtype=float)


def _rmse_mae(x_ref: np.ndarray, y_ref: np.ndarray, x_curve: np.ndarray, y_curve: np.ndarray) -> tuple[float, float]:
    y_interp = np.interp(x_ref, x_curve, y_curve)
    delta = y_interp - y_ref
    rmse = float(np.sqrt(np.mean(delta * delta)))
    mae = float(np.mean(np.abs(delta)))
    return rmse, mae


def run(data_path: Path, digitized_path: Path, output_path: Path, grid_points: int = 600) -> dict[str, float]:
    m_exp, _, _ = _load_exp_data(data_path)
    m_upper = float(np.max(m_exp))
    m_grid = np.linspace(0.0, m_upper, int(grid_points))
    x_grid = _salt_mole_fraction_from_molality(m_grid)

    variants = {
        "rule1": {
            "label": "Rule 1: species-fraction linear",
            "color": "#1f77b4",
            "linestyle": "-",
            "user_options": {"elec_model": {"rel_perm": {"rule": 1}}},
        },
        "rule1a": {
            "label": "Rule 1a: salt-fraction linear",
            "color": "#d62728",
            "linestyle": "--",
            "user_options": {"elec_model": {"rel_perm": {"rule": 7}}},
        },
        "rule4_num": {
            "label": "Rule 4: empirical + numerical d eps/dx",
            "color": "#2ca02c",
            "linestyle": "-.",
            "user_options": {"elec_model": {"rel_perm": {"rule": 4, "differential_mode": "numerical"}}},
        },
    }

    x_dh, y_dh = _load_digitized_dh(digitized_path)
    metrics: dict[str, float] = {"x_max": float(np.max(x_grid))}
    curves: dict[str, np.ndarray] = {}

    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    ax.scatter(x_dh, y_dh, facecolors="none", edgecolors="black", s=42, linewidths=1.2, label="Digitized DH (paper)")

    for key, cfg in variants.items():
        curve = _curve_for_user_options(cfg["user_options"], m_grid)
        curves[key] = curve
        rmse, mae = _rmse_mae(x_dh, y_dh, x_grid, curve)
        metrics[f"rmse_{key}"] = rmse
        metrics[f"mae_{key}"] = mae
        ax.plot(x_grid, curve, color=cfg["color"], linewidth=2.0, linestyle=cfg["linestyle"], label=cfg["label"])

    ax.set_xlim(0.0, max(0.2, float(np.max(x_grid))))
    ax.set_ylim(-3.0, 0.5)
    ax.set_xlabel(r"salt mole fraction, $x_{salt}$")
    ax.set_ylabel(r"DH contribution to $\ln(\gamma_{\pm}^{*})$")
    ax.legend(frameon=True, facecolor="white", edgecolor="black")
    ax.grid(False)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)

    metrics["delta_rmse_rule1a_vs_rule1"] = metrics["rmse_rule1a"] - metrics["rmse_rule1"]
    metrics["delta_rmse_rule4_num_vs_rule1"] = metrics["rmse_rule4_num"] - metrics["rmse_rule1"]
    metrics["delta_mae_rule1a_vs_rule1"] = metrics["mae_rule1a"] - metrics["mae_rule1"]
    metrics["delta_mae_rule4_num_vs_rule1"] = metrics["mae_rule4_num"] - metrics["mae_rule1"]
    return metrics


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=REPO_ROOT / "data" / "MIAC" / "ethanol" / "ethanol-LiBr.csv")
    parser.add_argument("--digitized", type=Path, default=Path(r"C:\Users\Tanner\Downloads\DH_from_article_2020.csv"))
    parser.add_argument("--output", type=Path, default=REPO_ROOT / "scripts" / "2020_model_analysis" / "output" / "figure6b_libr_ethanol_rule_compare_dh.png")
    parser.add_argument("--grid-points", type=int, default=600)
    args = parser.parse_args()

    metrics = run(args.data, args.digitized, args.output, grid_points=args.grid_points)
    for key in sorted(metrics):
        print(f"{key}={metrics[key]:.10f}")
    print(f"output={args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
