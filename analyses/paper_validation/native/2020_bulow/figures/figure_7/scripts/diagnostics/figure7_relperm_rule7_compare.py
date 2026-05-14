from __future__ import annotations

import csv
import math
import sys


from pathlib import Path
import sys as _bootstrap_sys
from pathlib import Path as _BootstrapPath

for _candidate in _BootstrapPath(__file__).resolve().parents:
    if (_candidate / "scripts" / "plot_outputs.py").is_file():
        if str(_candidate) not in _bootstrap_sys.path:
            _bootstrap_sys.path.insert(0, str(_candidate))
        break
else:
    raise ModuleNotFoundError("Could not locate repo root containing scripts/plot_outputs.py")
from scripts.plot_outputs import REPO_ROOT

import matplotlib
import numpy as np

FIGURE_DIR = Path(__file__).resolve().parents[1]
if str(FIGURE_DIR) not in sys.path:
    sys.path.insert(0, str(FIGURE_DIR))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts._env import require_epcsaft_install
from scripts.plot_outputs import paper_validation_dir

require_epcsaft_install()

import plot_figure_7 as fig7
from epcsaft.parameters import get_prop_dict
from scripts._epcsaft_oop import epcsaft_activity_coefficient, epcsaft_density
from _plot_common import configure_style, save_figure

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT_DIR = paper_validation_dir(__file__)


def _read_paper_points(panel: dict[str, object]) -> tuple[np.ndarray, np.ndarray]:
    path = Path(panel["data"])
    solvent = str(panel["solvent"])
    xs: list[float] = []
    ys: list[float] = []
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    use_bulow = any(str(row.get("source", "")).strip() == "Bulow 2020" for row in rows)
    for row in rows:
        if use_bulow and str(row.get("source", "")).strip() != "Bulow 2020":
            continue
        try:
            x_val = float(str(row.get("mole_fraction", "")).strip())
        except ValueError:
            x_val = math.nan
        try:
            y_val = float(str(row.get("miac", "")).strip())
        except ValueError:
            y_val = math.nan
        if not math.isfinite(y_val):
            try:
                m_val = float(str(row.get("molality", "")).strip())
                y_m = float(str(row.get("miac_m", row.get("gamma", ""))).strip())
                y_val = y_m * (1.0 + fig7.SOLVENT_MW[solvent] * m_val * 2.0)
            except ValueError:
                y_val = math.nan
        if not math.isfinite(x_val):
            try:
                m_val = float(str(row.get("molality", "")).strip())
                x_val = float(fig7._molality_for_salt_mole_fraction(np.asarray([m_val], dtype=float), solvent)[0])
            except ValueError:
                x_val = math.nan
        if math.isfinite(x_val) and math.isfinite(y_val):
            xs.append(float(x_val))
            ys.append(float(y_val))
    x_arr = np.asarray(xs, dtype=float)
    y_arr = np.asarray(ys, dtype=float)
    order = np.argsort(x_arr)
    return x_arr[order], y_arr[order]


def _curve(dataset: str, salt: str, solvent: str, x_vals: np.ndarray, user_options: dict | None = None) -> np.ndarray:
    x_ref = fig7._molality_to_species_molefraction(1e-8, salt, solvent)
    params = get_prop_dict(
        dataset, fig7._species_for_combo(salt, solvent), x_ref, fig7.T_REF, user_options=user_options or {}
    )
    species = fig7._species_for_combo(salt, solvent)
    pair_key = None
    out = np.empty_like(x_vals, dtype=float)
    m_grid = fig7._molality_for_salt_mole_fraction(x_vals, solvent)
    for idx, m in enumerate(m_grid):
        x = fig7._molality_to_species_molefraction(float(max(m, 1.0e-12)), salt, solvent)
        rho = epcsaft_density(fig7.T_REF, fig7.P_REF, x, params, phase="liq")
        vals = epcsaft_activity_coefficient(
            fig7.T_REF, rho, x, params, species=species, mean_ionic_form=True, basis="mole"
        )
        if pair_key is None:
            pair_key = fig7._resolve_pair_key(vals, salt)
        out[idx] = float(vals[pair_key])
    return out


def _metrics(y_model: np.ndarray, y_ref: np.ndarray) -> dict[str, float]:
    delta = np.asarray(y_model, dtype=float) - np.asarray(y_ref, dtype=float)
    return {
        "rmse": float(np.sqrt(np.mean(delta * delta))),
        "mae": float(np.mean(np.abs(delta))),
        "bias": float(np.mean(delta)),
        "max_abs": float(np.max(np.abs(delta))),
    }


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    configure_style()
    rows: list[dict[str, object]] = []
    for panel in fig7.PANELS:
        x_data, y_data = _read_paper_points(panel)
        visible = y_data <= float(panel["ylim"][1]) + 1e-12
        x_eval = x_data[visible]
        y_eval = y_data[visible]
        x_grid = np.linspace(float(panel["xlim"][0]), float(panel["xlim"][1]), 801)
        y_rev = _curve("2014_Held", str(panel["salt"]), str(panel["solvent"]), x_grid)
        y_adv = _curve("2020_Bulow", str(panel["salt"]), str(panel["solvent"]), x_grid)
        y_rule7 = _curve(
            "2020_Bulow",
            str(panel["salt"]),
            str(panel["solvent"]),
            x_grid,
            user_options={"elec_model": {"rel_perm": {"rule": 7}}},
        )
        y_adv_eval = np.interp(x_eval, x_grid, y_adv)
        y_rule7_eval = np.interp(x_eval, x_grid, y_rule7)
        stats_adv = _metrics(y_adv_eval, y_eval)
        stats_rule7 = _metrics(y_rule7_eval, y_eval)
        rows.append(
            {
                "panel": panel["id"],
                "salt": panel["salt"],
                "solvent": panel["solvent"],
                "baseline_rmse": stats_adv["rmse"],
                "rule7_rmse": stats_rule7["rmse"],
                "baseline_bias": stats_adv["bias"],
                "rule7_bias": stats_rule7["bias"],
                "rmse_delta_rule7_minus_baseline": float(stats_rule7["rmse"] - stats_adv["rmse"]),
                "note": "Metrics use visible paper points only; rule7 is the salt-basis dielectric variant.",
            }
        )

        fig, ax = plt.subplots(figsize=(5.4, 4.1))
        ax.scatter(
            x_eval, y_eval, s=26, facecolor="white", edgecolor="black", linewidth=0.8, label="Paper points", zorder=5
        )
        ax.plot(x_grid, y_rev, color="black", linewidth=1.8, label="ePC-SAFT revised", zorder=3)
        ax.plot(x_grid, y_adv, color="#228b22", linewidth=1.9, label="advanced rule 1", zorder=4)
        ax.plot(x_grid, y_rule7, color="#b22222", linewidth=1.9, linestyle="--", label="advanced rule 7", zorder=4)
        ax.set_xlim(*panel["xlim"])
        ax.set_ylim(*panel["ylim"])
        ax.set_xlabel(r"salt mole fraction, $x_{salt}$")
        ax.set_ylabel(r"$\gamma_{\pm}^{*}$ / -")
        ax.set_title(f"Figure {panel['id']} dielectric-rule comparison")
        ax.grid(True, alpha=0.25)
        ax.legend(loc="best", fontsize=8)
        save_figure(fig, OUT_DIR / f"figure_{panel['id']}_relperm_rule7.png")
        plt.close(fig)

    _write_csv(OUT_DIR / "figure7_relperm_rule7_summary.csv", rows)


if __name__ == "__main__":
    main()

