from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from data.epcsaft_properties import get_prop_dict
from pcsaft import pcsaft_dielc_eval

matplotlib.use("Agg")
import matplotlib.pyplot as plt

T_REF = 298.15
SPECIES = ["Li+", "Br-", "Ethanol"]
REPRESENTATIVE_X = np.asarray([0.0, 0.05, 0.10, 0.16732, 0.20], dtype=float)
ALPHA_RULE4 = 7.01

RULE_CONFIGS = [
    {
        "rule": 1,
        "label": "Rule 1",
        "color": "black",
        "linestyle": "-",
        "user_options": {"elec_model": {"rel_perm": {"rule": 1, "differential_mode": "analytical"}}},
    },
    {
        "rule": 7,
        "label": "Rule 1a",
        "color": "tab:green",
        "linestyle": "--",
        "user_options": {"elec_model": {"rel_perm": {"rule": 7, "differential_mode": "analytical"}}},
    },
    {
        "rule": 4,
        "label": "Rule 4",
        "color": "tab:orange",
        "linestyle": "-.",
        "user_options": {"elec_model": {"rel_perm": {"rule": 4, "differential_mode": "numerical"}}},
    },
]


def _salt_molefraction_to_species(x_salt: np.ndarray) -> np.ndarray:
    xs = np.asarray(x_salt, dtype=float)
    if np.any(xs < 0.0) or np.any(xs >= 1.0):
        raise ValueError("Salt mole fraction grid must satisfy 0 <= x_salt < 1.")
    x_ion_each = xs / (1.0 + xs)
    x_solvent = (1.0 - xs) / (1.0 + xs)
    return np.column_stack((x_ion_each, x_ion_each, x_solvent))


def _dxspecies_dxsalt(x_salt: np.ndarray) -> np.ndarray:
    xs = np.asarray(x_salt, dtype=float)
    denom = (1.0 + xs) ** 2
    return np.column_stack((1.0 / denom, 1.0 / denom, -2.0 / denom))


def _build_params(user_options: dict) -> dict:
    x_ref = np.asarray([1.0e-12, 1.0e-12, 1.0 - 2.0e-12], dtype=float)
    return get_prop_dict("bulow_2020", SPECIES, x_ref, T_REF, user_options=user_options)


def _eval_runtime(x_species: np.ndarray, params: dict) -> tuple[np.ndarray, np.ndarray]:
    eps = np.empty(x_species.shape[0], dtype=float)
    deps = np.empty_like(x_species, dtype=float)
    for idx, x in enumerate(x_species):
        eps_i, deps_i = pcsaft_dielc_eval(x, params)
        eps[idx] = float(eps_i)
        deps[idx, :] = np.asarray(deps_i, dtype=float)
    return eps, deps


def _formula_rule1(x_species: np.ndarray, dielc: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    eps = x_species @ dielc
    deps = np.tile(dielc, (x_species.shape[0], 1))
    return eps, deps


def _formula_rule1a(x_species: np.ndarray, dielc: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x_ion = x_species[:, 0] + x_species[:, 1]
    denom = 2.0 - x_ion
    eps_salt = 0.5 * (dielc[0] + dielc[1])
    solvent_sum = x_species[:, 2] * dielc[2]
    eps_num = 2.0 * solvent_sum + x_ion * eps_salt
    eps = eps_num / denom
    deps = np.empty_like(x_species)
    common_ion = (eps_salt * denom + eps_num) / (denom * denom)
    deps[:, 0] = common_ion
    deps[:, 1] = common_ion
    deps[:, 2] = 2.0 * dielc[2] / denom
    return eps, deps


def _formula_rule4(x_species: np.ndarray, dielc: np.ndarray, mw: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x_ion = x_species[:, 0] + x_species[:, 1]
    solvent_mass = x_species[:, 2] * mw[2]
    eps_sf = (x_species[:, 2] * mw[2] * dielc[2]) / solvent_mass
    den = 1.0 + ALPHA_RULE4 * x_ion
    eps = eps_sf / den
    deps = np.empty_like(x_species)
    deps[:, 0] = -ALPHA_RULE4 * eps_sf / (den * den)
    deps[:, 1] = -ALPHA_RULE4 * eps_sf / (den * den)
    deps[:, 2] = (1.0 / den) * (mw[2] / solvent_mass) * (dielc[2] - eps_sf)
    return eps, deps


def _formula_eval(rule: int, x_species: np.ndarray, dielc: np.ndarray, mw: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if rule == 1:
        return _formula_rule1(x_species, dielc)
    if rule == 7:
        return _formula_rule1a(x_species, dielc)
    if rule == 4:
        return _formula_rule4(x_species, dielc, mw)
    raise ValueError(f"Unsupported rule for closed-form evaluation: {rule}")


def _depsalt_from_deps(x_salt: np.ndarray, deps_dx: np.ndarray) -> np.ndarray:
    return np.einsum("ij,ij->i", deps_dx, _dxspecies_dxsalt(x_salt))


def _print_summary(name: str, x_points: np.ndarray, eps: np.ndarray, deps: np.ndarray, deps_salt: np.ndarray, formula_eps: np.ndarray, formula_deps: np.ndarray, formula_deps_salt: np.ndarray) -> None:
    print(f"\n{name}")
    print("  x_salt      eps_r     d eps/dx_Li   d eps/dx_Br   d eps/dx_EtOH   d eps/dx_LiBr")
    for x_val, eps_val, dep_row, dep_salt in zip(x_points, eps, deps, deps_salt):
        print(
            f"  {x_val:7.5f}  {eps_val:9.6f}  {dep_row[0]:12.6f}  {dep_row[1]:12.6f}  {dep_row[2]:14.6f}  {dep_salt:14.6f}"
        )
    print(f"  max |eps_runtime - eps_formula|      = {np.max(np.abs(eps - formula_eps)):.6e}")
    print(f"  max |deps_runtime - deps_formula|    = {np.max(np.abs(deps - formula_deps)):.6e}")
    print(f"  max |d eps/dx_s - formula|           = {np.max(np.abs(deps_salt - formula_deps_salt)):.6e}")



def _save_primary_plot(output_path: Path, x_salt_grid: np.ndarray, curves: list[dict]) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(8.4, 10.2), sharex=True)
    for curve in curves:
        axes[0].plot(x_salt_grid, curve["eps"], color=curve["color"], linestyle=curve["linestyle"], lw=1.8, label=curve["label"])
        axes[1].plot(x_salt_grid, curve["deps"][:, 0], color=curve["color"], linestyle=curve["linestyle"], lw=1.8, label=curve["label"])
        axes[2].plot(x_salt_grid, curve["deps"][:, 2], color=curve["color"], linestyle=curve["linestyle"], lw=1.8, label=curve["label"])

    axes[0].set_ylabel(r"$\varepsilon_r$")
    axes[1].set_ylabel(r"$\partial \varepsilon_r / \partial x_{Li^+}$")
    axes[2].set_ylabel(r"$\partial \varepsilon_r / \partial x_{EtOH}$")
    axes[2].set_xlabel(r"$x_{LiBr}$")
    axes[0].set_title("Ethanol-LiBr dielectric diagnostics from the C++ dielectric engine")
    for ax in axes:
        ax.grid(True, alpha=0.25)
        ax.legend(frameon=False, loc="best")
    axes[0].set_xlim(float(x_salt_grid[0]), float(x_salt_grid[-1]))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)



def _save_salt_slope_plot(output_path: Path, x_salt_grid: np.ndarray, curves: list[dict]) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    for curve in curves:
        ax.plot(x_salt_grid, curve["deps_salt"], color=curve["color"], linestyle=curve["linestyle"], lw=1.9, label=curve["label"])
    ax.set_title(r"Ethanol-LiBr dielectric slope on the paper axis: $d \varepsilon_r / d x_{LiBr}$")
    ax.set_xlabel(r"$x_{LiBr}$")
    ax.set_ylabel(r"$d \varepsilon_r / d x_{LiBr}$")
    ax.set_xlim(float(x_salt_grid[0]), float(x_salt_grid[-1]))
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)



def run(output_path: Path, salt_slope_output_path: Path, xmin: float, xmax: float, points: int) -> Path:
    x_salt_grid = np.linspace(float(xmin), float(xmax), int(points))
    x_species_grid = _salt_molefraction_to_species(x_salt_grid)
    x_species_rep = _salt_molefraction_to_species(REPRESENTATIVE_X)

    metrics: dict[str, dict[str, object]] = {}
    curves: list[dict] = []

    for cfg in RULE_CONFIGS:
        params = _build_params(cfg["user_options"])
        dielc = np.asarray(params["dielc"], dtype=float)
        mw = np.asarray(params["MW"], dtype=float)

        eps_runtime, deps_runtime = _eval_runtime(x_species_grid, params)
        eps_formula, deps_formula = _formula_eval(int(cfg["rule"]), x_species_grid, dielc, mw)
        deps_salt_runtime = _depsalt_from_deps(x_salt_grid, deps_runtime)
        deps_salt_formula = _depsalt_from_deps(x_salt_grid, deps_formula)

        eps_rep, deps_rep = _eval_runtime(x_species_rep, params)
        eps_rep_formula, deps_rep_formula = _formula_eval(int(cfg["rule"]), x_species_rep, dielc, mw)
        deps_salt_rep = _depsalt_from_deps(REPRESENTATIVE_X, deps_rep)
        deps_salt_rep_formula = _depsalt_from_deps(REPRESENTATIVE_X, deps_rep_formula)

        metrics[str(cfg["label"])] = {
            "dielc": dielc.tolist(),
            "max_abs_eps_diff": float(np.max(np.abs(eps_runtime - eps_formula))),
            "max_abs_deps_diff": float(np.max(np.abs(deps_runtime - deps_formula))),
            "max_abs_deps_salt_diff": float(np.max(np.abs(deps_salt_runtime - deps_salt_formula))),
            "representative_x_salt": REPRESENTATIVE_X.tolist(),
            "representative_eps": eps_rep.tolist(),
            "representative_deps_dx": deps_rep.tolist(),
            "representative_deps_dxsalt": deps_salt_rep.tolist(),
        }
        _print_summary(
            str(cfg["label"]),
            REPRESENTATIVE_X,
            eps_rep,
            deps_rep,
            deps_salt_rep,
            eps_rep_formula,
            deps_rep_formula,
            deps_salt_rep_formula,
        )

        curves.append(
            {
                "label": str(cfg["label"]),
                "color": cfg["color"],
                "linestyle": cfg["linestyle"],
                "eps": eps_runtime,
                "deps": deps_runtime,
                "deps_salt": deps_salt_runtime,
            }
        )

    _save_primary_plot(output_path, x_salt_grid, curves)
    _save_salt_slope_plot(salt_slope_output_path, x_salt_grid, curves)

    metrics_path = output_path.with_suffix(".json")
    metrics_path.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    print(f"\nSaved plot to {output_path}")
    print(f"Saved plot to {salt_slope_output_path}")
    print(f"Saved metrics to {metrics_path}")
    return output_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot dielectric constant and composition derivatives for ethanol-LiBr.")
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "scripts" / "2020_model_analysis" / "output" / "figure6b_libr_ethanol_dielc_diagnostics.png",
    )
    parser.add_argument(
        "--salt-slope-output",
        type=Path,
        default=REPO_ROOT / "scripts" / "2020_model_analysis" / "output" / "figure6b_libr_ethanol_dielc_salt_slope.png",
    )
    parser.add_argument("--xmin", type=float, default=0.0)
    parser.add_argument("--xmax", type=float, default=0.20)
    parser.add_argument("--points", type=int, default=500)
    args = parser.parse_args()
    run(args.output, args.salt_slope_output, args.xmin, args.xmax, args.points)
