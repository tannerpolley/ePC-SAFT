from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import matplotlib
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
FIGURE6_DIR = SCRIPT_DIR.parent
ANALYSIS_ROOT = FIGURE6_DIR.parent
REPO_ROOT = ANALYSIS_ROOT.parents[2]
FIGURE6B_DIAG_DIR = FIGURE6_DIR / "figure_6b" / "diagnostics"

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(FIGURE6B_DIAG_DIR) not in sys.path:
    sys.path.insert(0, str(FIGURE6B_DIAG_DIR))

from data.epcsaft_properties import get_prop_dict
from figure6b_libr_ethanol_contributions import (
    _calc_ln_miac_contributions,
    _load_exp_data,
    _molality_to_species_molefraction,
    _salt_mole_fraction_from_molality,
    P_REF,
    T_REF,
)

try:
    import feos
    import si_units as si
except ImportError as exc:  # pragma: no cover
    raise RuntimeError("feos and si_units must be importable in the PC-SAFT environment.") from exc

matplotlib.use("Agg")
import matplotlib.pyplot as plt

R_GAS = 8.31446261815324
RT = R_GAS * T_REF
DATA_PATH = REPO_ROOT / "data" / "MIAC" / "ethanol" / "ethanol-LiBr.csv"
PURE_JSON = SCRIPT_DIR / "feos_libr_ethanol_pure.json"
BINARY_JSON = SCRIPT_DIR / "feos_libr_ethanol_binary.json"
CURVES_CSV = SCRIPT_DIR / "figure6_feos_curves.csv"
STATS_CSV = SCRIPT_DIR / "figure6_feos_vs_pcsaft_stats.csv"
FIG6A_PNG = SCRIPT_DIR / "figure_6a_feos_vs_pcsaft.png"
FIG6B_PNG = SCRIPT_DIR / "figure_6b_feos_vs_pcsaft_mu.png"
BOOKKEEPING_PNG = SCRIPT_DIR / "figure_6b_feos_bookkeeping.png"
NOTES_MD = SCRIPT_DIR / "feos_analysis_notes.md"

REPO_SPECIES = ["Li+", "Br-", "Ethanol"]
FEOS_COMPONENTS = ["ethanol", "lithium ion", "bromide ion"]
FEOS_CATION = 1
FEOS_ANION = 2
MW_ETHANOL = 46.068e-3
MW_LI = 6.941e-3
MW_BR = 79.904e-3
CONTRIBUTIONS = ("hc", "disp", "assoc", "dh", "born")
COMPARE_KEYS = ("total",) + CONTRIBUTIONS
LABEL_TO_TERM = {
    "Hard Sphere": "hc",
    "Hard Chain": "hc",
    "Dispersion": "disp",
    "Association": "assoc",
    "Ionic": "dh",
    "Born": "born",
}
TERM_COLORS = {
    "hc": "tab:blue",
    "disp": "tab:orange",
    "assoc": "tab:green",
    "dh": "tab:red",
    "born": "tab:purple",
    "total": "black",
}
PCSAFT_MARKER = "x"
PCSAFT_MARK_EVERY = 25
FEOS_LINESTYLE = ":"


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _mean_ionic(a: np.ndarray, b: np.ndarray) -> float:
    return 0.5 * ((float(a[FEOS_CATION]) - float(b[FEOS_CATION])) + (float(a[FEOS_ANION]) - float(b[FEOS_ANION])))


def _feos_dimless_energy(quantity) -> float:
    return float(quantity / (si.JOULE / si.MOL)) / RT


def _feos_dimless_pressure(quantity) -> float:
    return float(quantity / si.PASCAL)


def _map_feos_pairs(pairs, converter) -> dict[str, float]:
    out = {term: 0.0 for term in CONTRIBUTIONS}
    for label, quantity in pairs:
        term = LABEL_TO_TERM.get(label)
        if term is None:
            continue
        out[term] += converter(quantity)
    return out


def _build_custom_feos_files() -> tuple[Path, Path]:
    x_ref = _molality_to_species_molefraction(1.0e-12)
    params = get_prop_dict("bulow_2020", REPO_SPECIES, x_ref, T_REF)
    pure_records = [
        {
            "identifier": {"name": "ethanol", "formula": "C2H6O"},
            "m": float(params["m"][2]),
            "sigma": float(params["s"][2]),
            "epsilon_k": float(params["e"][2]),
            "association_sites": [
                {
                    "na": 1.0,
                    "nb": 1.0,
                    "kappa_ab": float(params["vol_a"][2]),
                    "epsilon_k_ab": float(params["e_assoc"][2]),
                }
            ],
            "permittivity_record": {"ExperimentalData": {"data": [[float(T_REF), float(params["dielc"][2])]]}},
            "molarweight": MW_ETHANOL,
        },
        {
            "identifier": {"name": "lithium ion", "formula": "Li+"},
            "m": float(params["m"][0]),
            "sigma": float(params["s"][0]),
            "epsilon_k": float(params["e"][0]),
            "z": 1.0,
            "permittivity_record": {"ExperimentalData": {"data": [[float(T_REF), float(params["dielc"][0])]]}},
            "molarweight": MW_LI,
        },
        {
            "identifier": {"name": "bromide ion", "formula": "Br-"},
            "m": float(params["m"][1]),
            "sigma": float(params["s"][1]),
            "epsilon_k": float(params["e"][1]),
            "z": -1.0,
            "permittivity_record": {"ExperimentalData": {"data": [[float(T_REF), float(params["dielc"][1])]]}},
            "molarweight": MW_BR,
        },
    ]
    binary_records = [
        {"id1": {"name": "ethanol"}, "id2": {"name": "lithium ion"}, "k_ij": [float(params["k_ij"][2][0]), 0.0, 0.0, 0.0]},
        {"id1": {"name": "ethanol"}, "id2": {"name": "bromide ion"}, "k_ij": [float(params["k_ij"][2][1]), 0.0, 0.0, 0.0]},
        {"id1": {"name": "lithium ion"}, "id2": {"name": "bromide ion"}, "k_ij": [float(params["k_ij"][0][1]), 0.0, 0.0, 0.0]},
    ]
    _write_json(PURE_JSON, pure_records)
    _write_json(BINARY_JSON, binary_records)
    return PURE_JSON, BINARY_JSON


def _build_feos_eos() -> feos.EquationOfState:
    pure_path, binary_path = _build_custom_feos_files()
    parameters = feos.Parameters.from_json(FEOS_COMPONENTS, str(pure_path), str(binary_path))
    return feos.EquationOfState.epcsaft(parameters, epcsaft_variant="advanced")


def _state_for_molality(eos, molality: float):
    m_eval = float(molality) if molality > 0.0 else 1.0e-12
    x_repo = _molality_to_species_molefraction(m_eval)
    x_feos = np.asarray([x_repo[2], x_repo[0], x_repo[1]], dtype=float)
    state = feos.State(
        eos,
        temperature=T_REF * si.KELVIN,
        pressure=P_REF * si.PASCAL,
        molefracs=x_feos,
        total_moles=si.MOL,
    )
    return state


def _feos_mean_mu_terms(state) -> dict[str, float]:
    mu_c = _map_feos_pairs(state.chemical_potential_contributions(FEOS_CATION, feos.Contributions.Residual), _feos_dimless_energy)
    mu_a = _map_feos_pairs(state.chemical_potential_contributions(FEOS_ANION, feos.Contributions.Residual), _feos_dimless_energy)
    return {term: 0.5 * (mu_c[term] + mu_a[term]) for term in CONTRIBUTIONS}


def _feos_mean_lnfug_terms(state) -> tuple[dict[str, float], float]:
    mu_terms = _feos_mean_mu_terms(state)
    p_terms = _map_feos_pairs(state.pressure_contributions(), _feos_dimless_pressure)
    pressure_pa = float(state.pressure() / si.PASCAL)
    volume_m3 = float(state.volume / (si.METER**3))
    z_total = pressure_pa * volume_m3 / RT
    z_residual = z_total - 1.0
    corrections: dict[str, float] = {}
    for term in CONTRIBUTIONS:
        if abs(z_residual) <= 1.0e-14 or z_total <= 0.0:
            corrections[term] = 0.0
        else:
            z_alpha = p_terms[term] * volume_m3 / RT
            corrections[term] = -(z_alpha / z_residual) * math.log(z_total)
    return ({term: mu_terms[term] + corrections[term] for term in CONTRIBUTIONS}, z_total)


def _compute_feos_curves(m_grid: np.ndarray) -> dict[str, np.ndarray]:
    eos = _build_feos_eos()
    ref = _state_for_molality(eos, 0.0)
    ref_lnphi = np.asarray(ref.ln_phi(), dtype=float)
    ref_mu = _feos_mean_mu_terms(ref)
    ref_lnfug, _ = _feos_mean_lnfug_terms(ref)

    out = {
        "total": np.empty_like(m_grid, dtype=float),
        "mu_sum": np.empty_like(m_grid, dtype=float),
        "lnfug_sum": np.empty_like(m_grid, dtype=float),
        "z_total": np.empty_like(m_grid, dtype=float),
        "closure_total_minus_mu_sum": np.empty_like(m_grid, dtype=float),
        "closure_total_minus_lnfug_sum": np.empty_like(m_grid, dtype=float),
    }
    for term in CONTRIBUTIONS:
        out[f"{term}_mu"] = np.empty_like(m_grid, dtype=float)
        out[f"{term}_lnfug"] = np.empty_like(m_grid, dtype=float)

    for idx, molality in enumerate(m_grid):
        state = _state_for_molality(eos, float(molality))
        lnphi = np.asarray(state.ln_phi(), dtype=float)
        out["total"][idx] = _mean_ionic(lnphi, ref_lnphi)
        mu_terms = _feos_mean_mu_terms(state)
        lnfug_terms, z_total = _feos_mean_lnfug_terms(state)
        mu_sum = 0.0
        lnfug_sum = 0.0
        for term in CONTRIBUTIONS:
            mu_delta = mu_terms[term] - ref_mu[term]
            lnfug_delta = lnfug_terms[term] - ref_lnfug[term]
            out[f"{term}_mu"][idx] = mu_delta
            out[f"{term}_lnfug"][idx] = lnfug_delta
            mu_sum += mu_delta
            lnfug_sum += lnfug_delta
        out["mu_sum"][idx] = mu_sum
        out["lnfug_sum"][idx] = lnfug_sum
        out["z_total"][idx] = z_total
        out["closure_total_minus_mu_sum"][idx] = out["total"][idx] - mu_sum
        out["closure_total_minus_lnfug_sum"][idx] = out["total"][idx] - lnfug_sum
    return out


def _rmse(x_exp: np.ndarray, y_exp: np.ndarray, x_grid: np.ndarray, y_grid: np.ndarray) -> float:
    model = np.interp(x_exp, x_grid, y_grid)
    delta = model - y_exp
    return float(np.sqrt(np.mean(delta * delta)))


def _marker_plot(ax, x: np.ndarray, y: np.ndarray, color: str, label: str, zorder: int) -> None:
    ax.plot(
        x,
        y,
        linestyle="None",
        marker=PCSAFT_MARKER,
        markevery=PCSAFT_MARK_EVERY,
        markersize=4.2,
        markeredgewidth=1.0,
        color=color,
        label=label,
        zorder=zorder,
    )


def _build_stats_rows(pcsaft_curves: dict[str, np.ndarray], feos_curves: dict[str, np.ndarray]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for key in COMPARE_KEYS:
        pc = np.asarray(pcsaft_curves["total"] if key == "total" else pcsaft_curves[key], dtype=float)
        fe = np.asarray(feos_curves["total"] if key == "total" else feos_curves[f"{key}_mu"], dtype=float)
        delta = fe - pc
        denom = np.maximum(np.abs(pc), 1.0e-12)
        rows.append(
            {
                "series": key,
                "pcsaft_min": float(np.min(pc)),
                "pcsaft_max": float(np.max(pc)),
                "feos_min": float(np.min(fe)),
                "feos_max": float(np.max(fe)),
                "mean_signed_delta_feos_minus_pcsaft": float(np.mean(delta)),
                "mean_abs_delta": float(np.mean(np.abs(delta))),
                "max_abs_delta": float(np.max(np.abs(delta))),
                "rmse": float(np.sqrt(np.mean(delta * delta))),
                "mean_abs_percent_delta_vs_pcsaft": float(np.mean(np.abs(delta) / denom) * 100.0),
            }
        )
    rows.append(
        {
            "series": "feos_bookkeeping",
            "pcsaft_min": None,
            "pcsaft_max": None,
            "feos_min": None,
            "feos_max": None,
            "mean_signed_delta_feos_minus_pcsaft": None,
            "mean_abs_delta": float(np.mean(np.abs(feos_curves["closure_total_minus_mu_sum"]))),
            "max_abs_delta": float(np.max(np.abs(feos_curves["closure_total_minus_mu_sum"]))),
            "rmse": float(np.sqrt(np.mean(np.square(feos_curves["closure_total_minus_mu_sum"])))),
            "mean_abs_percent_delta_vs_pcsaft": None,
        }
    )
    return rows


def _write_notes(
    x_exp: np.ndarray,
    y_exp: np.ndarray,
    x_grid: np.ndarray,
    pcsaft_curves: dict[str, np.ndarray],
    feos_curves: dict[str, np.ndarray],
    stats_rows: list[dict[str, object]],
) -> None:
    stats_lookup = {row["series"]: row for row in stats_rows}
    notes = [
        "# Figure 6 feos Analysis Notes",
        "",
        "## Source provenance",
        "",
        "- This Figure 6 `feos` run is not independent of this repo's parameter set.",
        "- The custom `feos` ePC-SAFT pure/binary JSON files are built directly from `get_prop_dict(\"bulow_2020\", [\"Li+\", \"Br-\", \"Ethanol\"], ...)` in this repo.",
        "- That means the parameter source for the current Figure 6 `feos` comparison is the `PC-SAFT` repo's `bulow_2020` values, not stock `feos` parameter files.",
        "- The comparison is still useful for checking implementation behavior and bookkeeping across molality, but it is not an external parameter-set validation in the same sense as the water-only Figure 3 work.",
        "",
        "## Outputs",
        "",
        "- `figure_6a_feos_vs_pcsaft.png` compares Bulow-2020 experimental points, the current repo total curve, and the `feos` total curve.",
        "- `figure_6b_feos_vs_pcsaft_mu.png` compares the current repo Figure 6b-style $\\mu$-basis contributions against `feos` $\\mu$-basis contributions.",
        "- `figure_6b_feos_bookkeeping.png` shows how the `feos` total compares with the summed `feos` $\\mu$ and reconstructed per-term $\\ln\\varphi$ contributions across molality.",
        "- `figure6_feos_vs_pcsaft_stats.csv` reports per-series deltas between the `PC-SAFT` and `feos` curves.",
        "",
        "## Important caveat",
        "",
        "- `feos` Born uses the hard-sphere diameter path in its ePC-SAFT implementation rather than the repo's explicit `d_born` values, so the Born comparison is not a strict apples-to-apples port even though the other pure/binary values are ported from `bulow_2020`.",
        "",
        "## Fit summary",
        "",
        f"- Repo total RMSE vs Bulow-2020 data: {_rmse(x_exp, y_exp, x_grid, pcsaft_curves['total']):.4f}",
        f"- feos total RMSE vs Bulow-2020 data: {_rmse(x_exp, y_exp, x_grid, feos_curves['total']):.4f}",
        f"- Total mean |feos - pcsaft|: {float(stats_lookup['total']['mean_abs_delta']):.6f}",
        f"- Total max |feos - pcsaft|: {float(stats_lookup['total']['max_abs_delta']):.6f}",
        f"- Max feos closure gap |total - mu_sum|: {float(np.max(np.abs(feos_curves['closure_total_minus_mu_sum']))):.4f}",
        f"- Max feos closure gap |total - lnfug_sum|: {float(np.max(np.abs(feos_curves['closure_total_minus_lnfug_sum']))):.4f}",
    ]
    NOTES_MD.write_text("\n".join(notes) + "\n", encoding="utf-8")


def _plot_figure6a(x_exp: np.ndarray, y_exp: np.ndarray, x_grid: np.ndarray, pcsaft_curves: dict[str, np.ndarray], feos_curves: dict[str, np.ndarray]) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    ax.scatter(x_exp, y_exp, color="black", s=34, label="Experimental data (Bulow 2020)", zorder=6)
    _marker_plot(ax, x_grid, pcsaft_curves["total"], TERM_COLORS["total"], "PC-SAFT total", 4)
    ax.plot(x_grid, feos_curves["total"], color=TERM_COLORS["total"], linewidth=2.2, linestyle=FEOS_LINESTYLE, label="feos total", zorder=5)
    ax.set_xlim(0.0, 0.2)
    ax.set_ylim(0.0, 4.0)
    ax.set_xlabel(r"salt mole fraction, $x_{salt}$")
    ax.set_ylabel(r"$\ln(\gamma_{\pm}^{*})$")
    ax.set_title("Figure 6a comparison: LiBr in ethanol at 298.15 K and 1 bar")
    ax.grid(True, alpha=0.3, color="0.7")
    ax.legend(fontsize=9)
    ax.text(
        0.98,
        0.97,
        f"Repo RMSE={_rmse(x_exp, y_exp, x_grid, pcsaft_curves['total']):.3f}\nfeos RMSE={_rmse(x_exp, y_exp, x_grid, feos_curves['total']):.3f}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        bbox={"facecolor": "white", "edgecolor": "black", "alpha": 0.95, "boxstyle": "round,pad=0.25"},
    )
    fig.tight_layout()
    fig.savefig(FIG6A_PNG, dpi=220)
    plt.close(fig)


def _plot_figure6b(x_exp: np.ndarray, y_exp: np.ndarray, x_grid: np.ndarray, pcsaft_curves: dict[str, np.ndarray], feos_curves: dict[str, np.ndarray]) -> None:
    fig, ax = plt.subplots(figsize=(8.8, 6.2))
    ax.scatter(x_exp, y_exp, color="black", s=26, label="Experimental total (Bulow 2020)", zorder=7)
    for term in ("born", "dh", "hc", "disp", "assoc", "total"):
        _marker_plot(
            ax,
            x_grid,
            pcsaft_curves["total"] if term == "total" else pcsaft_curves[term],
            TERM_COLORS[term],
            f"PC-SAFT {term}",
            4,
        )
        ax.plot(
            x_grid,
            feos_curves["total"] if term == "total" else feos_curves[f"{term}_mu"],
            color=TERM_COLORS[term],
            linestyle=FEOS_LINESTYLE,
            linewidth=2.0 if term == "total" else 1.9,
            label=f"feos {term}",
            zorder=5,
        )
    ax.set_xlim(0.0, 0.2)
    ax.set_ylim(-3.0, 4.0)
    ax.set_xlabel(r"salt mole fraction, $x_{salt}$")
    ax.set_ylabel(r"$\ln(\gamma_{\pm}^{*})$")
    ax.set_title("Figure 6b comparison: $\\mu$-basis contributions")
    ax.grid(True, alpha=0.3, color="0.7")
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(FIG6B_PNG, dpi=220)
    plt.close(fig)


def _plot_bookkeeping(x_grid: np.ndarray, feos_curves: dict[str, np.ndarray]) -> None:
    fig, ax = plt.subplots(figsize=(8.4, 5.8))
    ax.plot(x_grid, feos_curves["total"], color="black", linewidth=2.2, label="feos total")
    ax.plot(x_grid, feos_curves["mu_sum"], color="crimson", linewidth=2.0, linestyle="--", label="feos summed $\\mu$ contributions")
    ax.plot(x_grid, feos_curves["lnfug_sum"], color="tab:blue", linewidth=2.0, linestyle="-.", label="feos summed reconstructed $\\ln\\varphi^\\alpha$")
    ax.set_xlim(0.0, 0.2)
    ax.set_xlabel(r"salt mole fraction, $x_{salt}$")
    ax.set_ylabel(r"$\ln(\gamma_{\pm}^{*})$")
    ax.set_title("feos Figure 6 bookkeeping across molality")
    ax.grid(True, alpha=0.3, color="0.7")
    ax.legend(fontsize=9)
    ax.text(
        0.98,
        0.97,
        f"max |total-mu_sum|={float(np.max(np.abs(feos_curves['closure_total_minus_mu_sum']))):.3f}\nmax |total-lnfug_sum|={float(np.max(np.abs(feos_curves['closure_total_minus_lnfug_sum']))):.3f}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        bbox={"facecolor": "white", "edgecolor": "black", "alpha": 0.95, "boxstyle": "round,pad=0.25"},
    )
    fig.tight_layout()
    fig.savefig(BOOKKEEPING_PNG, dpi=220)
    plt.close(fig)


def run_analysis() -> None:
    m_exp, x_exp, y_exp = _load_exp_data(DATA_PATH)
    m_grid = np.linspace(0.0, float(np.max(m_exp)), 1201)
    x_grid = _salt_mole_fraction_from_molality(m_grid)
    x_ref = _molality_to_species_molefraction(1.0e-12)
    pcsaft_params = get_prop_dict("bulow_2020", REPO_SPECIES, x_ref, T_REF)
    pcsaft_curves = _calc_ln_miac_contributions(m_grid, pcsaft_params, method="mu")
    feos_curves = _compute_feos_curves(m_grid)

    rows: list[dict[str, object]] = []
    for idx, molality in enumerate(m_grid):
        row: dict[str, object] = {
            "molality": float(molality),
            "x_salt": float(x_grid[idx]),
            "pcsaft_total": float(pcsaft_curves["total"][idx]),
            "feos_total": float(feos_curves["total"][idx]),
            "feos_mu_sum": float(feos_curves["mu_sum"][idx]),
            "feos_lnfug_sum": float(feos_curves["lnfug_sum"][idx]),
            "feos_z_total": float(feos_curves["z_total"][idx]),
            "feos_closure_total_minus_mu_sum": float(feos_curves["closure_total_minus_mu_sum"][idx]),
            "feos_closure_total_minus_lnfug_sum": float(feos_curves["closure_total_minus_lnfug_sum"][idx]),
        }
        for term in CONTRIBUTIONS:
            row[f"pcsaft_{term}_mu"] = float(pcsaft_curves[term][idx])
            row[f"feos_{term}_mu"] = float(feos_curves[f"{term}_mu"][idx])
            row[f"feos_{term}_lnfug"] = float(feos_curves[f"{term}_lnfug"][idx])
        rows.append(row)

    stats_rows = _build_stats_rows(pcsaft_curves, feos_curves)
    _write_csv(CURVES_CSV, rows)
    _write_csv(STATS_CSV, stats_rows)
    _plot_figure6a(x_exp, y_exp, x_grid, pcsaft_curves, feos_curves)
    _plot_figure6b(x_exp, y_exp, x_grid, pcsaft_curves, feos_curves)
    _plot_bookkeeping(x_grid, feos_curves)
    _write_notes(x_exp, y_exp, x_grid, pcsaft_curves, feos_curves, stats_rows)

    print(f"Wrote {CURVES_CSV}")
    print(f"Wrote {STATS_CSV}")
    print(f"Wrote {FIG6A_PNG}")
    print(f"Wrote {FIG6B_PNG}")
    print(f"Wrote {BOOKKEEPING_PNG}")
    print(f"Wrote {NOTES_MD}")


if __name__ == "__main__":
    run_analysis()
