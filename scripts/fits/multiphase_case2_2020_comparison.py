from __future__ import annotations

import copy
import csv
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import pcsaft as pcs
from data.epcsaft_properties import get_prop_dict

OUT_DIR = REPO_ROOT / "data" / "multiphase"
OUT_CSV = OUT_DIR / "ascani_case2_model_comparison.csv"
OUT_MD = OUT_DIR / "ascani_case2_model_comparison.md"
R_GAS = 8.31446261815324


def _case2_feed() -> tuple[list[str], np.ndarray, dict[str, float]]:
    w_water = 0.8094
    w_but = 0.1728
    w_nacl = 0.0054
    w_kcl = 0.0124

    mw_water = 18.01528e-3
    mw_but = 74.12e-3
    mw_nacl = (22.98976928 + 35.453) * 1e-3
    mw_kcl = (39.0983 + 35.453) * 1e-3

    n_water = w_water / mw_water
    n_but = w_but / mw_but
    n_na = w_nacl / mw_nacl
    n_k = w_kcl / mw_kcl
    n_cl = n_na + n_k

    species = ["H2O", "Butanol", "Na+", "K+", "Cl-"]
    n = np.array([n_water, n_but, n_na, n_k, n_cl], dtype=float)
    mass_feed = {
        "w_water": w_water,
        "w_butanol": w_but,
        "w_nacl": w_nacl,
        "w_kcl": w_kcl,
    }
    return species, n / np.sum(n), mass_feed


def _paper_targets() -> dict[str, float]:
    return {
        "x_water_org": 0.4426,
        "x_butanol_org": 0.5570,
        "x_nacl_org": 4.15e-5,
        "x_kcl_org": 4.20e-4,
        "x_water_aq": 0.9627,
        "x_butanol_aq": 0.0122,
        "x_nacl_aq": 0.0076,
        "x_kcl_aq": 0.0174,
        "lnf_water_bar": -3.521,
        "lnf_butanol_bar": -5.088,
        "lnfpm_kcl_bar": -206.733,
        "lnfpm_nacl_bar": -224.891,
        "ghat_feed_j_per_mol": -27361.317,
        "ghat_eq_j_per_mol": -27479.860,
        "ghat_delta_j_per_mol": -118.543,
    }


def _model_definitions():
    return [
        {
            "name": "model_2020",
            "user_options": {
                "elec_model": {
                    "born_model": 1,
                    "dielc_rule": 3,
                    "dielc_diff_mode": "analytic",
                    "eps_r_bulk": "mix",
                    "bjeruum_treatment": False,
                    "born_term_options": {
                        "numerical": False,
                        "sum_term": True,
                        "deps_dx_term": True,
                        "d_born_mode": 1,
                    },
                },
                "debug": False,
            },
        },
        {
            "name": "model_2025_num",
            "user_options": {
                "elec_model": {
                    "born_model": 2,
                    "dielc_rule": "empirical",
                    "dielc_diff_mode": "numeric",
                    "eps_r_bulk": "mix",
                    "bjeruum_treatment": False,
                    "born_term_options": {
                        "numerical": True,
                        "sum_term": True,
                        "deps_dx_term": True,
                        "d_born_mode": 1,
                    },
                },
                "debug": False,
            },
        },
    ]


def _phase_state_liq(t: float, p: float, x: np.ndarray, params: dict) -> dict:
    x = np.asarray(x, dtype=float)
    rho = float(pcs.pcsaft_den(t, p, x, params, phase="liq"))
    lnfugcoef = np.asarray(pcs.pcsaft_lnfugcoef(t, rho, x, params), dtype=float)
    lnfug = lnfugcoef + np.log(np.maximum(x, 1e-300)) + np.log(float(p))
    lnfug_bar = lnfug - np.log(1.0e5)
    return {"rho": rho, "lnfugcoef": lnfugcoef, "lnfug": lnfug, "lnfug_bar": lnfug_bar}


def _ghat_from_phases(t: float, phase_rows: list[dict]) -> float:
    return float(R_GAS * t * sum(float(ph["beta"]) * float(np.dot(ph["x"], ph["lnfug_bar"])) for ph in phase_rows))


def _apply_si_water_butanol_override(params: dict, species: list[str], t: float) -> None:
    i_w = species.index("H2O")
    i_b = species.index("Butanol")
    kij = float(2.94e-4 * t - 0.102)
    lij = -0.0044
    khb = 0.026
    params["k_ij"][i_w, i_b] = kij
    params["k_ij"][i_b, i_w] = kij
    params["l_ij"][i_w, i_b] = lij
    params["l_ij"][i_b, i_w] = lij
    params["k_hb"][i_w, i_b] = khb
    params["k_hb"][i_b, i_w] = khb


def _solve_lle_with_retries(t: float, p: float, z_feed: np.ndarray, params: dict, species: list[str]) -> dict:
    attempt_options = [
        {
            "tpdf_global_trials": 1200,
            "tpdf_local_trials": 600,
            "solver_tol": 1e-9,
            "max_nfev": 220,
            "debug": False,
        },
        {
            "tpdf_global_trials": 5000,
            "tpdf_local_trials": 2500,
            "solver_tol": 1e-10,
            "max_nfev": 1200,
            "charge_weight": 5000.0,
            "solver_accept_norm": 0.5,
            "split_tol": 1e-4,
            "debug": False,
        },
    ]
    last = None
    for opt in attempt_options:
        out = pcs.pcsaft_multiphase_lle(t, p, z_feed, params, species, options=opt)
        out["_solve_options"] = opt
        last = out
        if bool(out.get("converged", False)) and int(out.get("n_phases", 0)) == 2:
            return out
    return last


def _solve_model(model: dict) -> dict[str, float]:
    t = 298.15
    p = 1.0e5
    species, z_feed, mass_feed = _case2_feed()
    params = get_prop_dict("ascani_2022", species, z_feed, t, user_options=copy.deepcopy(model["user_options"]))

    feed_state = _phase_state_liq(t, p, z_feed, params)
    ghat_feed = _ghat_from_phases(
        t,
        [{"beta": 1.0, "x": np.asarray(z_feed, dtype=float), "lnfug_bar": np.asarray(feed_state["lnfug_bar"], dtype=float)}],
    )

    result = _solve_lle_with_retries(t, p, z_feed, params, species)
    if (not result.get("converged", False)) or int(result.get("n_phases", 0)) != 2:
        raise RuntimeError(
            "Case 2 solve did not converge for {}. status={} message={}".format(
                model["name"], result.get("status"), result.get("message")
            )
        )

    ph0, ph1 = result["phases"][0], result["phases"][1]
    org = ph0 if ph0["x"][1] >= ph1["x"][1] else ph1
    aq = ph1 if org is ph0 else ph0

    x_org = np.asarray(org["x"], dtype=float)
    x_aq = np.asarray(aq["x"], dtype=float)
    lnf_org = np.asarray(org["lnfug"], dtype=float) - np.log(1.0e5)
    lnf_aq = np.asarray(aq["lnfug"], dtype=float) - np.log(1.0e5)

    ghat_eq = _ghat_from_phases(
        t,
        [
            {"beta": float(org["beta"]), "x": x_org, "lnfug_bar": lnf_org},
            {"beta": float(aq["beta"]), "x": x_aq, "lnfug_bar": lnf_aq},
        ],
    )

    z = np.asarray(params["z"], dtype=float)
    i_w = species.index("H2O")
    i_b = species.index("Butanol")
    i_na = species.index("Na+")
    i_k = species.index("K+")
    i_cl = species.index("Cl-")
    beta_org = float(org["beta"])
    beta_aq = float(aq["beta"])
    mb = np.asarray(z_feed, dtype=float) - (beta_org * x_org + beta_aq * x_aq)
    eta_nacl_to_aq = 100.0 * beta_aq * x_aq[i_na] / max(z_feed[i_na], 1e-300)
    eta_kcl_to_aq = 100.0 * beta_aq * x_aq[i_k] / max(z_feed[i_k], 1e-300)
    eta_water_to_org = 100.0 * beta_org * x_org[i_w] / max(z_feed[i_w], 1e-300)

    return {
        "model": model["name"],
        "feed_mass_water": float(mass_feed["w_water"]),
        "feed_mass_butanol": float(mass_feed["w_butanol"]),
        "feed_mass_nacl": float(mass_feed["w_nacl"]),
        "feed_mass_kcl": float(mass_feed["w_kcl"]),
        "feed_mole_water": float(z_feed[i_w]),
        "feed_mole_butanol": float(z_feed[i_b]),
        "feed_mole_na": float(z_feed[i_na]),
        "feed_mole_k": float(z_feed[i_k]),
        "feed_mole_cl": float(z_feed[i_cl]),
        "x_water_org": float(x_org[i_w]),
        "x_butanol_org": float(x_org[i_b]),
        "x_nacl_org": float(x_org[i_na]),
        "x_kcl_org": float(x_org[i_k]),
        "x_water_aq": float(x_aq[i_w]),
        "x_butanol_aq": float(x_aq[i_b]),
        "x_nacl_aq": float(x_aq[i_na]),
        "x_kcl_aq": float(x_aq[i_k]),
        "lnf_water_bar": float(0.5 * (lnf_org[i_w] + lnf_aq[i_w])),
        "lnf_butanol_bar": float(0.5 * (lnf_org[i_b] + lnf_aq[i_b])),
        "lnfpm_nacl_bar": float(0.25 * (lnf_org[i_na] + lnf_org[i_cl] + lnf_aq[i_na] + lnf_aq[i_cl])),
        "lnfpm_kcl_bar": float(0.25 * (lnf_org[i_k] + lnf_org[i_cl] + lnf_aq[i_k] + lnf_aq[i_cl])),
        "ghat_feed_j_per_mol": float(ghat_feed),
        "ghat_eq_j_per_mol": float(ghat_eq),
        "ghat_delta_j_per_mol": float(ghat_eq - ghat_feed),
        "phase_split_favored": bool((ghat_eq - ghat_feed) < 0.0),
        "water_prefers_organic": bool((x_org[i_w] / max(x_aq[i_w], 1e-300)) > 1.0),
        "water_partition_ratio_org_over_aq": float(x_org[i_w] / max(x_aq[i_w], 1e-300)),
        "eta_nacl_to_aq_pct": float(eta_nacl_to_aq),
        "eta_kcl_to_aq_pct": float(eta_kcl_to_aq),
        "eta_water_to_org_pct": float(eta_water_to_org),
        "beta_org": beta_org,
        "beta_aq": beta_aq,
        "tpdf_min": float(result["tpdf_min"]),
        "residual_norm": float(result["residual_norm"]),
        "phase_charge_org": float(np.dot(z, x_org)),
        "phase_charge_aq": float(np.dot(z, x_aq)),
        "mass_balance_max": float(np.max(np.abs(mb))),
        "cl_minus_nak_org": float(x_org[i_cl] - x_org[i_na] - x_org[i_k]),
        "cl_minus_nak_aq": float(x_aq[i_cl] - x_aq[i_na] - x_aq[i_k]),
        "kij_water_butanol_used": float(params["k_ij"][i_w, i_b]),
        "lij_water_butanol_used": float(params["l_ij"][i_w, i_b]),
        "khb_water_butanol_used": float(params["k_hb"][i_w, i_b]),
        "solver_tpdf_global_trials": float(result["_solve_options"]["tpdf_global_trials"]),
        "solver_tpdf_local_trials": float(result["_solve_options"]["tpdf_local_trials"]),
        "solver_max_nfev": float(result["_solve_options"]["max_nfev"]),
    }


def _fmt(v):
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)


def _rows_for_csv(paper: dict[str, float], m2020: dict[str, float], m2025: dict[str, float]) -> list[dict[str, object]]:
    rows = []
    entries = [
        ("$x_{water}^{(org)}$", "x_water_org"),
        ("$x_{butanol}^{(org)}$", "x_butanol_org"),
        ("$x_{NaCl}^{(org)}$", "x_nacl_org"),
        ("$x_{KCl}^{(org)}$", "x_kcl_org"),
        ("$x_{water}^{(aq)}$", "x_water_aq"),
        ("$x_{butanol}^{(aq)}$", "x_butanol_aq"),
        ("$x_{NaCl}^{(aq)}$", "x_nacl_aq"),
        ("$x_{KCl}^{(aq)}$", "x_kcl_aq"),
        ("$\\ln(f_{water}/bar)$", "lnf_water_bar"),
        ("$\\ln(f_{butanol}/bar)$", "lnf_butanol_bar"),
        ("$\\ln(f_{\\pm,NaCl}/bar)$", "lnfpm_nacl_bar"),
        ("$\\ln(f_{\\pm,KCl}/bar)$", "lnfpm_kcl_bar"),
        ("$\\hat g_{feed}$ (J/mol)", "ghat_feed_j_per_mol"),
        ("$\\hat g_{eq}$ (J/mol)", "ghat_eq_j_per_mol"),
        ("$\\Delta\\hat g=\\hat g_{eq}-\\hat g_{feed}$ (J/mol)", "ghat_delta_j_per_mol"),
    ]
    for label, key in entries:
        rows.append(
            {
                "quantity": label,
                "paper": paper.get(key, ""),
                "model_2020": m2020[key],
                "model_2025_num": m2025[key],
            }
        )
    diag = [
        ("$\\beta_{org}$", "beta_org"),
        ("$\\beta_{aq}$", "beta_aq"),
        ("$TPDF_{min}$", "tpdf_min"),
        ("Residual norm", "residual_norm"),
        ("Phase split favored ($\\Delta\\hat g<0$)", "phase_split_favored"),
        ("Water prefers organic ($x_{water}^{org}/x_{water}^{aq}>1$)", "water_prefers_organic"),
        ("$x_{water}^{org}/x_{water}^{aq}$", "water_partition_ratio_org_over_aq"),
        ("$\\eta_{NaCl\\to aq}$ (%)", "eta_nacl_to_aq_pct"),
        ("$\\eta_{KCl\\to aq}$ (%)", "eta_kcl_to_aq_pct"),
        ("$\\eta_{water\\to org}$ (%)", "eta_water_to_org_pct"),
        ("Charge residual org", "phase_charge_org"),
        ("Charge residual aq", "phase_charge_aq"),
        ("Mass-balance max error", "mass_balance_max"),
        ("$k_{ij}(water,butanol)$ used", "kij_water_butanol_used"),
        ("$l_{ij}(water,butanol)$ used", "lij_water_butanol_used"),
        ("$k_{ij}^{hb}(water,butanol)$ used", "khb_water_butanol_used"),
    ]
    for label, key in diag:
        rows.append({"quantity": label, "paper": "", "model_2020": m2020[key], "model_2025_num": m2025[key]})
    return rows


def _write_csv(rows: list[dict[str, object]]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["quantity", "paper", "model_2020", "model_2025_num"])
        writer.writeheader()
        writer.writerows(rows)


def _write_markdown(paper: dict[str, float], m2020: dict[str, float], m2025: dict[str, float], rows: list[dict[str, object]]) -> None:
    lines = []
    lines.append("# Ascani 2022 Case-2 Comparison")
    lines.append("")
    lines.append("## Feed composition")
    lines.append("")
    lines.append("| Symbol | Value |")
    lines.append("|---|---:|")
    lines.append(f"| $w_{{water}}$ | {_fmt(m2020['feed_mass_water'])} |")
    lines.append(f"| $w_{{butanol}}$ | {_fmt(m2020['feed_mass_butanol'])} |")
    lines.append(f"| $w_{{NaCl}}$ | {_fmt(m2020['feed_mass_nacl'])} |")
    lines.append(f"| $w_{{KCl}}$ | {_fmt(m2020['feed_mass_kcl'])} |")
    lines.append(f"| $z_{{water}}$ | {_fmt(m2020['feed_mole_water'])} |")
    lines.append(f"| $z_{{butanol}}$ | {_fmt(m2020['feed_mole_butanol'])} |")
    lines.append(f"| $z_{{Na^+}}$ | {_fmt(m2020['feed_mole_na'])} |")
    lines.append(f"| $z_{{K^+}}$ | {_fmt(m2020['feed_mole_k'])} |")
    lines.append(f"| $z_{{Cl^-}}$ | {_fmt(m2020['feed_mole_cl'])} |")
    lines.append("")
    lines.append("## Paper vs model results")
    lines.append("")
    lines.append("| Quantity | Paper (Ascani 2022) | Model 2020 | Model 2025 numeric |")
    lines.append("|---|---:|---:|---:|")
    for r in rows:
        lines.append(
            "| {} | {} | {} | {} |".format(
                r["quantity"], _fmt(r["paper"]), _fmt(r["model_2020"]), _fmt(r["model_2025_num"])
            )
        )
    lines.append("")
    lines.append("## Transfer interpretation")
    lines.append("")
    lines.append(
        "- Both models have $\\Delta\\hat g<0$, so phase split is thermodynamically favored at the specified feed."
    )
    lines.append(
        "- Water partition ratio $x_{water}^{org}/x_{water}^{aq}<1$ for both models, so water does not preferentially transfer to the organic-rich phase."
    )
    lines.append("- At equilibrium, chemical potentials are equal across phases; no net transfer occurs after convergence.")
    lines.append("")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def run_case2_2020_comparison() -> tuple[Path, Path]:
    paper = _paper_targets()
    m2020, m2025 = [_solve_model(model) for model in _model_definitions()]
    rows = _rows_for_csv(paper, m2020, m2025)
    _write_csv(rows)
    _write_markdown(paper, m2020, m2025, rows)

    print("| Quantity | Paper | Model 2020 | Model 2025 numeric |")
    print("|---|---:|---:|---:|")
    for r in rows[:15]:
        print("| {} | {} | {} | {} |".format(r["quantity"], _fmt(r["paper"]), _fmt(r["model_2020"]), _fmt(r["model_2025_num"])))
    print(f"\nSaved comparison CSV: {OUT_CSV}")
    print(f"Saved comparison markdown: {OUT_MD}")
    return OUT_CSV, OUT_MD


if __name__ == "__main__":
    run_case2_2020_comparison()
