from __future__ import annotations

import argparse
import copy
import csv
import json
import sys
from pathlib import Path

import matplotlib
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts._env import require_pcsaft_install

require_pcsaft_install()

import pcsaft as pcs
from pcsaft.parameters import get_prop_dict
from pcsaft.parameters import DATASET_ROOT

matplotlib.use("Agg")
import matplotlib.pyplot as plt


T_REF = 298.15
P_REF = 1.0e5
R_GAS = 8.31446261815324
PARAMETER_DATASET = "ascani_2022"
SPECIES = ["H2O", "Butanol", "Na+", "K+", "Cl-"]
PAPER_TARGETS = {
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
PHASE_TITLES = {
    "organic": "Organic-rich phase",
    "aqueous": "Aqueous-rich phase",
}
RULE_LABELS = {
    0: "constant",
    1: "linear-mole",
    2: "linear-mass",
    3: "combined",
    4: "empirical",
    5: "rule5",
    6: "rule6",
    7: "linear-salt",
    8: "aqueous-organic",
}


def _load_user_options(dataset_name: str) -> dict:
    path = DATASET_ROOT / dataset_name / "user_options.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(payload, dict) and "canonical_user_options" in payload:
        payload = payload["canonical_user_options"]
    if not isinstance(payload, dict):
        raise TypeError(f"user_options.json for '{dataset_name}' did not resolve to a dict.")
    return payload


def _default_model_configs() -> list[dict]:
    return [
        {
            "key": "ascani2022_params_bulow2020_opts",
            "parameter_dataset": PARAMETER_DATASET,
            "options_dataset": "bulow_2020",
            "label": "Ascani 2022 params + current Bulow 2020 options",
            "color": "#2b6cb0",
            "coverage_note": (
                "Pure-component and binary-interaction parameters are fixed to ascani_2022. "
                "Runtime/electrolyte options come from the current bulow_2020 user_options.json, with solvent-specific ion sigma/dispersion precomputes disabled because ascani_2022 only provides pure/any_solvent.csv."
            ),
            "user_options": _load_user_options("bulow_2020"),
        },
        {
            "key": "ascani2022_params_figiel2025_opts",
            "parameter_dataset": PARAMETER_DATASET,
            "options_dataset": "figiel_2025",
            "label": "Ascani 2022 params + current Figiel 2025 options",
            "color": "#8b1e3f",
            "coverage_note": (
                "Pure-component and binary-interaction parameters are fixed to ascani_2022. "
                "Runtime/electrolyte options come from the current figiel_2025 user_options.json, with solvent-specific ion sigma/dispersion precomputes disabled because ascani_2022 only provides pure/any_solvent.csv."
            ),
            "user_options": _load_user_options("figiel_2025"),
        },
    ]


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

    n = np.array([n_water, n_but, n_na, n_k, n_cl], dtype=float)
    mass_feed = {
        "w_water": w_water,
        "w_butanol": w_but,
        "w_nacl": w_nacl,
        "w_kcl": w_kcl,
    }
    return SPECIES, n / np.sum(n), mass_feed


def _phase_state_liq(t: float, p: float, x: np.ndarray, params: dict) -> dict:
    x = np.asarray(x, dtype=float)
    rho = float(pcs.pcsaft_den(t, p, x, params, phase="liq"))
    lnfugcoef = np.asarray(pcs.pcsaft_lnfugcoef(t, rho, x, params), dtype=float)
    lnfug = lnfugcoef + np.log(np.maximum(x, 1e-300)) + np.log(float(p))
    return {
        "rho": rho,
        "lnfugcoef": lnfugcoef,
        "lnfug": lnfug,
        "lnfug_bar": lnfug - np.log(1.0e5),
    }


def _ghat_from_phases(t: float, phase_rows: list[dict]) -> float:
    return float(R_GAS * t * sum(float(ph["beta"]) * float(np.dot(ph["x"], ph["lnfug_bar"])) for ph in phase_rows))


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


def _pair_lnfug_bar(lnfug_bar: np.ndarray, species: list[str], cation: str, anion: str) -> float:
    return float(0.5 * (lnfug_bar[species.index(cation)] + lnfug_bar[species.index(anion)]))


def _format_value(value: object) -> str:
    if isinstance(value, bool):
        return "True" if value else "False"
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def _elec_summary(params: dict) -> str:
    elec = params["elec_model"]
    rel = elec["rel_perm"]
    born = elec["born_model"]
    rule = int(rel["rule"])
    return (
        f"dielc rule {rule} ({RULE_LABELS.get(rule, f'rule{rule}')}); "
        f"rel_perm diff={int(rel['differential_mode'])}; "
        f"Born={'on' if elec['include_born_model'] else 'off'}; "
        f"d_Born_mode={int(born['d_Born_mode'])}; "
        f"shell={bool(born['solvation_shell_model'])}; "
        f"sat={bool(born['dielectric_saturation'])}"
    )


def _build_params_for_config(config: dict, species: list[str], z_feed: np.ndarray) -> dict:
    user_options = copy.deepcopy(config["user_options"])
    user_options["solvated_ion_diameter_mixing_rule"] = False
    user_options["ion_dispersion_mixing_rule"] = False
    params = get_prop_dict(
        config["parameter_dataset"],
        species,
        z_feed,
        T_REF,
        user_options=user_options,
    )
    params["parameter_dataset"] = config["parameter_dataset"]
    params["user_options_dataset"] = config["options_dataset"]
    return params


def _solve_dataset(config: dict) -> dict:
    species, z_feed, mass_feed = _case2_feed()
    params = _build_params_for_config(config, species, z_feed)

    feed_state = _phase_state_liq(T_REF, P_REF, z_feed, params)
    ghat_feed = _ghat_from_phases(
        T_REF,
        [{"beta": 1.0, "x": np.asarray(z_feed, dtype=float), "lnfug_bar": np.asarray(feed_state["lnfug_bar"], dtype=float)}],
    )

    result = _solve_lle_with_retries(T_REF, P_REF, z_feed, params, species)
    if (not result.get("converged", False)) or int(result.get("n_phases", 0)) != 2:
        raise RuntimeError(
            "Case 2 solve did not converge for {}. status={} message={}".format(
                config["key"], result.get("status"), result.get("message")
            )
        )

    ph0, ph1 = result["phases"][0], result["phases"][1]
    org = ph0 if ph0["x"][1] >= ph1["x"][1] else ph1
    aq = ph1 if org is ph0 else ph0

    x_org = np.asarray(org["x"], dtype=float)
    x_aq = np.asarray(aq["x"], dtype=float)
    lnf_org = np.asarray(org["lnfug"], dtype=float) - np.log(1.0e5)
    lnf_aq = np.asarray(aq["lnfug"], dtype=float) - np.log(1.0e5)
    beta_org = float(org["beta"])
    beta_aq = float(aq["beta"])
    n_org = beta_org * x_org
    n_aq = beta_aq * x_aq
    z = np.asarray(params["z"], dtype=float)
    neutral_idx = np.where(np.abs(z) <= 1.0e-12)[0].astype(int)
    charged_idx = np.asarray(result["charged_species_indices"], dtype=int)
    e_matrix = np.asarray(result["e_matrix"], dtype=float)
    lnf_delta = lnf_org - lnf_aq
    neutral_gap = lnf_delta[neutral_idx] if neutral_idx.size else np.zeros(0, dtype=float)
    ionic_gap = e_matrix.dot(lnf_delta[charged_idx]) if charged_idx.size else np.zeros(0, dtype=float)

    ghat_eq = _ghat_from_phases(
        T_REF,
        [
            {"beta": beta_org, "x": x_org, "lnfug_bar": lnf_org},
            {"beta": beta_aq, "x": x_aq, "lnfug_bar": lnf_aq},
        ],
    )

    phase_rows = {}
    for phase_key, beta, x_phase, lnf_phase, n_phase in (
        ("organic", beta_org, x_org, lnf_org, n_org),
        ("aqueous", beta_aq, x_aq, lnf_aq, n_aq),
    ):
        phase_rows[phase_key] = {
            "beta": float(beta),
            "x": {sp: float(x_phase[i]) for i, sp in enumerate(species)},
            "lnfug_bar": {sp: float(lnf_phase[i]) for i, sp in enumerate(species)},
            "share_of_feed_pct": {
                sp: float(100.0 * n_phase[i] / max(z_feed[i], 1e-300))
                for i, sp in enumerate(species)
            },
            "lnfpm_bar": {
                "NaCl": _pair_lnfug_bar(lnf_phase, species, "Na+", "Cl-"),
                "KCl": _pair_lnfug_bar(lnf_phase, species, "K+", "Cl-"),
            },
        }

    mb = np.asarray(z_feed, dtype=float) - (n_org + n_aq)
    i_w = species.index("H2O")
    i_b = species.index("Butanol")
    i_na = species.index("Na+")
    i_k = species.index("K+")
    i_cl = species.index("Cl-")

    return {
        "key": config["key"],
        "parameter_dataset": config["parameter_dataset"],
        "options_dataset": config["options_dataset"],
        "label": config["label"],
        "color": config["color"],
        "coverage_note": config["coverage_note"],
        "elec_summary": _elec_summary(params),
        "feed_mass": dict(mass_feed),
        "feed_z": {sp: float(z_feed[i]) for i, sp in enumerate(species)},
        "phases": phase_rows,
        "paper_compare": {
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
            "beta_org": beta_org,
            "beta_aq": beta_aq,
            "tpdf_min": float(result["tpdf_min"]),
            "residual_norm": float(result["residual_norm"]),
            "phase_charge_org": float(np.dot(z, x_org)),
            "phase_charge_aq": float(np.dot(z, x_aq)),
            "mass_balance_max": float(np.max(np.abs(mb))),
            "neutral_gap_max": float(np.max(np.abs(neutral_gap))) if neutral_gap.size else 0.0,
            "mean_ionic_gap_max": float(np.max(np.abs(ionic_gap))) if ionic_gap.size else 0.0,
            "e_matrix_rank": float(np.linalg.matrix_rank(e_matrix)),
            "mean_ionic_pair_count": float(e_matrix.shape[0]),
            "eta_water_to_org_pct": float(100.0 * n_org[i_w] / max(z_feed[i_w], 1e-300)),
            "eta_butanol_to_org_pct": float(100.0 * n_org[i_b] / max(z_feed[i_b], 1e-300)),
            "eta_na_to_aq_pct": float(100.0 * n_aq[i_na] / max(z_feed[i_na], 1e-300)),
            "eta_k_to_aq_pct": float(100.0 * n_aq[i_k] / max(z_feed[i_k], 1e-300)),
            "eta_cl_to_aq_pct": float(100.0 * n_aq[i_cl] / max(z_feed[i_cl], 1e-300)),
        },
    }


def _summary_rows(results: list[dict]) -> list[dict[str, object]]:
    rows = []
    keys = [
        ("x_water_org", "$x_{water}^{(org)}$"),
        ("x_butanol_org", "$x_{butanol}^{(org)}$"),
        ("x_nacl_org", "$x_{NaCl}^{(org)}$"),
        ("x_kcl_org", "$x_{KCl}^{(org)}$"),
        ("x_water_aq", "$x_{water}^{(aq)}$"),
        ("x_butanol_aq", "$x_{butanol}^{(aq)}$"),
        ("x_nacl_aq", "$x_{NaCl}^{(aq)}$"),
        ("x_kcl_aq", "$x_{KCl}^{(aq)}$"),
        ("lnf_water_bar", "$\\ln(f_{water}/bar)$"),
        ("lnf_butanol_bar", "$\\ln(f_{butanol}/bar)$"),
        ("lnfpm_nacl_bar", "$\\ln(f_{\\pm,NaCl}/bar)$"),
        ("lnfpm_kcl_bar", "$\\ln(f_{\\pm,KCl}/bar)$"),
        ("ghat_feed_j_per_mol", "$\\hat g_{feed}$ (J/mol)"),
        ("ghat_eq_j_per_mol", "$\\hat g_{eq}$ (J/mol)"),
        ("ghat_delta_j_per_mol", "$\\Delta\\hat g$ (J/mol)"),
        ("beta_org", "$\\beta_{org}$"),
        ("beta_aq", "$\\beta_{aq}$"),
        ("tpdf_min", "$TPDF_{min}$"),
        ("residual_norm", "Residual norm"),
        ("phase_charge_org", "Charge residual org"),
        ("phase_charge_aq", "Charge residual aq"),
        ("mass_balance_max", "Mass balance max error"),
        ("neutral_gap_max", "Neutral fugacity gap max"),
        ("mean_ionic_gap_max", "Mean-ionic gap max"),
        ("e_matrix_rank", "E-matrix rank"),
        ("mean_ionic_pair_count", "Independent ionic-pair count"),
        ("eta_water_to_org_pct", "$\\eta_{water\\to org}$ (%)"),
        ("eta_butanol_to_org_pct", "$\\eta_{butanol\\to org}$ (%)"),
        ("eta_na_to_aq_pct", "$\\eta_{Na^+\\to aq}$ (%)"),
        ("eta_k_to_aq_pct", "$\\eta_{K^+\\to aq}$ (%)"),
        ("eta_cl_to_aq_pct", "$\\eta_{Cl^-\\to aq}$ (%)"),
    ]
    for key, label in keys:
        row = {"quantity": label, "paper": PAPER_TARGETS.get(key, "")}
        for result in results:
            row[result["key"]] = result["paper_compare"][key]
        rows.append(row)
    return rows


def _phase_detail_rows(results: list[dict]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for result in results:
        for phase_key, phase_title in PHASE_TITLES.items():
            phase = result["phases"][phase_key]
            for species in SPECIES:
                rows.append(
                    {
                        "model_key": result["key"],
                        "parameter_dataset": result["parameter_dataset"],
                        "options_dataset": result["options_dataset"],
                        "label": result["label"],
                        "phase": phase_key,
                        "phase_title": phase_title,
                        "beta": phase["beta"],
                        "species": species,
                        "mole_fraction": phase["x"][species],
                        "share_of_feed_pct": phase["share_of_feed_pct"][species],
                        "lnfug_bar": phase["lnfug_bar"][species],
                    }
                )
            for salt in ("NaCl", "KCl"):
                rows.append(
                    {
                        "model_key": result["key"],
                        "parameter_dataset": result["parameter_dataset"],
                        "options_dataset": result["options_dataset"],
                        "label": result["label"],
                        "phase": phase_key,
                        "phase_title": phase_title,
                        "beta": phase["beta"],
                        "species": f"{salt}_pm",
                        "mole_fraction": "",
                        "share_of_feed_pct": "",
                        "lnfug_bar": phase["lnfpm_bar"][salt],
                    }
                )
    return rows


def _write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_markdown(path: Path, results: list[dict], summary_rows: list[dict[str, object]]) -> None:
    lines = ["# Ascani Case 2 Dataset Comparison", ""]
    lines.append("## Validation basis")
    lines.append("")
    lines.append(f"- Fixed pure-component and binary-interaction parameter basis: `{PARAMETER_DATASET}`.")
    lines.append("- Only the electrolyte/runtime user options are swapped between the two runs.")
    lines.append("")
    lines.append("## Model presets")
    lines.append("")
    for result in results:
        lines.append(
            f"- `{result['key']}`: params=`{result['parameter_dataset']}`, "
            f"options=`{result['options_dataset']}`, {result['elec_summary']}"
        )
        lines.append(f"- Coverage: {result['coverage_note']}")
    lines.append("")
    lines.append("## Feed composition")
    lines.append("")
    lines.append("| Symbol | Value |")
    lines.append("|---|---:|")
    feed_mass = results[0]["feed_mass"]
    feed_z = results[0]["feed_z"]
    lines.append(f"| $w_{{water}}$ | {_format_value(feed_mass['w_water'])} |")
    lines.append(f"| $w_{{butanol}}$ | {_format_value(feed_mass['w_butanol'])} |")
    lines.append(f"| $w_{{NaCl}}$ | {_format_value(feed_mass['w_nacl'])} |")
    lines.append(f"| $w_{{KCl}}$ | {_format_value(feed_mass['w_kcl'])} |")
    lines.append(f"| $z_{{water}}$ | {_format_value(feed_z['H2O'])} |")
    lines.append(f"| $z_{{butanol}}$ | {_format_value(feed_z['Butanol'])} |")
    lines.append(f"| $z_{{Na^+}}$ | {_format_value(feed_z['Na+'])} |")
    lines.append(f"| $z_{{K^+}}$ | {_format_value(feed_z['K+'])} |")
    lines.append(f"| $z_{{Cl^-}}$ | {_format_value(feed_z['Cl-'])} |")
    lines.append("")
    lines.append("## Paper vs model results")
    lines.append("")
    header = "| Quantity | Paper (Ascani 2022) | " + " | ".join(result["key"] for result in results) + " |"
    divider = "|---|---:|" + "".join("---:|" for _ in results)
    lines.append(header)
    lines.append(divider)
    for row in summary_rows:
        cells = [row["quantity"], _format_value(row["paper"])]
        for result in results:
            cells.append(_format_value(row[result["key"]]))
        lines.append("| " + " | ".join(cells) + " |")
    lines.append("")
    lines.append("## Algorithm checks")
    lines.append("")
    lines.append(
        "- `Neutral fugacity gap max` is the maximum absolute phase-to-phase difference in $\\ln(f_i/bar)$ for neutral species."
    )
    lines.append(
        "- `Mean-ionic gap max` is the maximum absolute residual in $E(\\ln f^{org}-\\ln f^{aq})$ for the charged-species system."
    )
    lines.append(
        "- `E-matrix rank` should match the number of independent ionic-pair equations for the Ascani 2022 construction."
    )
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def _plot_results(results: list[dict], out_path: Path) -> None:
    plt.rcParams.update(
        {
            "font.size": 10,
            "font.family": "DejaVu Serif",
            "axes.linewidth": 1.0,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.2))
    fig.patch.set_facecolor("white")

    paper_bar = {"label": "Paper (Ascani 2022)", "color": "white", "edgecolor": "black", "hatch": "//"}
    width = 0.24

    composition_panels = [
        (
            axes[0, 0],
            "Organic-rich phase composition",
            ["H2O", "Butanol", "NaCl", "KCl"],
            ["x_water_org", "x_butanol_org", "x_nacl_org", "x_kcl_org"],
        ),
        (
            axes[0, 1],
            "Aqueous-rich phase composition",
            ["H2O", "Butanol", "NaCl", "KCl"],
            ["x_water_aq", "x_butanol_aq", "x_nacl_aq", "x_kcl_aq"],
        ),
    ]

    for ax, title, labels, keys in composition_panels:
        x = np.arange(len(labels), dtype=float)
        paper_vals = np.array([float(PAPER_TARGETS[key]) for key in keys], dtype=float)
        panel_vals = [paper_vals]
        ax.bar(
            x - width,
            paper_vals,
            width=width,
            color=paper_bar["color"],
            edgecolor=paper_bar["edgecolor"],
            linewidth=0.9,
            hatch=paper_bar["hatch"],
            label=paper_bar["label"] if ax is axes[0, 0] else None,
        )
        for idx, result in enumerate(results):
            vals = np.array([result["paper_compare"][key] for key in keys], dtype=float)
            panel_vals.append(vals)
            ax.bar(
                x + idx * width,
                vals,
                width=width,
                color=result["color"],
                edgecolor="black",
                linewidth=0.7,
                label=result["label"] if ax is axes[0, 0] else None,
            )
        all_vals = np.concatenate(panel_vals)
        positive = all_vals[all_vals > 0.0]
        ymin = max(float(np.min(positive)) * 0.5, 1e-6)
        ymax = min(float(np.max(all_vals)) * 1.8, 2.0)
        ax.set_yscale("log")
        ax.set_ylim(ymin, ymax)
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_ylabel("mole fraction")
        ax.set_title(title)
        ax.set_facecolor("white")

    fug_ax = axes[1, 0]
    fug_labels = ["H2O", "Butanol", "NaCl_pm", "KCl_pm"]
    fug_keys = ["lnf_water_bar", "lnf_butanol_bar", "lnfpm_nacl_bar", "lnfpm_kcl_bar"]
    x = np.arange(len(fug_labels), dtype=float)
    paper_vals = np.array([float(PAPER_TARGETS[key]) for key in fug_keys], dtype=float)
    fug_ax.bar(
        x - width,
        paper_vals,
        width=width,
        color=paper_bar["color"],
        edgecolor=paper_bar["edgecolor"],
        linewidth=0.9,
        hatch=paper_bar["hatch"],
    )
    for idx, result in enumerate(results):
        vals = np.array([result["paper_compare"][key] for key in fug_keys], dtype=float)
        fug_ax.bar(
            x + idx * width,
            vals,
            width=width,
            color=result["color"],
            edgecolor="black",
            linewidth=0.7,
        )
    fug_all = [paper_vals]
    for result in results:
        fug_all.append(np.array([result["paper_compare"][key] for key in fug_keys], dtype=float))
    fug_stack = np.concatenate(fug_all)
    fug_ax.set_yscale("symlog", linthresh=1.0, linscale=1.0)
    fug_ax.set_ylim(float(np.min(fug_stack)) * 1.15, 1.0)
    fug_ax.set_xticks(x)
    fug_ax.set_xticklabels(fug_labels)
    fug_ax.set_ylabel(r"$\ln(f/bar)$")
    fug_ax.set_title("Equilibrium fugacity comparison")
    fug_ax.set_facecolor("white")
    fug_ax.axhline(0.0, color="black", linewidth=0.9)

    share_ax = axes[1, 1]
    share_labels = ["Na+", "K+", "Cl-"]
    share_keys = ["eta_na_to_aq_pct", "eta_k_to_aq_pct", "eta_cl_to_aq_pct"]
    x = np.arange(len(share_labels), dtype=float)
    offsets = np.linspace(-0.5 * width, 0.5 * width, len(results))
    for idx, result in enumerate(results):
        vals = np.array([result["paper_compare"][key] for key in share_keys], dtype=float)
        share_ax.bar(
            x + offsets[idx],
            vals,
            width=width,
            color=result["color"],
            edgecolor="black",
            linewidth=0.7,
        )
    share_ax.set_xticks(x)
    share_ax.set_xticklabels(share_labels)
    share_ax.set_ylabel("share to aqueous phase (%)")
    share_ax.set_title("Ion partitioning (model only)")
    share_ax.set_facecolor("white")

    legend = axes[0, 0].legend(loc="upper right", fontsize=9, frameon=True)
    legend.get_frame().set_facecolor("white")
    legend.get_frame().set_edgecolor("black")
    legend.get_frame().set_alpha(1.0)

    note_lines = [
        "Hatched bars are paper values; pure/binary parameters stay fixed at ascani_2022 in both runs.",
        "Color changes only reflect the swapped current user-option sets.",
    ]
    fig.suptitle("Ascani case 2: fixed 2022 parameters with swapped runtime options", fontsize=14, y=0.98)
    fig.text(0.5, 0.02, " | ".join(note_lines), ha="center", va="bottom", fontsize=8)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0.02, 0.05, 0.98, 0.95))
    fig.savefig(out_path, dpi=240, bbox_inches="tight")
    plt.close(fig)


def run_analysis(out_dir: Path, model_configs: list[dict] | None = None) -> dict[str, Path]:
    results = [_solve_dataset(config) for config in (model_configs or _default_model_configs())]
    summary_rows = _summary_rows(results)
    detail_rows = _phase_detail_rows(results)

    out_dir.mkdir(parents=True, exist_ok=True)
    summary_csv = out_dir / "ascani_case2_dataset_summary.csv"
    details_csv = out_dir / "ascani_case2_phase_details.csv"
    summary_md = out_dir / "ascani_case2_dataset_summary.md"
    figure_png = out_dir / "ascani_case2_dataset_comparison.png"

    _write_csv(summary_csv, summary_rows, ["quantity", "paper", *[result["key"] for result in results]])
    _write_csv(
        details_csv,
        detail_rows,
        [
            "model_key",
            "parameter_dataset",
            "options_dataset",
            "label",
            "phase",
            "phase_title",
            "beta",
            "species",
            "mole_fraction",
            "share_of_feed_pct",
            "lnfug_bar",
        ],
    )
    _write_markdown(summary_md, results, summary_rows)
    _plot_results(results, figure_png)

    print("Validation basis:")
    print(f"- fixed parameter dataset: {PARAMETER_DATASET}")
    print("- swapped option datasets only:")
    for result in results:
        print(f"  - {result['key']}: options={result['options_dataset']}; {result['elec_summary']}")
        print(f"    coverage: {result['coverage_note']}")
    print("")
    print("| Quantity | Paper | " + " | ".join(result["key"] for result in results) + " |")
    print("|---|---:|" + "".join("---:|" for _ in results))
    for row in summary_rows[:20]:
        cells = [row["quantity"], _format_value(row["paper"])]
        for result in results:
            cells.append(_format_value(row[result["key"]]))
        print("| " + " | ".join(cells) + " |")
    print("")
    print(f"Saved summary CSV: {summary_csv}")
    print(f"Saved phase details CSV: {details_csv}")
    print(f"Saved markdown summary: {summary_md}")
    print(f"Saved figure: {figure_png}")

    return {
        "summary_csv": summary_csv,
        "details_csv": details_csv,
        "summary_md": summary_md,
        "figure_png": figure_png,
    }


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Ascani case 2 fixed-parameter multiphase comparison")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=REPO_ROOT / "scripts" / "multiphase_model_analysis" / "output",
        help="Directory for CSV, markdown, and PNG outputs.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    run_analysis(Path(args.outdir))


if __name__ == "__main__":
    main()

