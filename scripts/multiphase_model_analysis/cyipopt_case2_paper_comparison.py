from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts._env import require_pcsaft_install

require_pcsaft_install()

from scripts.multiphase_model_analysis import ascani_case2_dataset_comparison as case2
from scripts.multiphase_model_analysis.cyipopt_two_phase_experiment import solve_two_phase_lle_cyipopt


def _format_value(value: object) -> str:
    if isinstance(value, bool):
        return "True" if value else "False"
    if isinstance(value, float):
        return f"{value:.6g}"
    if value is None:
        return ""
    return str(value)


def _attempt_options() -> list[dict]:
    return [
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


def _solve_with_backend(backend_key: str, t: float, p: float, z_feed: np.ndarray, params: dict, species: list[str]) -> dict:
    last = None
    for opt in _attempt_options():
        if backend_key == "current":
            out = case2.pcs.pcsaft_multiphase_lle(t, p, z_feed, params, species, options=opt)
        elif backend_key == "cyipopt":
            out = solve_two_phase_lle_cyipopt(t, p, z_feed, params, species, options=opt)
        else:
            raise ValueError(f"Unknown backend_key: {backend_key}")
        out["_solve_options"] = dict(opt)
        last = out
        if bool(out.get("converged", False)) and int(out.get("n_phases", 0)) == 2:
            return out
    if last is None:
        raise RuntimeError(f"No solve executed for backend {backend_key}.")
    return last


def _analyze_solution(config: dict, backend_key: str, backend_label: str, result: dict, params: dict, species: list[str], z_feed: np.ndarray, feed_state: dict, ghat_feed: float, mass_feed: dict) -> dict:
    if (not result.get("converged", False)) or int(result.get("n_phases", 0)) != 2:
        raise RuntimeError(
            f"Case 2 solve did not converge for {config['key']} backend={backend_key}. "
            f"status={result.get('status')} message={result.get('message')}"
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
    ghat_eq = case2._ghat_from_phases(
        case2.T_REF,
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
                "NaCl": case2._pair_lnfug_bar(lnf_phase, species, "Na+", "Cl-"),
                "KCl": case2._pair_lnfug_bar(lnf_phase, species, "K+", "Cl-"),
            },
        }

    mb = np.asarray(z_feed, dtype=float) - (n_org + n_aq)
    i_w = species.index("H2O")
    i_b = species.index("Butanol")
    i_na = species.index("Na+")
    i_k = species.index("K+")
    i_cl = species.index("Cl-")
    paper_compare = {
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
    }
    abs_error = {}
    for key, paper_val in case2.PAPER_TARGETS.items():
        abs_error[key] = abs(float(paper_compare[key]) - float(paper_val))
    return {
        "config_key": config["key"],
        "label": config["label"],
        "backend_key": backend_key,
        "backend_label": backend_label,
        "parameter_dataset": config["parameter_dataset"],
        "options_dataset": config["options_dataset"],
        "elec_summary": case2._elec_summary(params),
        "coverage_note": config["coverage_note"],
        "feed_mass": dict(mass_feed),
        "feed_z": {sp: float(z_feed[i]) for i, sp in enumerate(species)},
        "feed_state": {
            "rho": float(feed_state["rho"]),
            "lnfug_bar": {sp: float(np.asarray(feed_state["lnfug_bar"], dtype=float)[i]) for i, sp in enumerate(species)},
        },
        "phases": phase_rows,
        "paper_compare": paper_compare,
        "paper_abs_error": abs_error,
        "solver_info": result.get("solver_info", {}),
        "solve_options": result.get("_solve_options", {}),
        "status": int(result.get("status", -1)),
        "message": str(result.get("message", "")),
    }


def _collect_results() -> list[dict]:
    rows: list[dict] = []
    for config in case2._default_model_configs():
        species, z_feed, mass_feed = case2._case2_feed()
        params = case2._build_params_for_config(config, species, z_feed)
        feed_state = case2._phase_state_liq(case2.T_REF, case2.P_REF, z_feed, params)
        ghat_feed = case2._ghat_from_phases(
            case2.T_REF,
            [{"beta": 1.0, "x": np.asarray(z_feed, dtype=float), "lnfug_bar": np.asarray(feed_state["lnfug_bar"], dtype=float)}],
        )
        current = _solve_with_backend("current", case2.T_REF, case2.P_REF, z_feed, params, species)
        rows.append(_analyze_solution(config, "current", "Current least_squares solver", current, params, species, z_feed, feed_state, ghat_feed, mass_feed))
        cy_result = _solve_with_backend("cyipopt", case2.T_REF, case2.P_REF, z_feed, params, species)
        rows.append(_analyze_solution(config, "cyipopt", "cyipopt IPOPT experiment", cy_result, params, species, z_feed, feed_state, ghat_feed, mass_feed))
    return rows


SUMMARY_KEYS = [
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
    ("eta_na_to_aq_pct", "$\\eta_{Na^+\\to aq}$ (%)"),
    ("eta_k_to_aq_pct", "$\\eta_{K^+\\to aq}$ (%)"),
    ("eta_cl_to_aq_pct", "$\\eta_{Cl^-\\to aq}$ (%)"),
]


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _phase_detail_rows(results: list[dict]) -> list[dict]:
    rows: list[dict] = []
    for result in results:
        for phase_key, phase_title in case2.PHASE_TITLES.items():
            phase = result["phases"][phase_key]
            for species in case2.SPECIES:
                rows.append(
                    {
                        "config_key": result["config_key"],
                        "backend_key": result["backend_key"],
                        "backend_label": result["backend_label"],
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
                        "config_key": result["config_key"],
                        "backend_key": result["backend_key"],
                        "backend_label": result["backend_label"],
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


def _write_markdown(path: Path, results: list[dict]) -> None:
    lines = ["# Ascani Case 2: cyipopt vs Current Solver vs Paper", ""]
    lines.append("## Validation basis")
    lines.append("")
    lines.append(f"- Fixed pure-component and binary-interaction parameter basis: `{case2.PARAMETER_DATASET}`.")
    lines.append("- Two backends are compared for each option set: the current package `least_squares` solver and the experimental `cyipopt` IPOPT backend.")
    lines.append("- Paper values are the case-2 targets already encoded in the repo and correspond only to the quantities tabulated in the paper basis.")
    lines.append("- Phase-specific fugacity values below are implementation outputs; the paper comparison for fugacity uses the same phase-averaged basis as the existing case-study script.")
    lines.append("")
    lines.append("## Feed composition")
    lines.append("")
    sample = results[0]
    lines.append("| Symbol | Value |")
    lines.append("|---|---:|")
    lines.append(f"| $w_{{water}}$ | {_format_value(sample['feed_mass']['w_water'])} |")
    lines.append(f"| $w_{{butanol}}$ | {_format_value(sample['feed_mass']['w_butanol'])} |")
    lines.append(f"| $w_{{NaCl}}$ | {_format_value(sample['feed_mass']['w_nacl'])} |")
    lines.append(f"| $w_{{KCl}}$ | {_format_value(sample['feed_mass']['w_kcl'])} |")
    lines.append(f"| $z_{{water}}$ | {_format_value(sample['feed_z']['H2O'])} |")
    lines.append(f"| $z_{{butanol}}$ | {_format_value(sample['feed_z']['Butanol'])} |")
    lines.append(f"| $z_{{Na^+}}$ | {_format_value(sample['feed_z']['Na+'])} |")
    lines.append(f"| $z_{{K^+}}$ | {_format_value(sample['feed_z']['K+'])} |")
    lines.append(f"| $z_{{Cl^-}}$ | {_format_value(sample['feed_z']['Cl-'])} |")
    lines.append("")

    grouped: dict[str, list[dict]] = {}
    for result in results:
        grouped.setdefault(result["config_key"], []).append(result)

    for config_key, group in grouped.items():
        group = sorted(group, key=lambda row: row["backend_key"])
        lines.append(f"## {config_key}")
        lines.append("")
        lines.append(f"- Option dataset: `{group[0]['options_dataset']}`")
        lines.append(f"- Electrolyte settings: {group[0]['elec_summary']}")
        lines.append(f"- Coverage note: {group[0]['coverage_note']}")
        lines.append("")
        lines.append("### Paper-comparable summary")
        lines.append("")
        lines.append("| Quantity | Paper | Current | |Error| current | cyipopt | |Error| cyipopt |")
        lines.append("|---|---:|---:|---:|---:|---:|")
        current = next(item for item in group if item["backend_key"] == "current")
        cyipopt = next(item for item in group if item["backend_key"] == "cyipopt")
        for key, label in SUMMARY_KEYS:
            paper_val = case2.PAPER_TARGETS.get(key, "")
            current_val = current["paper_compare"][key]
            cy_val = cyipopt["paper_compare"][key]
            err_current = current["paper_abs_error"].get(key, "")
            err_cy = cyipopt["paper_abs_error"].get(key, "")
            lines.append(
                "| "
                + " | ".join(
                    [
                        label,
                        _format_value(paper_val),
                        _format_value(current_val),
                        _format_value(err_current),
                        _format_value(cy_val),
                        _format_value(err_cy),
                    ]
                )
                + " |"
            )
        lines.append("")
        lines.append("### Solver diagnostics")
        lines.append("")
        lines.append("| Backend | Converged | Status | Message | Nit | Objective | Residual norm |")
        lines.append("|---|---:|---:|---|---:|---:|---:|")
        for item in group:
            info = item.get("solver_info", {})
            lines.append(
                "| "
                + " | ".join(
                    [
                        item["backend_label"],
                        _format_value(True),
                        _format_value(float(item["status"])),
                        _format_value(item["message"]),
                        _format_value(info.get("nit", "")),
                        _format_value(info.get("objective_value", "")),
                        _format_value(item["paper_compare"]["residual_norm"]),
                    ]
                )
                + " |"
            )
        lines.append("")
        lines.append("### Phase-resolved implementation values")
        lines.append("")
        for item in group:
            lines.append(f"#### {item['backend_label']}")
            lines.append("")
            for phase_key in ("organic", "aqueous"):
                phase = item["phases"][phase_key]
                lines.append(f"**{case2.PHASE_TITLES[phase_key]}**")
                lines.append("")
                lines.append(f"- Phase fraction $\\beta$: {_format_value(phase['beta'])}")
                lines.append("| Species | Mole fraction | Share of feed (%) | $\\ln(f/bar)$ |")
                lines.append("|---|---:|---:|---:|")
                for species in case2.SPECIES:
                    lines.append(
                        "| "
                        + " | ".join(
                            [
                                species,
                                _format_value(phase["x"][species]),
                                _format_value(phase["share_of_feed_pct"][species]),
                                _format_value(phase["lnfug_bar"][species]),
                            ]
                        )
                        + " |"
                    )
                lines.append(f"| NaCl_pm |  |  | {_format_value(phase['lnfpm_bar']['NaCl'])} |")
                lines.append(f"| KCl_pm |  |  | {_format_value(phase['lnfpm_bar']['KCl'])} |")
                lines.append("")
        lines.append("### Notes")
        lines.append("")
        lines.append("- The paper-comparable fugacity entries are the phase-averaged neutral and mean-ionic values already used in the repo’s original Ascani case-2 comparison.")
        lines.append("- The per-phase $\\ln(f/bar)$ rows are implementation diagnostics and are not directly tabulated in the paper target set.")
        lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def run_report(out_dir: Path) -> dict[str, Path]:
    results = _collect_results()
    details_csv = out_dir / "cyipopt_case2_phase_details.csv"
    payload_json = out_dir / "cyipopt_case2_paper_comparison.json"
    summary_md = out_dir / "cyipopt_case2_paper_comparison.md"
    _write_csv(
        details_csv,
        _phase_detail_rows(results),
        ["config_key", "backend_key", "backend_label", "phase", "phase_title", "beta", "species", "mole_fraction", "share_of_feed_pct", "lnfug_bar"],
    )
    _write_json(payload_json, results)
    _write_markdown(summary_md, results)
    print(f"Saved phase details CSV: {details_csv}")
    print(f"Saved JSON comparison: {payload_json}")
    print(f"Saved markdown summary: {summary_md}")
    return {"details_csv": details_csv, "payload_json": payload_json, "summary_md": summary_md}


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare current and cyipopt multiphase results against Ascani 2022 case 2.")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=REPO_ROOT / "scripts" / "multiphase_model_analysis" / "output",
        help="Directory for markdown, JSON, and CSV outputs.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    run_report(Path(args.outdir))


if __name__ == "__main__":
    main()
