from __future__ import annotations

import argparse
import json
from pathlib import Path

R_GAS = 8.31446261815324
T_REF = 298.15

REPO_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = Path(__file__).resolve().parent / "output"

FEED_MASS = {
    "water": 0.8094,
    "butanol": 0.1728,
    "NaCl": 0.0054,
    "KCl": 0.0124,
}

REPORTED_PHASE_DATA = {
    "organic": {
        "water": 0.4426,
        "butanol": 0.5570,
        "NaCl": 4.15e-5,
        "KCl": 4.20e-4,
    },
    "aqueous": {
        "water": 0.9627,
        "butanol": 0.0122,
        "NaCl": 0.0076,
        "KCl": 0.0174,
    },
}

TABLE5_LNF_BAR = {
    "water": -3.521,
    "butanol": -5.088,
    "NaCl": -224.891,
    "KCl": -206.733,
}

MW_KG_PER_MOL = {
    "water": 18.01528e-3,
    "butanol": 74.12e-3,
    "NaCl": (22.98976928 + 35.453) * 1.0e-3,
    "KCl": (39.0983 + 35.453) * 1.0e-3,
}

PAPER_TARGETS = {
    "ghat_feed_j_per_mol": -27361.317,
    "ghat_eq_j_per_mol": -27479.860,
}


def _formula_molefractions_from_massfractions(massfractions: dict[str, float]) -> dict[str, float]:
    moles = {name: massfractions[name] / MW_KG_PER_MOL[name] for name in MW_KG_PER_MOL}
    total = sum(moles.values())
    return {name: value / total for name, value in moles.items()}


def _componentwise_mass_split() -> dict[str, float]:
    organic = REPORTED_PHASE_DATA["organic"]
    aqueous = REPORTED_PHASE_DATA["aqueous"]
    return {
        name: (FEED_MASS[name] - aqueous[name]) / (organic[name] - aqueous[name])
        for name in MW_KG_PER_MOL
    }


def _phase_molar_data_if_reported_values_are_massfractions() -> dict[str, object]:
    beta_org = _componentwise_mass_split()["water"]
    beta_aq = 1.0 - beta_org
    moles_org = {
        name: beta_org * REPORTED_PHASE_DATA["organic"][name] / MW_KG_PER_MOL[name]
        for name in MW_KG_PER_MOL
    }
    moles_aq = {
        name: beta_aq * REPORTED_PHASE_DATA["aqueous"][name] / MW_KG_PER_MOL[name]
        for name in MW_KG_PER_MOL
    }
    total_org = sum(moles_org.values())
    total_aq = sum(moles_aq.values())
    alpha_org = total_org / (total_org + total_aq)
    alpha_aq = total_aq / (total_org + total_aq)
    x_org = {name: value / total_org for name, value in moles_org.items()}
    x_aq = {name: value / total_aq for name, value in moles_aq.items()}
    return {
        "beta_org_mass": beta_org,
        "beta_aq_mass": beta_aq,
        "alpha_org_molar": alpha_org,
        "alpha_aq_molar": alpha_aq,
        "x_org_formula_molefrac": x_org,
        "x_aq_formula_molefrac": x_aq,
    }


def _ion_molefractions_from_formula_unit_moles(formula_unit_moles: dict[str, float]) -> tuple[dict[str, float], float]:
    ion_moles = {
        "water": formula_unit_moles["water"],
        "butanol": formula_unit_moles["butanol"],
        "Na+": formula_unit_moles["NaCl"],
        "K+": formula_unit_moles["KCl"],
        "Cl-": formula_unit_moles["NaCl"] + formula_unit_moles["KCl"],
    }
    total = sum(ion_moles.values())
    return ({name: value / total for name, value in ion_moles.items()}, total)


def _formula_basis_ghat(alpha_org: float, alpha_aq: float, x_org: dict[str, float], x_aq: dict[str, float]) -> float:
    sum_org = sum(x_org[name] * TABLE5_LNF_BAR[name] for name in MW_KG_PER_MOL)
    sum_aq = sum(x_aq[name] * TABLE5_LNF_BAR[name] for name in MW_KG_PER_MOL)
    return R_GAS * T_REF * (alpha_org * sum_org + alpha_aq * sum_aq)


def _mass_fraction_ghat(double_count_salts: bool) -> float:
    organic = REPORTED_PHASE_DATA["organic"]
    aqueous = REPORTED_PHASE_DATA["aqueous"]
    beta_org = _componentwise_mass_split()["water"]
    beta_aq = 1.0 - beta_org

    def _phase_sum(phase: dict[str, float]) -> float:
        salt_factor = 2.0 if double_count_salts else 1.0
        return (
            phase["water"] * TABLE5_LNF_BAR["water"]
            + phase["butanol"] * TABLE5_LNF_BAR["butanol"]
            + salt_factor * phase["NaCl"] * TABLE5_LNF_BAR["NaCl"]
            + salt_factor * phase["KCl"] * TABLE5_LNF_BAR["KCl"]
        )

    return R_GAS * T_REF * (beta_org * _phase_sum(organic) + beta_aq * _phase_sum(aqueous))


def _feed_based_reference_ghat(double_count_salts: bool) -> float:
    feed_formula = _formula_molefractions_from_massfractions(FEED_MASS)
    if double_count_salts:
        sum_val = (
            FEED_MASS["water"] * TABLE5_LNF_BAR["water"]
            + FEED_MASS["butanol"] * TABLE5_LNF_BAR["butanol"]
            + 2.0 * FEED_MASS["NaCl"] * TABLE5_LNF_BAR["NaCl"]
            + 2.0 * FEED_MASS["KCl"] * TABLE5_LNF_BAR["KCl"]
        )
    else:
        sum_val = sum(feed_formula[name] * TABLE5_LNF_BAR[name] for name in MW_KG_PER_MOL)
    return R_GAS * T_REF * sum_val


def _general_ionic_basis_ghat_from_mass_consistent_phase_data() -> dict[str, object]:
    beta_org = _componentwise_mass_split()["water"]
    beta_aq = 1.0 - beta_org
    formula_unit_moles_org = {
        name: beta_org * REPORTED_PHASE_DATA["organic"][name] / MW_KG_PER_MOL[name]
        for name in MW_KG_PER_MOL
    }
    formula_unit_moles_aq = {
        name: beta_aq * REPORTED_PHASE_DATA["aqueous"][name] / MW_KG_PER_MOL[name]
        for name in MW_KG_PER_MOL
    }
    x_org, total_org = _ion_molefractions_from_formula_unit_moles(formula_unit_moles_org)
    x_aq, total_aq = _ion_molefractions_from_formula_unit_moles(formula_unit_moles_aq)
    alpha_org = total_org / (total_org + total_aq)
    alpha_aq = total_aq / (total_org + total_aq)

    # For 1:1 salts, the printed mean ionic fugacities imply
    # ln f_Na+ = 2 ln fpm_NaCl - ln f_Cl-
    # ln f_K+  = 2 ln fpm_KCl  - ln f_Cl-
    # Under electroneutrality, the arbitrary chloride gauge cancels identically.
    phase_sum_org = (
        x_org["water"] * TABLE5_LNF_BAR["water"]
        + x_org["butanol"] * TABLE5_LNF_BAR["butanol"]
        + 2.0 * x_org["Na+"] * TABLE5_LNF_BAR["NaCl"]
        + 2.0 * x_org["K+"] * TABLE5_LNF_BAR["KCl"]
    )
    phase_sum_aq = (
        x_aq["water"] * TABLE5_LNF_BAR["water"]
        + x_aq["butanol"] * TABLE5_LNF_BAR["butanol"]
        + 2.0 * x_aq["Na+"] * TABLE5_LNF_BAR["NaCl"]
        + 2.0 * x_aq["K+"] * TABLE5_LNF_BAR["KCl"]
    )
    ghat_eq = R_GAS * T_REF * (alpha_org * phase_sum_org + alpha_aq * phase_sum_aq)

    feed_formula_unit_moles = {
        name: FEED_MASS[name] / MW_KG_PER_MOL[name]
        for name in MW_KG_PER_MOL
    }
    z_feed, _ = _ion_molefractions_from_formula_unit_moles(feed_formula_unit_moles)
    feed_sum = (
        z_feed["water"] * TABLE5_LNF_BAR["water"]
        + z_feed["butanol"] * TABLE5_LNF_BAR["butanol"]
        + 2.0 * z_feed["Na+"] * TABLE5_LNF_BAR["NaCl"]
        + 2.0 * z_feed["K+"] * TABLE5_LNF_BAR["KCl"]
    )
    ghat_feed = R_GAS * T_REF * feed_sum

    return {
        "description": (
            "Use the full ionic convention with components (water, butanol, Na+, K+, Cl-) and the printed mean ionic fugacities. "
            "For the two 1:1 salts, the single-ion chloride reference cancels under electroneutrality, so the ion-basis ghat is unique."
        ),
        "alpha_org_ion_basis": alpha_org,
        "alpha_aq_ion_basis": alpha_aq,
        "x_org_ion_molefrac": x_org,
        "x_aq_ion_molefrac": x_aq,
        "ghat_eq_j_per_mol": ghat_eq,
        "ghat_feed_j_per_mol": ghat_feed,
        "ghat_delta_eq_minus_feed_j_per_mol": ghat_eq - ghat_feed,
        "delta_vs_paper_eq_j_per_mol": ghat_eq - PAPER_TARGETS["ghat_eq_j_per_mol"],
    }


def _least_squares_beta_if_reported_values_are_formula_molefractions() -> dict[str, float]:
    z_feed = _formula_molefractions_from_massfractions(FEED_MASS)
    organic = REPORTED_PHASE_DATA["organic"]
    aqueous = REPORTED_PHASE_DATA["aqueous"]
    feed_vec = [z_feed[name] for name in MW_KG_PER_MOL]
    org_vec = [organic[name] for name in MW_KG_PER_MOL]
    aq_vec = [aqueous[name] for name in MW_KG_PER_MOL]

    diff_org_aq = [org - aq for org, aq in zip(org_vec, aq_vec)]
    diff_feed_aq = [feed - aq for feed, aq in zip(feed_vec, aq_vec)]
    numerator = sum(a * b for a, b in zip(diff_feed_aq, diff_org_aq))
    denominator = sum(a * a for a in diff_org_aq)
    beta = numerator / denominator
    recon = [beta * org + (1.0 - beta) * aq for org, aq in zip(org_vec, aq_vec)]
    residual = [feed - model for feed, model in zip(feed_vec, recon)]
    return {
        "beta_least_squares": beta,
        "max_abs_residual": max(abs(value) for value in residual),
    }


def build_audit() -> dict[str, object]:
    componentwise_beta = _componentwise_mass_split()
    phase_molar_data = _phase_molar_data_if_reported_values_are_massfractions()
    ghat_literal = _formula_basis_ghat(
        phase_molar_data["alpha_org_molar"],
        phase_molar_data["alpha_aq_molar"],
        phase_molar_data["x_org_formula_molefrac"],
        phase_molar_data["x_aq_formula_molefrac"],
    )
    ghat_mass_once = _mass_fraction_ghat(double_count_salts=False)
    ghat_mass_twice = _mass_fraction_ghat(double_count_salts=True)
    ghat_feed_formula = _feed_based_reference_ghat(double_count_salts=False)
    ghat_feed_mass_twice = _feed_based_reference_ghat(double_count_salts=True)
    molefraction_fit = _least_squares_beta_if_reported_values_are_formula_molefractions()
    ionic_basis = _general_ionic_basis_ghat_from_mass_consistent_phase_data()

    cases = {
        "literal_stated_formula": {
            "description": (
                "Interpret the reported phase compositions as mass fractions, convert them to actual phase mole fractions, "
                "use the stated formula with molar phase fractions alpha^(k), and use the printed Table 5 mean ionic fugacities once per salt."
            ),
            "ghat_eq_j_per_mol": ghat_literal,
            "delta_vs_paper_j_per_mol": ghat_literal - PAPER_TARGETS["ghat_eq_j_per_mol"],
        },
        "mass_fraction_sum_once": {
            "description": (
                "Use reported phase compositions directly inside the sum as if they were the x_i values in the printed formula, "
                "with one mean ionic fugacity contribution per salt."
            ),
            "ghat_eq_j_per_mol": ghat_mass_once,
            "delta_vs_paper_j_per_mol": ghat_mass_once - PAPER_TARGETS["ghat_eq_j_per_mol"],
        },
        "mass_fraction_sum_double_count_salts": {
            "description": (
                "Use reported phase compositions directly in the sum and count each salt mean ionic fugacity twice to mimic "
                "a cation-plus-anion contribution."
            ),
            "ghat_eq_j_per_mol": ghat_mass_twice,
            "delta_vs_paper_j_per_mol": ghat_mass_twice - PAPER_TARGETS["ghat_eq_j_per_mol"],
        },
    }

    closest_key = min(cases, key=lambda key: abs(cases[key]["delta_vs_paper_j_per_mol"]))

    return {
        "paper_targets": dict(PAPER_TARGETS),
        "feed_massfractions": dict(FEED_MASS),
        "reported_phase_data": REPORTED_PHASE_DATA,
        "table5_lnf_bar": dict(TABLE5_LNF_BAR),
        "componentwise_mass_split_beta": componentwise_beta,
        "formula_molefraction_fit_if_reported_values_taken_literally": molefraction_fit,
        "molar_data_if_reported_values_are_massfractions": phase_molar_data,
        "feed_formula_unit_molefractions": _formula_molefractions_from_massfractions(FEED_MASS),
        "feed_reference_ghat_formula_basis_j_per_mol": ghat_feed_formula,
        "feed_reference_ghat_mass_basis_double_count_salts_j_per_mol": ghat_feed_mass_twice,
        "general_ionic_reconstruction": ionic_basis,
        "reconstructions": cases,
        "closest_simple_reconstruction_key": closest_key,
        "closest_simple_reconstruction_delta_j_per_mol": cases[closest_key]["delta_vs_paper_j_per_mol"],
        "bottom_line": (
            "The printed case-2 phase compositions behave like mass fractions, not mole fractions, and the printed "
            "ghat cannot be recovered from the printed numbers by the literal stated formula."
        ),
    }


def _format_float(value: float) -> str:
    value = float(value)
    if abs(value) >= 1.0e4:
        return f"{value:.3f}"
    if abs(value) >= 1.0:
        return f"{value:.6f}"
    return f"{value:.6g}"


def _build_markdown_report(audit: dict[str, object]) -> str:
    lines: list[str] = []
    lines.append("# Ascani 2022 Case-2 ghat Audit")
    lines.append("")
    lines.append("## Inputs")
    lines.append("")
    lines.append("| Quantity | Value |")
    lines.append("|---|---:|")
    for key, value in audit["feed_massfractions"].items():
        lines.append(f"| feed {key} mass fraction | {_format_float(value)} |")
    for key, value in audit["paper_targets"].items():
        lines.append(f"| paper {key} | {_format_float(value)} |")
    lines.append("")
    lines.append("## Basis check")
    lines.append("")
    lines.append("| Quantity | Value |")
    lines.append("|---|---:|")
    for key, value in audit["componentwise_mass_split_beta"].items():
        lines.append(f"| beta from {key} mass balance | {_format_float(value)} |")
    lines.append(
        f"| least-squares beta if reported values are treated as formula-unit mole fractions | "
        f"{_format_float(audit['formula_molefraction_fit_if_reported_values_taken_literally']['beta_least_squares'])} |"
    )
    lines.append(
        f"| max residual under literal formula-unit mole-fraction interpretation | "
        f"{_format_float(audit['formula_molefraction_fit_if_reported_values_taken_literally']['max_abs_residual'])} |"
    )
    lines.append("")
    lines.append("## ghat Reconstructions")
    lines.append("")
    lines.append("| Reconstruction | ghat_eq (J/mol) | Delta vs paper (J/mol) |")
    lines.append("|---|---:|---:|")
    for key, payload in audit["reconstructions"].items():
        lines.append(
            f"| {key} | {_format_float(payload['ghat_eq_j_per_mol'])} | {_format_float(payload['delta_vs_paper_j_per_mol'])} |"
        )
    lines.append(f"| feed_reference_formula_basis | {_format_float(audit['feed_reference_ghat_formula_basis_j_per_mol'])} |  |")
    lines.append(
        f"| feed_reference_mass_basis_double_count_salts | "
        f"{_format_float(audit['feed_reference_ghat_mass_basis_double_count_salts_j_per_mol'])} |  |"
    )
    lines.append("")
    lines.append("## General Ionic Convention")
    lines.append("")
    lines.append("| Quantity | Value |")
    lines.append("|---|---:|")
    lines.append(
        f"| ion-basis ghat_eq | {_format_float(audit['general_ionic_reconstruction']['ghat_eq_j_per_mol'])} |"
    )
    lines.append(
        f"| ion-basis ghat_feed | {_format_float(audit['general_ionic_reconstruction']['ghat_feed_j_per_mol'])} |"
    )
    lines.append(
        f"| ion-basis delta (eq-feed) | "
        f"{_format_float(audit['general_ionic_reconstruction']['ghat_delta_eq_minus_feed_j_per_mol'])} |"
    )
    lines.append(
        f"| ion-basis delta vs paper eq | "
        f"{_format_float(audit['general_ionic_reconstruction']['delta_vs_paper_eq_j_per_mol'])} |"
    )
    lines.append("")
    lines.append("## Conclusion")
    lines.append("")
    lines.append(f"- {audit['bottom_line']}")
    lines.append(
        f"- Closest simple reconstruction: `{audit['closest_simple_reconstruction_key']}` "
        f"with delta {_format_float(audit['closest_simple_reconstruction_delta_j_per_mol'])} J/mol."
    )
    return "\n".join(lines) + "\n"


def run_audit(outdir: Path) -> dict[str, Path]:
    audit = build_audit()
    outdir.mkdir(parents=True, exist_ok=True)
    json_path = outdir / "ascani_case2_ghat_audit.json"
    md_path = outdir / "ascani_case2_ghat_audit.md"

    json_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    md_path.write_text(_build_markdown_report(audit), encoding="utf-8")

    print("Ascani case-2 ghat audit")
    print(f"- json: {json_path}")
    print(f"- markdown: {md_path}")
    print(
        "- literal stated formula ghat_eq (J/mol): "
        f"{_format_float(audit['reconstructions']['literal_stated_formula']['ghat_eq_j_per_mol'])}"
    )
    print(
        "- closest simple reconstruction: "
        f"{audit['closest_simple_reconstruction_key']} "
        f"({_format_float(audit['reconstructions'][audit['closest_simple_reconstruction_key']]['ghat_eq_j_per_mol'])} J/mol)"
    )
    return {"json": json_path, "markdown": md_path}


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit Ascani 2022 case-2 ghat from printed paper values.")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=OUTPUT_DIR,
        help="Directory for JSON and Markdown outputs.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    run_audit(args.outdir)


if __name__ == "__main__":
    main()
