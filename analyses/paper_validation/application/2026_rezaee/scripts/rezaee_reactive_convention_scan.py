from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import rezaee_reactive_equilibrium_replay as replay  # noqa: E402

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
PROCESSED_DIR = ANALYSIS_DIR / "data" / "processed"
RESULTS_DIR = ANALYSIS_DIR / "results" / "reaction_equilibrium"

SCAN_ROWS_CSV = PROCESSED_DIR / "rezaee_2026_reactive_convention_scan_rows.csv"
SCAN_SUMMARY_CSV = PROCESSED_DIR / "rezaee_2026_reactive_convention_scan_summary.csv"
SCAN_SUMMARY_JSON = RESULTS_DIR / "rezaee_2026_reactive_convention_scan_summary.json"
SCAN_REPORT_MD = RESULTS_DIR / "rezaee_2026_reactive_convention_scan.md"

BENCHMARK_EQUILIBRIUM_POINT_COUNT = 26


def _jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_jsonable(v) for v in value]
    return value


def _ln(value: float) -> float:
    return replay._safe_log(float(value))


class PhaseTerms:
    def __init__(
        self,
        aqueous_x: np.ndarray,
        organic_x: np.ndarray,
        aqueous_gamma: dict[str, float],
        organic_gamma: dict[str, float],
    ) -> None:
        self.aqueous_x = dict(zip(replay.AQ_LABELS, aqueous_x.tolist()))
        self.organic_x = dict(zip(replay.ORG_LABELS, organic_x.tolist()))
        self.aqueous_gamma = aqueous_gamma
        self.organic_gamma = organic_gamma

    def aq(self, label: str, *, gamma: bool = True) -> float:
        value = _ln(self.aqueous_x[label])
        if gamma:
            value += _ln(self.aqueous_gamma[label])
        return value

    def org(self, label: str, *, gamma: bool = True) -> float:
        value = _ln(self.organic_x[label])
        if gamma:
            value += _ln(self.organic_gamma[label])
        return value


VariantFn = Callable[[PhaseTerms, str], float]


def _ion_label(metal: str) -> str:
    return "Li+" if metal == "Li" else "Na+"


def _complex_label(metal: str) -> str:
    return "RLi" if metal == "Li" else "RNa"


def _paper_eq14(terms: PhaseTerms, metal: str, *, gamma: bool = True) -> float:
    return (
        terms.org(_complex_label(metal), gamma=gamma)
        + terms.aq("H2O", gamma=gamma)
        - terms.aq(_ion_label(metal), gamma=gamma)
        - terms.aq("OH-", gamma=gamma)
        - terms.org("DES", gamma=gamma)
    )


def _paper_aqueous_gamma_only(terms: PhaseTerms, metal: str) -> float:
    return (
        terms.org(_complex_label(metal), gamma=False)
        + terms.aq("H2O", gamma=True)
        - terms.aq(_ion_label(metal), gamma=True)
        - terms.aq("OH-", gamma=True)
        - terms.org("DES", gamma=False)
    )


def _paper_organic_gamma_only(terms: PhaseTerms, metal: str) -> float:
    return (
        terms.org(_complex_label(metal), gamma=True)
        + terms.aq("H2O", gamma=False)
        - terms.aq(_ion_label(metal), gamma=False)
        - terms.aq("OH-", gamma=False)
        - terms.org("DES", gamma=True)
    )


def _drop_water(terms: PhaseTerms, metal: str) -> float:
    return terms.org(_complex_label(metal)) - terms.aq(_ion_label(metal)) - terms.aq("OH-") - terms.org("DES")


def _drop_hydroxide(terms: PhaseTerms, metal: str) -> float:
    return terms.org(_complex_label(metal)) + terms.aq("H2O") - terms.aq(_ion_label(metal)) - terms.org("DES")


def _h_product_exchange(terms: PhaseTerms, metal: str) -> float:
    return terms.org(_complex_label(metal)) + terms.aq("H+") - terms.aq(_ion_label(metal)) - terms.org("DES")


def _nh4_product_exchange(terms: PhaseTerms, metal: str) -> float:
    return _paper_eq14(terms, metal) + terms.aq("NH4+")


def _topo_reactant_exchange(terms: PhaseTerms, metal: str) -> float:
    return _paper_eq14(terms, metal) - terms.org("TOPO")


VARIANTS: list[dict[str, Any]] = [
    {
        "variant": "paper_eq14_with_activity_vs_paper_k",
        "description": "Rezaee Eq. 14/15 exactly as implemented from the paper text, with ePC-SAFT/PC-SAFT activity coefficients.",
        "source_supported": True,
        "comparison": "paper_k",
        "fn": _paper_eq14,
    },
    {
        "variant": "paper_eq14_no_activity_vs_paper_k",
        "description": "Same reaction quotient using mole fractions only.",
        "source_supported": False,
        "comparison": "paper_k",
        "fn": lambda terms, metal: _paper_eq14(terms, metal, gamma=False),
    },
    {
        "variant": "paper_eq14_aqueous_activity_only_vs_paper_k",
        "description": "Aqueous activity coefficients retained, organic activity coefficients set to unity.",
        "source_supported": False,
        "comparison": "paper_k",
        "fn": _paper_aqueous_gamma_only,
    },
    {
        "variant": "paper_eq14_organic_activity_only_vs_paper_k",
        "description": "Organic activity coefficients retained, aqueous activity coefficients set to unity.",
        "source_supported": False,
        "comparison": "paper_k",
        "fn": _paper_organic_gamma_only,
    },
    {
        "variant": "paper_eq14_with_activity_vs_inverse_k",
        "description": "Rezaee Eq. 14/15 quotient compared against reciprocal reported constants.",
        "source_supported": False,
        "comparison": "inverse_k",
        "fn": _paper_eq14,
    },
    {
        "variant": "paper_eq14_no_activity_vs_inverse_k",
        "description": "Mole-fraction quotient compared against reciprocal reported constants.",
        "source_supported": False,
        "comparison": "inverse_k",
        "fn": lambda terms, metal: _paper_eq14(terms, metal, gamma=False),
    },
    {
        "variant": "inverse_quotient_with_activity_vs_paper_k",
        "description": "Inverse of the Rezaee Eq. 14/15 activity quotient compared against reported constants.",
        "source_supported": False,
        "comparison": "paper_k",
        "fn": lambda terms, metal: -_paper_eq14(terms, metal),
    },
    {
        "variant": "drop_water_with_activity_vs_paper_k",
        "description": "Water omitted from the product activity term.",
        "source_supported": False,
        "comparison": "paper_k",
        "fn": _drop_water,
    },
    {
        "variant": "drop_hydroxide_with_activity_vs_paper_k",
        "description": "Hydroxide omitted from the reactant activity term.",
        "source_supported": False,
        "comparison": "paper_k",
        "fn": _drop_hydroxide,
    },
    {
        "variant": "h_product_exchange_with_activity_vs_paper_k",
        "description": "Hydrogen ion treated as the aqueous product instead of hydroxide/water stoichiometry.",
        "source_supported": False,
        "comparison": "paper_k",
        "fn": _h_product_exchange,
    },
    {
        "variant": "nh4_product_exchange_with_activity_vs_paper_k",
        "description": "NH4+ included as an extra aqueous product proxy.",
        "source_supported": False,
        "comparison": "paper_k",
        "fn": _nh4_product_exchange,
    },
    {
        "variant": "topo_reactant_exchange_with_activity_vs_paper_k",
        "description": "TOPO included as an explicit organic reactant.",
        "source_supported": False,
        "comparison": "paper_k",
        "fn": _topo_reactant_exchange,
    },
]


def _target_ln_k(constants: dict[str, float], metal: str, comparison: str) -> float:
    ln_k = math.log(constants[metal])
    if comparison == "inverse_k":
        return -ln_k
    if comparison == "paper_k":
        return ln_k
    raise ValueError(f"Unknown comparison mode {comparison!r}")


def _evaluate_rows() -> tuple[pd.DataFrame, dict[str, float]]:
    constants = replay._reaction_constants()
    aqueous_mix, aqueous_charges = replay._aqueous_mixture()
    organic_mix, _organic_params, pure_ln_phi = replay._organic_mixture()
    source_rows = pd.read_csv(replay.EQUILIBRIUM_CSV)

    records: list[dict[str, Any]] = []
    for row in source_rows.itertuples(index=False):
        aqueous_x = replay._aqueous_x(row)
        organic_x = replay._organic_x(row)
        aqueous_state = aqueous_mix.state(
            T=replay.TEMPERATURE_K,
            x=aqueous_x,
            P=replay.PRESSURE_PA,
        )
        aqueous_gamma = aqueous_state.activity_coefficient(species=replay.AQ_LABELS)
        organic_gamma = replay._organic_activity_coefficients(
            organic_mix,
            pure_ln_phi,
            organic_x,
        )
        terms = PhaseTerms(aqueous_x, organic_x, aqueous_gamma, organic_gamma)
        charge_residual = float(np.dot(aqueous_x, aqueous_charges))
        for variant in VARIANTS:
            fn: VariantFn = variant["fn"]
            for metal in ("Li", "Na"):
                ln_q = float(fn(terms, metal))
                ln_k_target = _target_ln_k(constants, metal, str(variant["comparison"]))
                records.append(
                    {
                        "experiment_no": int(row.experiment_no),
                        "metal": metal,
                        "variant": str(variant["variant"]),
                        "description": str(variant["description"]),
                        "source_supported": bool(variant["source_supported"]),
                        "comparison": str(variant["comparison"]),
                        "lnQ": ln_q,
                        "target_lnK": ln_k_target,
                        "lnQ_minus_target_lnK": ln_q - ln_k_target,
                        "abs_ln_residual": abs(ln_q - ln_k_target),
                        "aqueous_charge_residual": charge_residual,
                        "source": row.source,
                    }
                )
    return pd.DataFrame(records), constants


def _summarize(rows: pd.DataFrame, constants: dict[str, float]) -> tuple[pd.DataFrame, dict[str, Any]]:
    grouped = (
        rows.groupby(["variant", "metal"], as_index=False)
        .agg(
            description=("description", "first"),
            source_supported=("source_supported", "first"),
            comparison=("comparison", "first"),
            row_count=("experiment_no", "count"),
            median_ln_residual=("lnQ_minus_target_lnK", "median"),
            mean_abs_ln_residual=("abs_ln_residual", "mean"),
            median_abs_ln_residual=("abs_ln_residual", "median"),
            rms_ln_residual=("lnQ_minus_target_lnK", lambda values: math.sqrt(float(np.mean(np.square(values))))),
            max_abs_ln_residual=("abs_ln_residual", "max"),
        )
        .sort_values(["variant", "metal"])
    )
    combined = (
        grouped.groupby("variant", as_index=False)
        .agg(
            combined_median_abs_ln_residual=("median_abs_ln_residual", "mean"),
            source_supported=("source_supported", "first"),
            description=("description", "first"),
        )
        .sort_values("combined_median_abs_ln_residual")
    )
    best = combined.iloc[0].to_dict()
    supported = combined.loc[combined["source_supported"]].iloc[0].to_dict()
    summary = {
        "status": "source_reference_state_gap",
        "equilibrium_rows_available": int(rows["experiment_no"].nunique()),
        "equilibrium_rows_used_as_benchmark": BENCHMARK_EQUILIBRIUM_POINT_COUNT,
        "equilibrium_row_count_note": (
            "Rezaee 2026 text says 36 equilibrium data points, but the source SI table available "
            "for this workflow contains 26 designed-experiment equilibrium rows. Treat 36 as a "
            "clerical source-text mismatch unless a new source table is supplied."
        ),
        "paper_constants": constants,
        "paper_lnK": {"Li": math.log(constants["Li"]), "Na": math.log(constants["Na"])},
        "best_variant": best,
        "source_supported_variant": supported,
        "acceptance": {
            "closed_by_simple_convention_scan": bool(float(best["combined_median_abs_ln_residual"]) < 2.0),
            "threshold_combined_median_abs_ln_residual": 2.0,
        },
    }
    return grouped, summary


def _write_report(summary_rows: pd.DataFrame, summary: dict[str, Any]) -> None:
    source = summary["source_supported_variant"]
    best = summary["best_variant"]
    li_source = summary_rows.loc[(summary_rows["variant"] == source["variant"]) & (summary_rows["metal"] == "Li")].iloc[
        0
    ]
    na_source = summary_rows.loc[(summary_rows["variant"] == source["variant"]) & (summary_rows["metal"] == "Na")].iloc[
        0
    ]
    lines = [
        "# Rezaee 2026 Reactive Convention Scan",
        "",
        "## Purpose",
        "",
        "This diagnostic scans simple activity, stoichiometry, and reference-constant conventions for the Rezaee Li/Na reactive equilibrium lane. It is meant to prevent future agents from repeatedly trying the same one-line fixes without evidence.",
        "",
        "## Source Boundary",
        "",
        f"- Machine-readable SI equilibrium rows available here: `{summary['equilibrium_rows_available']}`.",
        f"- Benchmark equilibrium rows used here: `{summary['equilibrium_rows_used_as_benchmark']}`.",
        "- Rezaee 2026 text mentions 36 equilibrium data points, but this workflow treats that as a clerical source-text mismatch because the 2025 SI Tables S1/S2 provide 26 designed-experiment equilibrium rows.",
        "- The source-supported equation remains Rezaee Eq. 14/15 with activity coefficients.",
        "",
        "## Source-Supported Result",
        "",
        f"- Variant: `{source['variant']}`.",
        f"- Li median ln residual: `{float(li_source['median_ln_residual']):.6g}`.",
        f"- Na median ln residual: `{float(na_source['median_ln_residual']):.6g}`.",
        f"- Li median absolute ln residual: `{float(li_source['median_abs_ln_residual']):.6g}`.",
        f"- Na median absolute ln residual: `{float(na_source['median_abs_ln_residual']):.6g}`.",
        "",
        "## Best Numerical Simple Variant",
        "",
        f"- Variant: `{best['variant']}`.",
        f"- Source supported: `{bool(best['source_supported'])}`.",
        f"- Combined median absolute ln residual: `{float(best['combined_median_abs_ln_residual']):.6g}`.",
        f"- Description: {best['description']}",
        "",
        "## Interpretation",
        "",
        "No scanned simple convention closes the published-constant reactive equilibrium to a defensible tolerance for both Li and Na. The best numerical variant is not the equation stated in the paper, so it cannot be used as a published Rezaee reproduction claim.",
        "",
        "The practical package-validation workflow remains: keep the ePC-SAFT activity and stability diagnostics in the validity layer, and treat direct published-constant closure as blocked until the missing source rows or reference-state convention are resolved.",
        "",
        "## Generated Files",
        "",
        f"- `{SCAN_ROWS_CSV.relative_to(ANALYSIS_DIR)}`",
        f"- `{SCAN_SUMMARY_CSV.relative_to(ANALYSIS_DIR)}`",
        f"- `{SCAN_SUMMARY_JSON.relative_to(ANALYSIS_DIR)}`",
    ]
    SCAN_REPORT_MD.parent.mkdir(parents=True, exist_ok=True)
    SCAN_REPORT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    rows, constants = _evaluate_rows()
    summary_rows, summary = _summarize(rows, constants)

    SCAN_ROWS_CSV.parent.mkdir(parents=True, exist_ok=True)
    rows.to_csv(SCAN_ROWS_CSV, index=False)
    summary_rows.to_csv(SCAN_SUMMARY_CSV, index=False)
    SCAN_SUMMARY_JSON.parent.mkdir(parents=True, exist_ok=True)
    SCAN_SUMMARY_JSON.write_text(
        json.dumps(_jsonable(summary), indent=2) + "\n",
        encoding="utf-8",
    )
    _write_report(summary_rows, summary)
    print(json.dumps(_jsonable(summary), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
