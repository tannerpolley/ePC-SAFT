"""Fixed-phase electrolyte thermodynamic diagnostics for v4 solver validation."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from epcsaft._types import SolutionError
from epcsaft.epcsaft import ePCSAFTMixture
from epcsaft.equilibrium import EquilibriumOptions, _json_like
from epcsaft.equilibrium_core.electrolyte_basis import build_electrolyte_basis
from scripts.validation.equilibrium_core.phase_state import phase_state

REPO_ROOT = Path(__file__).resolve().parents[3]
KHUIDAIDA_ANALYSIS = REPO_ROOT / "analyses" / "paper_validation" / "application" / "2026_khudaida"
KHUIDAIDA_PARAMS = REPO_ROOT / "data" / "reference" / "epcsaft_parameters" / "2026_Khudaida"
DATASET = "2026_Khudaida"
SPECIES = ["H2O", "Ethanol", "Butanol", "Na+", "Cl-"]
FORMULA_SPECIES = ["H2O", "Ethanol", "Butanol", "NaCl"]
CHARGES = np.asarray([0.0, 0.0, 0.0, 1.0, -1.0], dtype=float)
PRESSURE_PA = 100000.0
RESIDUAL_TOL = 1.0e-6
KHUIDAIDA_CACHED_MATRIX_DIAGNOSTIC_RESIDUAL_ENVELOPE = 5.0e-2


def _fixed_phase_electrolyte_fugacity_residuals(
    aq_comp: np.ndarray,
    org_comp: np.ndarray,
    aq_state: dict[str, Any],
    org_state: dict[str, Any],
    basis: dict[str, Any],
    species: list[str],
) -> tuple[dict[str, float], dict[str, float]]:
    aq_lnf = np.log(aq_comp) + aq_state["ln_phi"]
    org_lnf = np.log(org_comp) + org_state["ln_phi"]
    neutral_residuals = {
        str(species[index]): float(org_lnf[int(index)] - aq_lnf[int(index)]) for index in basis["neutral_indices"]
    }
    salt_residuals = {}
    for pair in basis["salt_pairs"]:
        cation_i = int(pair["cation"])
        anion_i = int(pair["anion"])
        salt_residuals[str(pair["label"])] = float(
            (org_lnf[cation_i] + org_lnf[anion_i]) - (aq_lnf[cation_i] + aq_lnf[anion_i])
        )
    return neutral_residuals, salt_residuals


def _fixed_phase_gibbs_proxy(composition: np.ndarray, state_payload: dict[str, Any]) -> float:
    comp = np.asarray(composition, dtype=float)
    return float(np.sum(comp * (np.log(comp) + np.asarray(state_payload["ln_phi"], dtype=float))))


def load_khudaida_tieline_case(*, figure: int, tie_line: int) -> dict[str, Any]:
    """Load one Khudaida cached tie-line and expand NaCl to explicit ions."""
    figure_dir = _figure_data_dir(figure)
    feed_row = _row_by_tieline(_read_csv(figure_dir / "feed_compositions.csv"), tie_line)
    exp_rows = _phase_rows(_read_csv(figure_dir / "experimental_tielines.csv"), tie_line)
    package_rows = _phase_rows(_read_csv(figure_dir / "model_tielines.csv"), tie_line)

    reference_feed = _feed_row_to_formula(feed_row)
    return _json_like(
        {
            "dataset": DATASET,
            "parameter_source": str(KHUIDAIDA_PARAMS),
            "figure": int(figure),
            "tie_line": int(tie_line),
            "temperature_K": float(feed_row["temperature_K"]),
            "salt_wtfrac": float(feed_row["salt_wtfrac"]),
            "species": list(SPECIES),
            "formula_species": list(FORMULA_SPECIES),
            "reference_feed_formula": reference_feed,
            "reference_feed_composition": _formula_to_explicit(reference_feed),
            "experimental": _source_payload(exp_rows, reference_feed, beta=None),
            "package": _source_payload(package_rows, reference_feed, beta=_phase_beta(package_rows)),
        }
    )


def evaluate_khudaida_tieline(*, figure: int, tie_line: int, source: str) -> dict[str, Any]:
    """Compute current-runtime fixed-phase residual diagnostics for a cached tie-line."""
    case = load_khudaida_tieline_case(figure=figure, tie_line=tie_line)
    phase_payload = _source_from_case(case, source)
    temperature = float(case["temperature_K"])
    mixture = _khudaida_mixture(phase_payload["feed_composition"], temperature)
    feed_comp = np.asarray(phase_payload["feed_composition"], dtype=float)
    org_comp = np.asarray(phase_payload["organic_composition"], dtype=float)
    aq_comp = np.asarray(phase_payload["aqueous_composition"], dtype=float)
    feed_state = _state_payload(mixture, temperature, feed_comp, "feed")
    org_state = _state_payload(mixture, temperature, org_comp, "organic")
    aq_state = _state_payload(mixture, temperature, aq_comp, "aqueous")

    basis = build_electrolyte_basis(
        SPECIES, CHARGES, phase_payload["feed_composition"], salt_labels=("NaCl",)
    ).to_dict()
    neutral_residuals, mean_ionic_residuals = _fixed_phase_electrolyte_fugacity_residuals(
        aq_comp,
        org_comp,
        aq_state,
        org_state,
        basis,
        SPECIES,
    )
    residual_values = [*neutral_residuals.values(), *mean_ionic_residuals.values()]
    residual_norm = float(np.max(np.abs(residual_values))) if residual_values else 0.0

    beta = float(phase_payload["organic_phase_fraction"])
    g_feed = _fixed_phase_gibbs_proxy(feed_comp, feed_state)
    g_org = _fixed_phase_gibbs_proxy(org_comp, org_state)
    g_aq = _fixed_phase_gibbs_proxy(aq_comp, aq_state)
    g_split = float(beta * g_org + (1.0 - beta) * g_aq)
    charge_error = _max_charge_error(phase_payload)
    material_error = _material_balance_error(phase_payload)
    consistent = (
        residual_norm <= RESIDUAL_TOL and material_error <= 1.0e-10 and charge_error <= 1.0e-8 and g_split < g_feed
    )

    return _json_like(
        {
            "dataset": DATASET,
            "source_key": _source_key(source),
            "source": phase_payload["source"],
            "figure": int(figure),
            "tie_line": int(tie_line),
            "temperature_K": temperature,
            "pressure_Pa": PRESSURE_PA,
            "salt_wtfrac": float(case["salt_wtfrac"]),
            "species": list(SPECIES),
            "neutral_residuals": neutral_residuals,
            "mean_ionic_residuals": mean_ionic_residuals,
            "cached_model_residual_norm": phase_payload.get("cached_model_residual_norm"),
            "cached_model_objective": phase_payload.get("cached_model_objective"),
            "fixed_phase_residual_norm": residual_norm,
            "phase_charge_balance": {
                "feed": phase_payload["feed_charge_balance"],
                "organic": phase_payload["organic_charge_balance"],
                "aqueous": phase_payload["aqueous_charge_balance"],
            },
            "charge_balance_error": charge_error,
            "material_balance_error": material_error,
            "reference_feed_material_balance_error": phase_payload["reference_feed_material_balance_error"],
            "gibbs_feed": g_feed,
            "gibbs_organic": g_org,
            "gibbs_aqueous": g_aq,
            "gibbs_split": g_split,
            "gibbs_delta": float(g_split - g_feed),
            "organic_phase_fraction": beta,
            "density_branch": {
                "feed": feed_state["density_branch"],
                "organic": org_state["density_branch"],
                "aqueous": aq_state["density_branch"],
            },
            "fugacity_contribution_terms": {
                "feed": feed_state["fugacity_contribution_terms"],
                "organic": org_state["fugacity_contribution_terms"],
                "aqueous": aq_state["fugacity_contribution_terms"],
            },
            "dielectric_terms": {
                "feed": feed_state["dielectric_terms"],
                "organic": org_state["dielectric_terms"],
                "aqueous": aq_state["dielectric_terms"],
            },
            "basis_rank": basis["basis_rank"],
            "e_matrix": basis["e_matrix"],
            "decision": _fixed_tieline_decision(_source_key(source), consistent),
        }
    )


def compare_khudaida_aad_tables() -> dict[str, Any]:
    """Summarize package-vs-paper AAD gaps from Khudaida Tables 9 and 10."""
    table_dir = _khudaida_analysis_data_dir("tables_9_10")
    rows_compared = []
    max_package = 0.0
    max_paper = 0.0
    package_nonfinite = 0
    for table_name in ("table_9", "table_10"):
        rows = _read_csv(table_dir / f"{table_name}.csv")
        by_temp = {}
        for row in rows:
            key = float(row["T / K"])
            by_temp.setdefault(key, {})[row["Model"]] = _optional_finite_float(row["Grand AAD"])
        for temperature, models in sorted(by_temp.items()):
            package = models["ePC-SAFT (package)"]
            paper = models["ePC-SAFT (paper)"]
            if package is None:
                package_nonfinite += 1
            else:
                max_package = max(max_package, package)
            if paper is not None:
                max_paper = max(max_paper, paper)
            rows_compared.append(
                {
                    "table": table_name,
                    "temperature_K": temperature,
                    "package_grand_aad": package,
                    "paper_epcsaft_grand_aad": paper,
                    "package_minus_paper": package - paper if package is not None and paper is not None else None,
                    "package_aad_state": "finite" if package is not None else "nonfinite_or_rejected",
                }
            )

    return _json_like(
        {
            "dataset": DATASET,
            "tables": ["table_9", "table_10"],
            "rows_compared": len(rows_compared),
            "package_nonfinite_count": package_nonfinite,
            "max_package_grand_aad": max_package,
            "max_paper_epcsaft_grand_aad": max_paper,
            "rows": rows_compared,
            "decision": (
                "package_aad_exceeds_paper_epcsaft_reference"
                if max_package > max_paper
                else "package_aad_does_not_exceed_paper_epcsaft_reference"
            ),
        }
    )


def load_khudaida_digitized_paper_epcsaft(*, figure: int) -> list[dict[str, Any]]:
    """Load a figure's digitized paper ePC-SAFT organic-branch points.

    The imported ``paper_epcsaft_digitized.csv`` files store only the ternary
    salt-free coordinates visible in Figures 2-7. They are therefore a
    salt-free organic-phase plotting/shape reference, not full NaCl phase
    compositions.
    """
    figure_dir = _figure_data_dir(figure)
    raw_rows = _read_csv(figure_dir / "paper_epcsaft_digitized.csv")
    feed_by_tie = {int(row["tie_line"]): row for row in _read_csv(figure_dir / "feed_compositions.csv")}
    rows = []
    for raw in raw_rows:
        point_id = int(raw["point_id"])
        feed_row = feed_by_tie.get(point_id) or next(iter(feed_by_tie.values()))
        salt_free = _normalize_formula(
            np.asarray(
                [
                    float(raw["x_water_salt_free"]),
                    float(raw["x_ethanol_salt_free"]),
                    float(raw["x_isobutanol_salt_free"]),
                ],
                dtype=float,
            )
        )
        rows.append(
            {
                "dataset": DATASET,
                "figure": int(figure),
                "point_id": point_id,
                "tie_line": point_id,
                "temperature_K": float(feed_row["temperature_K"]),
                "salt_wtfrac": float(feed_row["salt_wtfrac"]),
                "salt_free_components": ["H2O", "Ethanol", "Butanol"],
                "salt_free_composition": salt_free,
                "source": raw.get("source", "digitized_user_supplied") or "digitized_user_supplied",
                "series_label": raw.get("series_label", "paper-epcsaft") or "paper-epcsaft",
            }
        )
    return _json_like(rows)


def compare_khudaida_digitized_paper_to_package(*, figure: int) -> dict[str, Any]:
    """Compare native/package organic predictions to digitized paper ePC-SAFT points."""
    figure_dir = _figure_data_dir(figure)
    digitized_rows = load_khudaida_digitized_paper_epcsaft(figure=figure)
    model_csv_rows = _read_csv(figure_dir / "model_tielines.csv")
    model_by_tie = {
        int(row["tie_line"]): row for row in model_csv_rows if str(row.get("phase", "")).strip().lower() == "organic"
    }
    rows = []
    errors = []
    invalid_count = 0
    for digitized in digitized_rows:
        tie_line = int(digitized["tie_line"])
        model_row = model_by_tie.get(tie_line)
        paper = np.asarray(digitized["salt_free_composition"], dtype=float)
        base_row = {
            "figure": int(figure),
            "paper_epcsaft_point_id": int(digitized["point_id"]),
            "paper_epcsaft_tie_line": tie_line,
            "temperature_K": float(digitized["temperature_K"]),
            "salt_wtfrac": float(digitized["salt_wtfrac"]),
            "paper_water_salt_free": float(paper[0]),
            "paper_ethanol_salt_free": float(paper[1]),
            "paper_isobutanol_salt_free": float(paper[2]),
        }
        if model_row is None:
            invalid_count += 1
            rows.append(
                {
                    **base_row,
                    "comparison_state": "package_model_row_absent",
                    "package_converged": None,
                    "package_residual_norm": None,
                }
            )
            continue
        try:
            package = _salt_free_from_formula(_phase_row_to_formula(model_row))
        except ValueError:
            invalid_count += 1
            rows.append(
                {
                    **base_row,
                    "comparison_state": "package_model_row_invalid",
                    "package_converged": _optional_bool(model_row.get("converged")),
                    "package_residual_norm": _optional_finite_float(model_row.get("residual_norm")),
                }
            )
            continue
        abs_error = np.abs(package - paper)
        errors.append(abs_error)
        rows.append(
            {
                **base_row,
                "comparison_state": "compared",
                "package_water_salt_free": float(package[0]),
                "package_ethanol_salt_free": float(package[1]),
                "package_isobutanol_salt_free": float(package[2]),
                "water_abs_error": float(abs_error[0]),
                "ethanol_abs_error": float(abs_error[1]),
                "isobutanol_abs_error": float(abs_error[2]),
                "max_abs_error": float(np.max(abs_error)),
                "aad": float(np.mean(abs_error)),
                "package_converged": _optional_bool(model_row.get("converged")),
                "package_residual_norm": _optional_finite_float(model_row.get("residual_norm")),
            }
        )
    error_arr = np.vstack(errors) if errors else np.empty((0, 3), dtype=float)
    component_aad = np.mean(error_arr, axis=0) if errors else np.full(3, np.nan)
    grand_aad = float(np.mean(error_arr)) if errors else None
    max_abs_error = float(np.max(error_arr)) if errors else None
    source = str(digitized_rows[0]["source"]) if digitized_rows else ""
    package_source = str(next(iter(model_by_tie.values())).get("source", "")) if model_by_tie else ""
    tolerance = 2.0e-2
    return _json_like(
        {
            "dataset": DATASET,
            "figure": int(figure),
            "source": source,
            "package_source": package_source,
            "basis": "organic_phase_salt_free_ternary",
            "digitized_rows": len(digitized_rows),
            "rows_compared": len(rows),
            "finite_rows_compared": len(errors),
            "package_unusable_rows": invalid_count,
            "component_aad": {
                "water": _finite_or_none(component_aad[0]),
                "ethanol": _finite_or_none(component_aad[1]),
                "isobutanol": _finite_or_none(component_aad[2]),
            },
            "organic_salt_free_grand_aad": grand_aad,
            "organic_salt_free_max_abs_error": max_abs_error,
            "rows": rows,
            "decision": (
                "package_matches_digitized_paper_epcsaft_on_salt_free_basis"
                if max_abs_error is not None and max_abs_error <= tolerance
                else "package_differs_from_digitized_paper_epcsaft_on_salt_free_basis"
            ),
        }
    )


def summarize_khudaida_digitized_paper_matrix() -> dict[str, Any]:
    """Summarize package-vs-digitized-paper comparisons over Figures 2-7."""
    figures = [2, 3, 4, 5, 6, 7]
    missing_data: list[str] = []
    comparisons = []
    for figure in figures:
        try:
            comparisons.append(compare_khudaida_digitized_paper_to_package(figure=figure))
        except Exception as exc:  # pragma: no cover - surfaced through diagnostics.
            missing_data.append(f"figure_{figure}: {exc}")
    max_grand = max(
        (
            float(item["organic_salt_free_grand_aad"])
            for item in comparisons
            if item["organic_salt_free_grand_aad"] is not None
        ),
        default=None,
    )
    max_abs = max(
        (
            float(item["organic_salt_free_max_abs_error"])
            for item in comparisons
            if item["organic_salt_free_max_abs_error"] is not None
        ),
        default=None,
    )
    rows_compared = int(sum(int(item["rows_compared"]) for item in comparisons))
    finite_rows_compared = int(sum(int(item["finite_rows_compared"]) for item in comparisons))
    invalid_rows = int(sum(int(item["package_unusable_rows"]) for item in comparisons))
    tolerance = 2.0e-2
    return _json_like(
        {
            "dataset": DATASET,
            "figures": figures,
            "missing_data": missing_data,
            "rows_compared": rows_compared,
            "finite_rows_compared": finite_rows_compared,
            "package_unusable_rows": invalid_rows,
            "max_organic_salt_free_grand_aad": max_grand,
            "max_organic_salt_free_max_abs_error": max_abs,
            "per_figure": comparisons,
            "decision": (
                "package_matches_digitized_paper_epcsaft_on_salt_free_basis"
                if not missing_data and max_abs is not None and max_abs <= tolerance
                else "package_differs_from_digitized_paper_epcsaft_on_salt_free_basis"
            ),
        }
    )


def summarize_khudaida_matrix() -> dict[str, Any]:
    """Summarize fixture completeness and charge neutrality over all Khudaida figures."""
    figures = [2, 3, 4, 5, 6, 7]
    case_count = 0
    missing_data: list[str] = []
    temperatures: set[float] = set()
    salt_wtfracs: set[float] = set()
    max_charge = 0.0
    package_cached_residuals: list[float] = []
    package_cached_residual_rows: list[dict[str, Any]] = []
    package_converged = 0
    package_invalid_count = 0
    digitized_feed_figures: list[int] = []
    digitized_paper_epcsaft_figures: list[int] = []
    for figure in figures:
        figure_dir = _figure_data_dir(figure)
        if (figure_dir / "feed_compositions_digitized.csv").is_file():
            digitized_feed_figures.append(figure)
        if (figure_dir / "paper_epcsaft_digitized.csv").is_file():
            digitized_paper_epcsaft_figures.append(figure)
        try:
            feed_rows = _read_csv(figure_dir / "feed_compositions.csv")
        except Exception as exc:  # pragma: no cover - surfaced through diagnostics.
            missing_data.append(f"figure_{figure}: {exc}")
            continue
        for feed_row in feed_rows:
            tie_line = int(feed_row["tie_line"])
            case_count += 1
            temperatures.add(float(feed_row["temperature_K"]))
            salt_wtfracs.add(float(feed_row["salt_wtfrac"]))
            try:
                exp_rows = _phase_rows(_read_csv(figure_dir / "experimental_tielines.csv"), tie_line)
                reference_feed = _feed_row_to_formula(feed_row)
                experimental = _source_payload(exp_rows, reference_feed, beta=None)
            except Exception as exc:  # pragma: no cover - surfaced through diagnostics.
                missing_data.append(f"figure_{figure} tie_line {tie_line}: {exc}")
                continue
            max_charge = max(
                max_charge,
                abs(float(experimental["feed_charge_balance"])),
                abs(float(experimental["organic_charge_balance"])),
                abs(float(experimental["aqueous_charge_balance"])),
            )
            try:
                package_rows = _phase_rows(_read_csv(figure_dir / "model_tielines.csv"), tie_line)
                package = _source_payload(package_rows, reference_feed, beta=_phase_beta(package_rows))
                max_charge = max(
                    max_charge,
                    abs(float(package["feed_charge_balance"])),
                    abs(float(package["organic_charge_balance"])),
                    abs(float(package["aqueous_charge_balance"])),
                )
                cached_residual = package.get("cached_model_residual_norm")
                if package.get("cached_model_converged"):
                    package_converged += 1
            except ValueError:
                package_invalid_count += 1
                package_rows = _phase_rows(_read_csv(figure_dir / "model_tielines.csv"), tie_line)
                cached_residual = _optional_finite_float(package_rows["organic"].get("residual_norm"))
                if _optional_bool(package_rows["organic"].get("converged")):
                    package_converged += 1
            if cached_residual is not None:
                cached_residual_float = float(cached_residual)
                package_cached_residuals.append(cached_residual_float)
                package_cached_residual_rows.append(
                    {
                        "figure": int(figure),
                        "tie_line": int(tie_line),
                        "residual_norm": cached_residual_float,
                    }
                )

    package_cached_residual_norm_max_case = (
        max(package_cached_residual_rows, key=lambda row: float(row["residual_norm"]))
        if package_cached_residual_rows
        else None
    )
    package_cached_strict_residual_pass_count = int(
        sum(residual <= RESIDUAL_TOL for residual in package_cached_residuals)
    )
    package_cached_diagnostic_residual_over_envelope_count = int(
        sum(residual > KHUIDAIDA_CACHED_MATRIX_DIAGNOSTIC_RESIDUAL_ENVELOPE for residual in package_cached_residuals)
    )

    return _json_like(
        {
            "dataset": DATASET,
            "figures": figures,
            "case_count": case_count,
            "temperatures_K": sorted(temperatures),
            "salt_wtfracs": sorted(salt_wtfracs),
            "max_charge_balance_error": max_charge,
            "missing_data": missing_data,
            "package_cached_converged_count": package_converged,
            "package_invalid_model_count": package_invalid_count,
            "package_cached_residual_norm_max": max(package_cached_residuals) if package_cached_residuals else None,
            "package_cached_residual_norm_min": min(package_cached_residuals) if package_cached_residuals else None,
            "package_cached_residual_norm_max_case": package_cached_residual_norm_max_case,
            "package_cached_strict_residual_pass_count": package_cached_strict_residual_pass_count,
            "package_cached_diagnostic_residual_envelope": KHUIDAIDA_CACHED_MATRIX_DIAGNOSTIC_RESIDUAL_ENVELOPE,
            "package_cached_diagnostic_residual_over_envelope_count": package_cached_diagnostic_residual_over_envelope_count,
            "digitized_feed_figures": digitized_feed_figures,
            "digitized_paper_epcsaft_figures": digitized_paper_epcsaft_figures,
        }
    )


def evaluate_khudaida_solver_gate(*, figure: int, tie_line: int, source: str = "package") -> dict[str, Any]:
    """Try the public v4 solver from a fixed tie-line feed and classify rejection."""
    fixed = evaluate_khudaida_tieline(figure=figure, tie_line=tie_line, source=source)
    case = load_khudaida_tieline_case(figure=figure, tie_line=tie_line)
    phase_payload = _source_from_case(case, source)
    mixture = _khudaida_mixture(phase_payload["feed_composition"], float(case["temperature_K"]))
    solver_diagnostics: dict[str, Any]
    try:
        result = mixture.equilibrium(
            kind="electrolyte_lle",
            T=float(case["temperature_K"]),
            P=PRESSURE_PA,
            z=phase_payload["feed_composition"],
            options=EquilibriumOptions(max_iterations=40),
        )
    except SolutionError as exc:
        solver_outcome = "rejected"
        solver_diagnostics = dict(getattr(exc, "diagnostics", None) or {})
        acceptance_gate = str(solver_diagnostics.get("acceptance_gate", "predictive_solve_failed"))
    else:
        solver_outcome = "accepted"
        solver_diagnostics = dict(result.diagnostics)
        acceptance_gate = str(solver_diagnostics.get("acceptance_gate", "predictive_nonlinear_solve"))

    if fixed["fixed_phase_residual_norm"] > RESIDUAL_TOL:
        decision = "current_surface_inconsistent_before_solver"
    elif solver_outcome == "accepted":
        decision = "solver_accepts_package_fixed_tieline_feed"
    else:
        decision = "fixed_tieline_consistent_solver_suspect"

    return _json_like(
        {
            "dataset": DATASET,
            "source_key": _source_key(source),
            "source": fixed["source"],
            "figure": int(figure),
            "tie_line": int(tie_line),
            "fixed_phase_residual_norm": fixed["fixed_phase_residual_norm"],
            "gibbs_delta": fixed["gibbs_delta"],
            "solver_outcome": solver_outcome,
            "acceptance_gate": acceptance_gate,
            "solver_diagnostics": _sanitize_solver_diagnostics(solver_diagnostics),
            "decision": decision,
        }
    )


def _figure_data_dir(figure: int) -> Path:
    return _khudaida_analysis_data_dir(f"figure_{int(figure)}")


def _khudaida_analysis_data_dir(result_set: str) -> Path:
    candidates = [
        KHUIDAIDA_ANALYSIS / "figures" / result_set / "output" / "data",
    ]
    for path in candidates:
        if path.is_dir():
            return path
    raise FileNotFoundError(f"Khudaida analysis data directory is missing: {candidates[0]}")


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"Khudaida diagnostic data file is missing: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def _row_by_tieline(rows: list[dict[str, str]], tie_line: int) -> dict[str, str]:
    for row in rows:
        if int(row["tie_line"]) == int(tie_line):
            return row
    raise KeyError(f"Could not find Khudaida tie-line {tie_line}.")


def _phase_rows(rows: list[dict[str, str]], tie_line: int) -> dict[str, dict[str, str]]:
    selected = {row["phase"].strip().lower(): row for row in rows if int(row["tie_line"]) == int(tie_line)}
    if set(selected) != {"organic", "aqueous"}:
        raise KeyError(f"Expected organic and aqueous rows for Khudaida tie-line {tie_line}.")
    return selected


def _phase_beta(rows: dict[str, dict[str, str]]) -> float:
    return float(rows["organic"].get("beta") or 0.5)


def _optional_float(value: Any) -> float | None:
    if value in (None, ""):
        return None
    return float(value)


def _finite_or_none(value: Any) -> float | None:
    numeric = float(value)
    return numeric if np.isfinite(numeric) else None


def _optional_finite_float(value: Any) -> float | None:
    if value in (None, ""):
        return None
    return _finite_or_none(value)


def _optional_bool(value: Any) -> bool | None:
    if value in (None, ""):
        return None
    return str(value).strip().lower() in {"true", "1", "yes"}


def _source_payload(
    rows: dict[str, dict[str, str]], reference_feed_formula: np.ndarray, beta: float | None
) -> dict[str, Any]:
    organic_formula = _phase_row_to_formula(rows["organic"])
    aqueous_formula = _phase_row_to_formula(rows["aqueous"])
    organic = _formula_to_explicit(organic_formula)
    aqueous = _formula_to_explicit(aqueous_formula)
    reference_feed = _formula_to_explicit(reference_feed_formula)
    phase_fraction = _best_phase_fraction(reference_feed, organic, aqueous) if beta is None else float(beta)
    phase_fraction = float(np.clip(phase_fraction, 1.0e-12, 1.0 - 1.0e-12))
    feed = phase_fraction * organic + (1.0 - phase_fraction) * aqueous
    return {
        "source": rows["organic"].get("source", ""),
        "organic_phase_fraction": phase_fraction,
        "feed_composition": feed,
        "reference_feed_composition": reference_feed,
        "organic_composition": organic,
        "aqueous_composition": aqueous,
        "organic_formula": organic_formula,
        "aqueous_formula": aqueous_formula,
        "feed_charge_balance": _charge(feed),
        "organic_charge_balance": _charge(organic),
        "aqueous_charge_balance": _charge(aqueous),
        "reference_feed_material_balance_error": float(np.max(np.abs(reference_feed - feed))),
        "cached_model_residual_norm": _optional_float(rows["organic"].get("residual_norm")),
        "cached_model_objective": _optional_float(rows["organic"].get("objective")),
        "cached_model_converged": _optional_bool(rows["organic"].get("converged")),
    }


def _feed_row_to_formula(row: dict[str, str]) -> np.ndarray:
    return _normalize_formula(
        np.asarray(
            [
                float(row["x_water_total"]),
                float(row["x_ethanol_total"]),
                float(row["x_isobutanol_total"]),
                float(row["x_nacl_total"]),
            ],
            dtype=float,
        )
    )


def _phase_row_to_formula(row: dict[str, str]) -> np.ndarray:
    return _normalize_formula(
        np.asarray(
            [
                float(row["x_water"]),
                float(row["x_ethanol"]),
                float(row["x_isobutanol"]),
                float(row["x_nacl"]),
            ],
            dtype=float,
        )
    )


def _salt_free_from_formula(row: np.ndarray) -> np.ndarray:
    formula = _normalize_formula(np.asarray(row, dtype=float))
    return _normalize_formula(formula[:3])


def _normalize_formula(value: np.ndarray) -> np.ndarray:
    total = float(np.sum(value))
    if not np.isfinite(total) or total <= 0.0:
        raise ValueError("Khudaida formula composition has non-positive total.")
    return np.asarray(value, dtype=float) / total


def _formula_to_explicit(formula: np.ndarray) -> np.ndarray:
    x = _normalize_formula(np.asarray(formula, dtype=float))
    expanded = np.asarray([x[0], x[1], x[2], x[3], x[3]], dtype=float)
    return expanded / float(np.sum(expanded))


def _best_phase_fraction(reference_feed: np.ndarray, organic: np.ndarray, aqueous: np.ndarray) -> float:
    direction = organic - aqueous
    denom = float(np.dot(direction, direction))
    if denom <= 0.0:
        return 0.5
    beta = float(np.dot(reference_feed - aqueous, direction) / denom)
    return float(np.clip(beta, 1.0e-12, 1.0 - 1.0e-12))


def _source_from_case(case: dict[str, Any], source: str) -> dict[str, Any]:
    key = _source_key(source)
    if key not in case:
        raise KeyError("Khudaida source must be 'experimental' or 'package'.")
    return case[key]


def _source_key(source: str) -> str:
    token = str(source).strip().lower()
    if token in {"experimental", "experiment", "paper", "paper_table"}:
        return "experimental"
    if token in {"package", "model", "epcsaft_package"}:
        return "package"
    raise KeyError("Khudaida source must be 'experimental' or 'package'.")


def _khudaida_mixture(feed: Any, temperature: float) -> ePCSAFTMixture:
    return ePCSAFTMixture.from_dataset(
        DATASET,
        SPECIES,
        np.asarray(feed, dtype=float),
        float(temperature),
        user_options=_khudaida_user_options(),
    )


def _khudaida_user_options() -> dict[str, Any]:
    options_path = KHUIDAIDA_PARAMS / "user_options.json"
    with options_path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _state_payload(
    mixture: ePCSAFTMixture,
    temperature: float,
    composition: np.ndarray,
    label: str,
) -> dict[str, Any]:
    payload = phase_state(
        mixture,
        float(temperature),
        PRESSURE_PA,
        composition,
        "liq",
        f"khudaida_{label}_diagnostic",
        include_diagnostics=True,
    )
    state = payload["state"]
    terms = state.fugacity_coefficient(return_contribution_terms=True)["terms"]
    epsilon, epsilon_derivatives = state.relative_permittivity()
    payload["density_branch"] = {
        "phase": "liq",
        "density_molar": float(state.molar_density()),
        "mass_density": float(state.mass_density()),
        "pressure": float(state.pressure()),
    }
    payload["fugacity_contribution_terms"] = {
        name: _json_like(np.asarray(value, dtype=float))
        for name, value in terms.items()
        if name in {"hc", "disp", "assoc", "ion", "born"}
    }
    payload["dielectric_terms"] = {
        "relative_permittivity": float(epsilon),
        "composition_derivatives": _json_like(np.asarray(epsilon_derivatives, dtype=float)),
    }
    return payload


def _material_balance_error(payload: dict[str, Any]) -> float:
    feed = np.asarray(payload["feed_composition"], dtype=float)
    organic = np.asarray(payload["organic_composition"], dtype=float)
    aqueous = np.asarray(payload["aqueous_composition"], dtype=float)
    beta = float(payload["organic_phase_fraction"])
    return float(np.max(np.abs(feed - (beta * organic + (1.0 - beta) * aqueous))))


def _max_charge_error(payload: dict[str, Any]) -> float:
    return float(
        max(
            abs(float(payload["feed_charge_balance"])),
            abs(float(payload["organic_charge_balance"])),
            abs(float(payload["aqueous_charge_balance"])),
        )
    )


def _charge(composition: Any) -> float:
    return float(np.dot(np.asarray(composition, dtype=float), CHARGES))


def _fixed_tieline_decision(source: str, consistent: bool) -> str:
    if source == "package":
        return (
            "package_fixed_tieline_internally_consistent"
            if consistent
            else "no_fixed_tieline_satisfies_current_surface"
        )
    return (
        "experimental_tieline_matches_current_surface"
        if consistent
        else "thermodynamic_surface_differs_from_reference_tieline"
    )


def _sanitize_solver_diagnostics(diagnostics: dict[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key in (
        "acceptance_gate",
        "rejection_reason",
        "solver_residual_norm",
        "gibbs_delta",
        "material_balance_error",
        "charge_balance_error",
        "phase_distance",
        "equilibrium_route",
        "phase_equilibrium_model",
        "phase_labels_swapped",
        "phase_label_basis",
    ):
        if key in diagnostics:
            out[key] = diagnostics[key]
    return out
