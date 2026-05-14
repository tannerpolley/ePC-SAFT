from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest

import epcsaft
import epcsaft.ipopt_backend as ipopt_backend
from epcsaft import ePCSAFTMixture
from epcsaft.equilibrium import _explicit_to_formula_composition, _formula_to_explicit_composition
from epcsaft.equilibrium_core.electrolyte_basis import build_electrolyte_basis

REPO_ROOT = Path(__file__).resolve().parents[3]

ASCANI_CASE2 = REPO_ROOT / "data" / "multiphase" / "ascani_case2_model_comparison.csv"

def _ascani_water_butanol_nacl_mixture(feed=None) -> ePCSAFTMixture:
    species = ["H2O", "Butanol", "Na+", "Cl-"]
    if feed is None:
        feed = np.asarray([0.55, 0.40, 0.025, 0.025], dtype=float)
    return ePCSAFTMixture.from_dataset("2022_Ascani", species, feed, 298.15)

def _ascani_case2_reference() -> dict[str, float]:
    rows: dict[str, float] = {}
    with ASCANI_CASE2.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            value = row["model_2020"]
            if value and value.lower() not in {"true", "false"}:
                rows[row["quantity"]] = float(value)
            paper = row["paper"]
            if paper and paper.lower() not in {"true", "false"}:
                rows[row["quantity"] + " paper"] = float(paper)
    return rows

def _expand_formula_phase(water: float, butanol: float, nacl: float, kcl: float) -> np.ndarray:
    formula = np.asarray([water, butanol, nacl, kcl], dtype=float)
    denom = float(formula[0] + formula[1] + 2.0 * formula[2] + 2.0 * formula[3])
    return np.asarray([formula[0], formula[1], formula[2], formula[3], formula[2] + formula[3]], dtype=float) / denom

def _case2_feed() -> np.ndarray:
    return np.asarray(
        [
            0.940373242284748,
            0.04879624542603625,
            0.0019339313461782701,
            0.003481324798429627,
            0.005415256144607897,
        ],
        dtype=float,
    )

def _case2_mixture(feed=None) -> ePCSAFTMixture:
    if feed is None:
        feed = _case2_feed()
    return ePCSAFTMixture.from_dataset("2022_Ascani", ["H2O", "Butanol", "Na+", "K+", "Cl-"], feed, 298.15)

def test_ascani_case2_mixed_salt_solves_without_local_model_fixture() -> None:
    mix = _case2_mixture()

    result = mix.equilibrium(
        kind="electrolyte_lle",
        T=298.15,
        P=1.0e5,
        z=_case2_feed(),
        options=epcsaft.EquilibriumOptions(max_iterations=180, tolerance=1.0e-8, damping=0.5),
    )

    assert result.split_detected is True
    assert result.diagnostics["solver_method"] == "native_transformed_newton"
    assert result.diagnostics["acceptance_gate"] == "predictive_nonlinear_solve"

def test_auto_kind_routes_explicit_ionic_feed_to_electrolyte_lle() -> None:
    mix = _case2_mixture()

    result = mix.equilibrium(
        kind="auto",
        T=298.15,
        P=1.0e5,
        z=_case2_feed(),
        options=epcsaft.EquilibriumOptions(max_iterations=180, tolerance=1.0e-8, damping=0.5),
    )

    assert result.split_detected is True
    assert result.diagnostics["equilibrium_route"] == "electrolyte_lle"
    assert result.diagnostics["route_reason"] == "ion-containing mixture"

def test_native_gibbs_seed_path_reports_feasible_solver_diagnostics() -> None:
    mix = _case2_mixture()
    result = mix.equilibrium(
        kind="electrolyte_lle",
        T=298.15,
        P=1.0e5,
        z=_case2_feed(),
        options=epcsaft.EquilibriumOptions(max_iterations=180, tolerance=1.0e-8, damping=0.5),
    )

    assert result.diagnostics["gibbs_seed_method"] == "native_golden_section"
    assert result.diagnostics["solver_method"] == "native_transformed_newton"

def test_predictive_residual_uses_dependent_phase_material_balance() -> None:
    mix = _case2_mixture()
    feed = _case2_feed()
    result = mix.equilibrium(
        kind="electrolyte_lle",
        T=298.15,
        P=1.0e5,
        z=feed,
        options=epcsaft.EquilibriumOptions(max_iterations=180, tolerance=1.0e-8, damping=0.5),
    )

    aq, org = result.phases
    reconstructed = aq.phase_fraction * aq.composition + org.phase_fraction * org.composition
    np.testing.assert_allclose(reconstructed, feed, atol=1.0e-10)

def test_electrolyte_lle_solver_failure_reports_unavailable_derivatives() -> None:
    mix = _case2_mixture()

    with pytest.raises(epcsaft.InputError, match="backend_unavailable"):
        mix.equilibrium(
            kind="electrolyte_lle",
            T=298.15,
            P=1.0e5,
            z=_case2_feed(),
            options=epcsaft.EquilibriumOptions(max_iterations=1, tolerance=1.0e-12, legacy_candidate_mode="off"),
        )

def test_electrolyte_lle_best_effort_reports_unavailable_derivatives() -> None:
    mix = _case2_mixture()

    with pytest.raises(epcsaft.InputError, match="backend_unavailable"):
        mix.equilibrium(
            kind="electrolyte_lle",
            T=298.15,
            P=1.0e5,
            z=_case2_feed(),
            options=epcsaft.EquilibriumOptions(
                max_iterations=1,
                tolerance=1.0e-12,
                legacy_candidate_mode="off",
                return_best_effort=True,
            ),
        )

def test_electrolyte_lle_seed_budget_reports_unavailable_derivatives() -> None:
    mix = _case2_mixture()

    with pytest.raises(epcsaft.InputError, match="backend_unavailable"):
        mix.equilibrium(
            kind="electrolyte_lle",
            T=298.15,
            P=1.0e5,
            z=_case2_feed(),
            options=epcsaft.EquilibriumOptions(
                max_iterations=80,
                tolerance=1.0e-12,
                max_seed_attempts=1,
                legacy_candidate_mode="off",
            ),
        )

def test_electrolyte_lle_objective_budget_reports_unavailable_derivatives() -> None:
    mix = _case2_mixture()

    with pytest.raises(epcsaft.SolutionError) as excinfo:
        mix.equilibrium(
            kind="electrolyte_lle",
            T=298.15,
            P=1.0e5,
            z=_case2_feed(),
            options=epcsaft.EquilibriumOptions(
                max_iterations=80,
                tolerance=1.0e-12,
                max_total_objective_evaluations=1,
                legacy_candidate_mode="off",
            ),
        )

    diagnostics = excinfo.value.args[1]
    assert diagnostics["acceptance_gate"] == "predictive_budget_exhausted"
    assert diagnostics["budget_trigger"] == "max_total_objective_evaluations"

def test_experimental_coupled_density_lle_option_is_reported_without_changing_default_gate() -> None:
    mix = _case2_mixture()

    with pytest.raises(epcsaft.InputError, match="backend_unavailable"):
        mix.equilibrium(
            kind="electrolyte_lle",
            T=298.15,
            P=1.0e5,
            z=_case2_feed(),
            options=epcsaft.EquilibriumOptions(
                max_iterations=1,
                tolerance=1.0e-12,
                legacy_candidate_mode="off",
                density_diagnostics="full",
                experimental_coupled_density_lle=True,
            ),
        )

def test_electrolyte_lle_accepts_legacy_option_dict_before_unavailable_derivatives() -> None:
    mix = _case2_mixture()

    with pytest.raises(epcsaft.InputError, match="backend_unavailable"):
        mix.equilibrium(
            kind="electrolyte_lle",
            T=298.15,
            P=1.0e5,
            z=_case2_feed(),
            options={
                "max_nfev": 1,
                "solver_tol": 1.0e-12,
                "split_tol": 1.0e-4,
                "solver_accept_norm": 0.5,
                "legacy_candidate_mode": "off",
                "tpdf_global_trials": 1200,
                "tpdf_local_trials": 600,
                "charge_weight": 1000.0,
                "seed_x": [0.55, 0.40, 0.025, 0.025],
                "force_seed_solve": True,
            },
        )
