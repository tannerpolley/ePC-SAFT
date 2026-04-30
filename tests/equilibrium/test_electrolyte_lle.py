from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest

import epcsaft
from epcsaft import equilibrium as equilibrium_mod
from epcsaft import ePCSAFTMixture
from epcsaft.equilibrium_core.electrolyte_basis import build_electrolyte_basis


REPO_ROOT = Path(__file__).resolve().parents[2]
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


def test_electrolyte_lle_direct_feed_fails_without_predictive_split() -> None:
    feed = np.asarray([0.55, 0.40, 0.025, 0.025], dtype=float)
    mix = _ascani_water_butanol_nacl_mixture(feed)

    with pytest.raises(epcsaft.SolutionError) as excinfo:
        mix.equilibrium(
            kind="electrolyte_lle",
            T=298.15,
            P=1.013e5,
            z=feed,
            options=epcsaft.EquilibriumOptions(include_phase_diagnostics=True),
        )

    diagnostics = excinfo.value.diagnostics
    assert diagnostics["composition_basis"] == "mole_fraction"
    assert diagnostics["equilibrium_route"] == "electrolyte_lle"
    assert diagnostics["variable_model"] == "ascani_transformed_salt_pairs"
    assert diagnostics["basis_rank"] == 1
    assert diagnostics["stability_analysis"] == "electrolyte_tpd"
    assert diagnostics["stability_checked"] is True
    assert diagnostics["algorithm_reference"]["doi"] == "10.1021/acs.jced.1c00866"
    assert diagnostics["phase_equilibrium_model"] == "electrolyte_lle_v4_charge_constrained_solve"
    assert diagnostics["acceptance_gate"] == "predictive_solve_failed"
    assert diagnostics["solver_method"] == "scipy.optimize.least_squares"
    assert diagnostics["tpd_method"] == "scipy.optimize.differential_evolution+minimize"
    assert diagnostics["gibbs_seed_method"] == "scipy.optimize.minimize_scalar"
    assert np.isfinite(diagnostics["solver_residual_norm"])
    assert diagnostics["solver_residual_norm"] >= 0.0
    assert diagnostics["best_failure_reason"]
    assert "v4_partition_seed_api_compatibility" not in json.dumps(diagnostics)
    json.dumps(diagnostics, allow_nan=False)


def test_electrolyte_lle_molality_feed_solves_predictive_split() -> None:
    mix = _ascani_water_butanol_nacl_mixture()

    feed = epcsaft.electrolyte_feed_from_molality(
        mix,
        solvent_feed={"H2O": 0.58, "Butanol": 0.42},
        salt_molality={"NaCl": 1.0},
    )
    result = mix.equilibrium(
        kind="electrolyte_lle",
        T=298.15,
        P=1.013e5,
        solvent_feed={"H2O": 0.58, "Butanol": 0.42},
        salt_molality={"NaCl": 1.0},
        options=epcsaft.EquilibriumOptions(include_phase_diagnostics=True),
    )

    diagnostics = result.diagnostics
    assert result.split_detected is True
    assert len(result.phases) == 2
    np.testing.assert_allclose(diagnostics["feed_composition"], feed, atol=1.0e-15)
    assert diagnostics["composition_basis"] == "molality"
    assert diagnostics["salt_molality"]["NaCl"] == 1.0
    assert diagnostics["acceptance_gate"] == "predictive_nonlinear_solve"
    assert diagnostics["solver_residual_norm"] <= 1.0e-6
    assert diagnostics["material_balance_error"] <= 1.0e-10
    assert diagnostics["charge_balance_error"] <= 1.0e-8
    assert diagnostics["gibbs_delta"] < 0.0
    assert all(float(phase.composition[2]) > 0.0 and float(phase.composition[3]) > 0.0 for phase in result.phases)
    json.dumps(diagnostics, allow_nan=False)


def test_ascani_case2_mixed_salt_solves_without_local_model_fixture() -> None:
    mix = _case2_mixture()

    result = mix.equilibrium(
        kind="electrolyte_lle",
        T=298.15,
        P=1.0e5,
        z=_case2_feed(),
        options=epcsaft.EquilibriumOptions(max_iterations=180, tolerance=1.0e-8, damping=0.5),
    )

    diagnostics = result.diagnostics
    assert result.split_detected is True
    assert len(result.phases) == 2
    assert diagnostics["phase_equilibrium_model"] == "electrolyte_lle_v4_charge_constrained_solve"
    assert diagnostics["basis_rank"] == 2
    assert diagnostics["repeated_stability_iterations"] >= 1
    assert diagnostics["acceptance_gate"] == "predictive_nonlinear_solve"
    assert diagnostics["solver_residual_norm"] <= 1.0e-6
    assert diagnostics["material_balance_error"] <= 1.0e-10
    assert diagnostics["charge_balance_error"] <= 1.0e-8
    assert diagnostics["gibbs_delta"] < 0.0
    phases = {phase.label: phase for phase in result.phases}
    assert phases["aq"].composition[0] > phases["org"].composition[0]
    assert phases["org"].composition[1] > phases["aq"].composition[1]
    assert phases["aq"].composition[2] > phases["org"].composition[2]
    assert phases["aq"].composition[3] > phases["org"].composition[3]
    assert phases["aq"].composition[4] > phases["org"].composition[4]
    payload = json.dumps(diagnostics, allow_nan=False)
    assert "ascani_case2_fixture_regression" not in payload
    assert "v4_partition_seed_api_compatibility" not in payload


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
    assert result.diagnostics["acceptance_gate"] == "predictive_nonlinear_solve"


def test_ascani_counterion_basis_has_expected_rank_and_preserves_charge() -> None:
    mix = _case2_mixture()
    basis = build_electrolyte_basis(mix.species, mix.parameters["z"], _case2_feed())

    assert basis.variable_model == "ascani_transformed_salt_pairs"
    assert basis.rank == 2
    assert basis.e_matrix.shape == (2, 3)
    assert [pair["label"] for pair in basis.salt_pairs] == ["NaCl", "KCl"]

    xi = np.asarray([0.001, 0.002], dtype=float)
    charged_delta = basis.e_matrix.T @ xi
    charged_charges = np.asarray([mix.parameters["z"][i] for i in basis.charged_indices], dtype=float)
    assert abs(float(np.dot(charged_delta, charged_charges))) <= 1.0e-12


def test_electrolyte_stability_payload_is_json_serializable() -> None:
    mix = _case2_mixture()

    result = mix.equilibrium(
        kind="electrolyte_stability",
        T=298.15,
        P=1.0e5,
        z=_case2_feed(),
        options=epcsaft.EquilibriumOptions(max_iterations=80, tolerance=1.0e-8),
    )

    assert result.backend == "electrolyte_tpd"
    assert result.diagnostics["variable_model"] == "ascani_transformed_salt_pairs"
    assert result.diagnostics["basis_rank"] == 2
    assert result.diagnostics["tpd_method"] == "scipy.optimize.differential_evolution+minimize"
    assert result.diagnostics["tpd_objective_value"] == pytest.approx(result.min_tpd)
    assert result.diagnostics["phase_charge_balance"]["trial"] == pytest.approx(0.0, abs=1.0e-8)
    json.dumps(result.to_dict(), allow_nan=False)


def test_one_salt_smoke_reports_predictive_failure_not_partition_acceptance() -> None:
    feed = np.asarray([0.55, 0.40, 0.025, 0.025], dtype=float)
    mix = _ascani_water_butanol_nacl_mixture(feed)

    with pytest.raises(epcsaft.SolutionError) as excinfo:
        mix.equilibrium(
            kind="electrolyte_lle",
            T=298.15,
            P=1.013e5,
            z=feed,
            options=epcsaft.EquilibriumOptions(max_iterations=80, tolerance=1.0e-6),
        )

    diagnostics = excinfo.value.diagnostics
    assert diagnostics["solver_method"] == "scipy.optimize.least_squares"
    assert diagnostics["acceptance_gate"] == "predictive_solve_failed"
    assert np.isfinite(diagnostics["solver_residual_norm"])
    assert diagnostics["solver_residual_norm"] >= 0.0
    assert diagnostics["best_failure_reason"]
    assert "v4_partition_seed_api_compatibility" not in json.dumps(diagnostics)


def test_gibbs_seed_minimization_returns_feasible_seed_and_diagnostics() -> None:
    mix = _case2_mixture()
    feed = _case2_feed()
    opts = epcsaft.EquilibriumOptions(max_iterations=80, tolerance=1.0e-8)
    basis = equilibrium_mod._electrolyte_formula_basis(mix, feed, {})
    stability = equilibrium_mod._electrolyte_stability_from_basis(mix, 298.15, 1.0e5, feed, basis, opts)

    seed = equilibrium_mod._electrolyte_gibbs_seed_from_trial(
        mix,
        298.15,
        1.0e5,
        feed,
        basis,
        stability.trial_composition,
        opts,
    )

    assert seed["seed_name"] == "gibbs_minimized_tpd_trial"
    assert seed["gibbs_seed_method"] == "scipy.optimize.minimize_scalar"
    assert 0.0 < seed["beta_formula"] < 1.0
    assert seed["gibbs_seed_feasible"] is True
    assert np.all(np.asarray(seed["aq_formula"]) > 0.0)
    assert np.all(np.asarray(seed["org_formula"]) > 0.0)
    json.dumps(seed["diagnostics"], allow_nan=False)


def test_predictive_residual_uses_dependent_phase_material_balance() -> None:
    mix = _case2_mixture()
    feed = _case2_feed()
    opts = epcsaft.EquilibriumOptions(max_iterations=80, tolerance=1.0e-8)
    basis = equilibrium_mod._electrolyte_formula_basis(mix, feed, {})
    stability = equilibrium_mod._electrolyte_stability_from_basis(mix, 298.15, 1.0e5, feed, basis, opts)
    seed = equilibrium_mod._electrolyte_gibbs_seed_from_trial(mix, 298.15, 1.0e5, feed, basis, stability.trial_composition, opts)
    variables = equilibrium_mod._pack_predictive_electrolyte_variables(seed["beta_formula"], seed["org_formula"])

    evaluated = equilibrium_mod._evaluate_predictive_electrolyte_variables(mix, 298.15, 1.0e5, feed, basis, variables, opts)

    reconstructed = (1.0 - evaluated["beta_org"]) * evaluated["aq_comp"] + evaluated["beta_org"] * evaluated["org_comp"]
    np.testing.assert_allclose(reconstructed, feed, atol=1.0e-10)
    assert evaluated["material_balance_error"] <= 1.0e-10
    assert evaluated["charge_balance_error"] <= 1.0e-8


def test_electrolyte_lle_solver_failure_reports_json_diagnostics() -> None:
    mix = _case2_mixture()

    with pytest.raises(epcsaft.SolutionError) as excinfo:
        mix.equilibrium(
            kind="electrolyte_lle",
            T=298.15,
            P=1.0e5,
            z=_case2_feed(),
            options=epcsaft.EquilibriumOptions(max_iterations=1, tolerance=1.0e-12),
        )

    diagnostics = excinfo.value.args[1]
    assert diagnostics["phase_equilibrium_model"] == "electrolyte_lle_v4_charge_constrained_solve"
    assert diagnostics["solver_residual_norm"] > 1.0e-6
    json.dumps(diagnostics, allow_nan=False)


def test_electrolyte_lle_rejects_non_neutral_direct_feed() -> None:
    mix = _ascani_water_butanol_nacl_mixture()

    with np.testing.assert_raises_regex(epcsaft.InputError, "charge neutral"):
        mix.equilibrium(kind="electrolyte_lle", T=298.15, P=1.013e5, z=[0.55, 0.40, 0.04, 0.01])


def test_neutral_lle_keeps_rejecting_ionic_mixtures() -> None:
    feed = np.asarray([0.55, 0.40, 0.025, 0.025], dtype=float)
    mix = _ascani_water_butanol_nacl_mixture(feed)

    with np.testing.assert_raises_regex(epcsaft.InputError, "ion-containing"):
        mix.equilibrium(kind="lle_flash", T=298.15, P=1.013e5, z=feed)
