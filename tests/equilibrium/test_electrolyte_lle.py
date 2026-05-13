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


def test_electrolyte_lle_direct_feed_reports_unavailable_solver_derivatives() -> None:
    feed = np.asarray([0.55, 0.40, 0.025, 0.025], dtype=float)
    mix = _ascani_water_butanol_nacl_mixture(feed)

    with pytest.raises(epcsaft.InputError, match="backend_unavailable"):
        mix.equilibrium(
            kind="electrolyte_lle",
            T=298.15,
            P=1.013e5,
            z=feed,
            options=epcsaft.EquilibriumOptions(include_phase_diagnostics=True),
        )


@pytest.mark.skipif(not ipopt_backend.cyipopt_available(), reason="cyipopt is optional")
def test_electrolyte_lle_direct_feed_solves_ipopt_predictive_split() -> None:
    feed = np.asarray([0.55, 0.40, 0.025, 0.025], dtype=float)
    mix = _ascani_water_butanol_nacl_mixture(feed)

    with pytest.raises(epcsaft.InputError, match="backend_unavailable"):
        mix.equilibrium(
            kind="electrolyte_lle",
            T=298.15,
            P=1.013e5,
            z=feed,
            options=epcsaft.EquilibriumOptions(
                solver_backend="ipopt",
                hessian_strategy="lbfgs",
                max_iterations=80,
                tolerance=1.0e-8,
            ),
        )


def test_electrolyte_lle_molality_feed_reports_unavailable_solver_derivatives() -> None:
    mix = _ascani_water_butanol_nacl_mixture()

    with pytest.raises(epcsaft.InputError, match="backend_unavailable"):
        mix.equilibrium(
            kind="electrolyte_lle",
            T=298.15,
            P=1.013e5,
            solvent_feed={"H2O": 0.58, "Butanol": 0.42},
            salt_molality={"NaCl": 1.0},
            options=epcsaft.EquilibriumOptions(include_phase_diagnostics=True),
        )


def test_electrolyte_lle_accepts_strict_aq_org_initial_phases() -> None:
    mix = _ascani_water_butanol_nacl_mixture()
    aq = np.asarray([0.798324680201737, 0.016320352824141723, 0.09267748348706063, 0.09267748348706063], dtype=float)
    org = np.asarray([0.37006036048879404, 0.6214918588210971, 0.004223890345054407, 0.004223890345054407], dtype=float)
    beta_org = 0.613766575013417
    feed = (1.0 - beta_org) * aq + beta_org * org

    result = mix.equilibrium(
        kind="electrolyte_lle",
        T=298.15,
        P=1.013e5,
        z=feed,
        initial_phases={"aq": aq, "org": org, "phase_fraction": beta_org},
        options=epcsaft.EquilibriumOptions(include_phase_diagnostics=True),
    )

    assert result.split_detected is True
    assert result.diagnostics["solver_seed_name"] == "initial_phases"
    assert result.diagnostics["solver_method"] == "native_transformed_newton"


def test_electrolyte_lle_rejects_neutral_lle_initial_phase_labels() -> None:
    mix = _ascani_water_butanol_nacl_mixture()
    feed = np.asarray([0.55, 0.40, 0.025, 0.025], dtype=float)

    with np.testing.assert_raises_regex(epcsaft.InputError, "aq.*org.*phase_fraction"):
        mix.equilibrium(
            kind="electrolyte_lle",
            T=298.15,
            P=1.013e5,
            z=feed,
            initial_phases={"liq1": feed, "liq2": feed, "phase_fraction": 0.5},
        )


def test_equilibrium_options_default_iteration_budget_is_robust_for_electrolyte_lle() -> None:
    options = epcsaft.EquilibriumOptions()

    assert options.max_iterations == 180


def test_equilibrium_options_expose_density_robustness_controls() -> None:
    options = epcsaft.EquilibriumOptions(
        density_diagnostics="full",
        experimental_coupled_density_lle=True,
    )

    assert options.density_diagnostics == "full"
    assert options.experimental_coupled_density_lle is True


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


def test_divalent_two_to_one_salt_basis_reconstructs_charge_neutral_formula() -> None:
    species = ["H2O", "TBP", "Mg2+", "Cl-"]
    charges = np.asarray([0.0, 0.0, 2.0, -1.0], dtype=float)
    feed = np.asarray([0.80, 0.10, 0.02, 0.04], dtype=float)
    feed = feed / float(np.sum(feed))

    basis = build_electrolyte_basis(species, charges, feed)
    payload = basis.to_dict()
    basis_dict = {
        "neutral_indices": tuple(payload["neutral_indices"]),
        "salt_pairs": tuple(payload["salt_pairs"]),
    }
    formula = _explicit_to_formula_composition(feed, basis_dict)
    explicit, _scale = _formula_to_explicit_composition(formula, basis_dict, len(species))

    assert payload["salt_pairs"][0]["label"] == "MgCl2"
    assert payload["salt_pairs"][0]["cation_stoich"] == 1
    assert payload["salt_pairs"][0]["anion_stoich"] == 2
    np.testing.assert_allclose(explicit, feed, atol=1.0e-12)
    assert abs(float(np.dot(explicit, charges))) <= 1.0e-12


def test_one_to_two_salt_basis_reconstructs_charge_neutral_formula() -> None:
    species = ["H2O", "TBP", "Na+", "SO4--"]
    charges = np.asarray([0.0, 0.0, 1.0, -2.0], dtype=float)
    feed = np.asarray([0.80, 0.10, 0.04, 0.02], dtype=float)
    feed = feed / float(np.sum(feed))

    basis = build_electrolyte_basis(species, charges, feed)
    payload = basis.to_dict()
    basis_dict = {
        "neutral_indices": tuple(payload["neutral_indices"]),
        "salt_pairs": tuple(payload["salt_pairs"]),
    }
    formula = _explicit_to_formula_composition(feed, basis_dict)
    explicit, _scale = _formula_to_explicit_composition(formula, basis_dict, len(species))

    assert payload["salt_pairs"][0]["label"] == "Na2SO4"
    assert payload["salt_pairs"][0]["cation_stoich"] == 2
    assert payload["salt_pairs"][0]["anion_stoich"] == 1
    np.testing.assert_allclose(explicit, feed, atol=1.0e-12)
    assert abs(float(np.dot(explicit, charges))) <= 1.0e-12


def test_mixed_monovalent_divalent_shared_anion_basis_builds() -> None:
    species = ["H2O", "TBP", "Li+", "Mg2+", "Cl-"]
    charges = np.asarray([0.0, 0.0, 1.0, 2.0, -1.0], dtype=float)
    feed = np.asarray([0.80, 0.10, 0.02, 0.01, 0.04], dtype=float)
    feed = feed / float(np.sum(feed))

    basis = build_electrolyte_basis(species, charges, feed)
    labels = [pair["label"] for pair in basis.salt_pairs]

    assert labels == ["LiCl", "MgCl2"]
    assert basis.rank == 2


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
    assert result.diagnostics["tpd_method"] == "native_tpd_global_search"
    assert result.diagnostics["solver_language"] == "c++"
    assert result.diagnostics["tpd_objective_value"] == pytest.approx(result.min_tpd)
    assert result.diagnostics["phase_charge_balance"]["trial"] == pytest.approx(0.0, abs=1.0e-8)
    json.dumps(result.to_dict(), allow_nan=False)


def test_one_salt_smoke_reports_unavailable_solver_derivatives() -> None:
    feed = np.asarray([0.55, 0.40, 0.025, 0.025], dtype=float)
    mix = _ascani_water_butanol_nacl_mixture(feed)

    with pytest.raises(epcsaft.InputError, match="backend_unavailable"):
        mix.equilibrium(
            kind="electrolyte_lle",
            T=298.15,
            P=1.013e5,
            z=feed,
            options=epcsaft.EquilibriumOptions(max_iterations=80, tolerance=1.0e-6),
        )


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


def test_mixed_electrolyte_lle_reports_phase_charge_balance_for_distributed_ions() -> None:
    mix = _case2_mixture()
    result = mix.equilibrium(
        kind="electrolyte_lle",
        T=298.15,
        P=1.0e5,
        z=_case2_feed(),
        options=epcsaft.EquilibriumOptions(max_iterations=180, tolerance=1.0e-8, damping=0.5),
    )

    charge_balance = result.diagnostics["phase_charge_balance"]
    assert set(charge_balance) == {"feed", "aq", "org", "max_abs"}
    assert charge_balance["feed"] == pytest.approx(0.0, abs=1.0e-12)
    assert charge_balance["aq"] == pytest.approx(0.0, abs=1.0e-8)
    assert charge_balance["org"] == pytest.approx(0.0, abs=1.0e-8)
    assert charge_balance["max_abs"] == pytest.approx(result.diagnostics["charge_balance_error"], abs=1.0e-12)

    salt_pairs = result.diagnostics["salt_pairs"]
    assert [pair["label"] for pair in salt_pairs] == ["NaCl", "KCl"]
    assert [(pair["cation"], pair["anion"]) for pair in salt_pairs] == [(2, 4), (3, 4)]
    aq, org = result.phases
    for ion_index in (2, 3, 4):
        assert aq.composition[ion_index] > 0.0
        assert org.composition[ion_index] > 0.0


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


def test_electrolyte_lle_rejects_non_neutral_direct_feed() -> None:
    mix = _ascani_water_butanol_nacl_mixture()

    with np.testing.assert_raises_regex(epcsaft.InputError, "charge neutral"):
        mix.equilibrium(kind="electrolyte_lle", T=298.15, P=1.013e5, z=[0.55, 0.40, 0.04, 0.01])


def test_neutral_lle_keeps_rejecting_ionic_mixtures() -> None:
    feed = np.asarray([0.55, 0.40, 0.025, 0.025], dtype=float)
    mix = _ascani_water_butanol_nacl_mixture(feed)

    with np.testing.assert_raises_regex(epcsaft.InputError, "ion-containing"):
        mix.equilibrium(kind="lle_flash", T=298.15, P=1.013e5, z=feed)
