from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from epcsaft import _core, ePCSAFTMixture
from epcsaft.equilibrium_core.electrolyte_basis import (
    build_electrolyte_basis,
    explicit_to_formula_composition,
    formula_to_explicit_composition,
)
from tests.equilibrium.core.test_stability import _assert_stability_native_ipopt_gate
from tests.helpers.numeric import assert_allclose


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
    formula = explicit_to_formula_composition(feed, basis_dict)
    explicit, _scale = formula_to_explicit_composition(formula, basis_dict, len(species))

    assert payload["salt_pairs"][0]["label"] == "MgCl2"
    assert payload["salt_pairs"][0]["cation_stoich"] == 1
    assert payload["salt_pairs"][0]["anion_stoich"] == 2
    assert_allclose(explicit, feed, atol=1.0e-12)
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
    formula = explicit_to_formula_composition(feed, basis_dict)
    explicit, _scale = formula_to_explicit_composition(formula, basis_dict, len(species))

    assert payload["salt_pairs"][0]["label"] == "Na2SO4"
    assert payload["salt_pairs"][0]["cation_stoich"] == 2
    assert payload["salt_pairs"][0]["anion_stoich"] == 1
    assert_allclose(explicit, feed, atol=1.0e-12)
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

def test_electrolyte_stability_builds_native_route_request_before_ipopt_gate(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mix = _case2_mixture()
    calls: list[dict[str, object]] = []

    def fake_route(
        _native,
        temperature,
        pressure,
        feed_composition,
        max_iterations,
        tolerance,
        timeout_seconds,
        hessian_mode,
        iteration_history_limit,
        stability_tolerance,
        trial_initial_composition,
        continuation_state,
    ):
        calls.append(
            {
                "temperature": temperature,
                "pressure": pressure,
                "feed_composition": feed_composition,
                "max_iterations": max_iterations,
                "tolerance": tolerance,
                "timeout_seconds": timeout_seconds,
                "stability_tolerance": stability_tolerance,
            }
        )
        return {
            "backend": "ipopt",
            "compiled": False,
            "ran": False,
            "accepted": False,
            "status": "ipopt_dependency_required",
        }

    monkeypatch.setattr(_core, "_native_electrolyte_stability_tpd_route_result", fake_route)

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="electrolyte_stability",
            T=298.15,
            P=1.0e5,
            z=_case2_feed(),
            options=epcsaft.EquilibriumOptions(max_iterations=80, tolerance=1.0e-8),
        )

    _assert_stability_native_ipopt_gate(excinfo, route="electrolyte_stability")
    assert len(calls) == 1
    call = calls[0]
    assert call["temperature"] == pytest.approx(298.15)
    assert call["pressure"] == pytest.approx(1.0e5)
    assert call["feed_composition"] == pytest.approx(_case2_feed())
    assert call["max_iterations"] == 80
    assert call["tolerance"] == pytest.approx(1.0e-8)
    assert call["timeout_seconds"] == pytest.approx(0.0)
    assert call["stability_tolerance"] == pytest.approx(1.0e-8)


def test_electrolyte_stability_converts_accepted_native_route_payload(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mix = _case2_mixture()
    feed = _case2_feed()

    def fake_route(
        _native,
        temperature,
        pressure,
        feed_composition,
        max_iterations,
        tolerance,
        timeout_seconds,
        hessian_mode,
        iteration_history_limit,
        stability_tolerance,
        trial_initial_composition,
        continuation_state,
    ):
        assert temperature == pytest.approx(298.15)
        assert pressure == pytest.approx(1.0e5)
        assert feed_composition == pytest.approx(feed)
        assert max_iterations > 0
        assert tolerance > 0.0
        assert timeout_seconds == pytest.approx(0.0)
        assert stability_tolerance > 0.0
        return {
            "backend": "ipopt",
            "compiled": True,
            "ran": True,
            "solver_accepted": True,
            "accepted": True,
            "stable": True,
            "status": "accepted",
            "solver_status": "Solve_Succeeded",
            "application_status": "Solve_Succeeded",
            "parent_phase": "liq",
            "trial_phase": "liq",
            "seed_name": "canonical_charge_neutral_feed",
            "min_tpd": 4.0e-6,
            "objective": 4.0e-6,
            "trial_composition": feed.tolist(),
            "constraints": [0.0, 0.0],
            "derivative_backend": "cppad_implicit",
        }

    monkeypatch.setattr(_core, "_native_electrolyte_stability_tpd_route_result", fake_route)

    result = mix.equilibrium(
        kind="electrolyte_stability",
        T=298.15,
        P=1.0e5,
        z=feed,
    )

    assert isinstance(result, epcsaft.StabilityResult)
    assert result.backend == "native_equilibrium_nlp"
    assert result.problem_kind == "electrolyte_stability"
    assert result.stable is True
    assert result.min_tpd == pytest.approx(4.0e-6)
    assert result.parent_phase == "liq"
    assert result.trial_phase == "liq"
    assert result.trial_composition == pytest.approx(feed)
    assert len(result.trials) == 1
    assert result.trials[0].seed_name == "canonical_charge_neutral_feed"
    assert result.trials[0].diagnostics["constraints"] == pytest.approx([0.0, 0.0])
    assert result.diagnostics["route_count"] == 1
    assert result.diagnostics["derivative_backend"] == "cppad_implicit"
