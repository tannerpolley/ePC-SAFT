from __future__ import annotations

import math

import numpy as np
import pytest

import epcsaft
from tests.api.reactive.test_reactive_speciation_options import _assert_reactive_speciation_route_pending


def _salt_speciation_mixture() -> epcsaft.ePCSAFTMixture:
    params = {
        "m": np.asarray([1.2047, 1.0, 1.0, 1.0]),
        "s": np.asarray([2.7927, 3.0, 2.8232, 2.7560]),
        "e": np.asarray([353.95, 200.0, 230.0, 170.0]),
        "z": np.asarray([0.0, 0.0, 1.0, -1.0]),
        "dielc": np.asarray([78.09, 8.0, 8.0, 8.0]),
        "d_born": np.asarray([0.0, 0.0, 3.445, 4.1]),
        "MW": np.asarray([18.01528e-3, 58.44e-3, 22.989e-3, 35.45e-3]),
    }
    return epcsaft.ePCSAFTMixture.from_params(params, species=["H2O", "NaCl", "Na+", "Cl-"])


def _salt_speciation_request(standard_state: str = "mole_fraction_activity") -> dict[str, object]:
    species = ["H2O", "NaCl", "Na+", "Cl-"]
    mix = _salt_speciation_mixture()
    initial_x = np.asarray([0.998, 0.001, 0.0005, 0.0005], dtype=float)
    state = mix.state(T=298.15, P=1.0e5, x=initial_x, phase="liq")
    if standard_state == "concentration":
        density = state.molar_density()
        log_k = math.log(density * initial_x[2]) + math.log(density * initial_x[3])
        log_k -= math.log(density * initial_x[1])
    else:
        gamma = state.activity_coefficient(species=species)
        log_k = math.log(initial_x[2] * gamma["Na+"]) + math.log(initial_x[3] * gamma["Cl-"])
        log_k -= math.log(initial_x[1] * gamma["NaCl"])
    return {
        "species": species,
        "mixture_factory": lambda x, T, P: mix,
        "T": 298.15,
        "P": 1.0e5,
        "balances": {
            "water_total": {"H2O": 1.0},
            "sodium_total": {"NaCl": 1.0, "Na+": 1.0},
            "chloride_total": {"NaCl": 1.0, "Cl-": 1.0},
        },
        "totals": {"water_total": 0.998, "sodium_total": 0.0015, "chloride_total": 0.0015},
        "reactions": [
            epcsaft.ReactionDefinition(
                stoichiometry={"NaCl": -1.0, "Na+": 1.0, "Cl-": 1.0},
                log_equilibrium_constant=log_k,
                name="salt_dissociation",
                standard_state=standard_state,
            )
        ],
        "initial_x": initial_x,
    }


def test_solve_reactive_speciation_activity_coupled_state_requires_native_ipopt_route() -> None:
    with pytest.raises(epcsaft.InputError) as excinfo:
        epcsaft.solve_reactive_speciation(
            **_salt_speciation_request(),
            options=epcsaft.ReactiveSpeciationOptions(max_iterations=8, tolerance=1.0e-8),
        )

    _assert_reactive_speciation_route_pending(excinfo)


def test_solve_reactive_speciation_concentration_standard_state_requires_native_ipopt_route() -> None:
    with pytest.raises(epcsaft.InputError) as excinfo:
        epcsaft.solve_reactive_speciation(
            **_salt_speciation_request("concentration"),
            options=epcsaft.ReactiveSpeciationOptions(max_iterations=50, tolerance=1.0e-8),
        )

    _assert_reactive_speciation_route_pending(excinfo)


def test_concentration_standard_state_can_validate_activity_output_before_route_gate() -> None:
    with pytest.raises(epcsaft.InputError) as excinfo:
        epcsaft.solve_reactive_speciation(
            **_salt_speciation_request("concentration"),
            options=epcsaft.ReactiveSpeciationOptions(activity_output="never"),
        )

    _assert_reactive_speciation_route_pending(excinfo)


def test_reactive_speciation_sweep_returns_failed_results_for_pending_auto_route() -> None:
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )

    results = epcsaft.solve_reactive_speciation_sweep(
        species=["A", "B"],
        mixture_factory=lambda x, T, P: mix,
        points=[
            {"T": 298.15, "P": 1.0e5, "totals": {"total": 1.0}, "initial_x": [0.5, 0.5]},
            {"T": 298.15, "P": 1.0e5, "totals": {"total": 1.0}, "initial_x": [0.9, 0.1]},
        ],
        balances={"total": {"A": 1.0, "B": 1.0}},
        reactions=[
            epcsaft.ReactionDefinition(
                {"A": -1.0, "B": 1.0},
                log_equilibrium_constant=math.log(3.0),
                standard_state="ideal_mole_fraction",
            )
        ],
        options=epcsaft.ReactiveSpeciationOptions(error_mode="result"),
        continuation="auto",
    )

    assert len(results) == 2
    assert all(isinstance(result, epcsaft.ReactiveSpeciationResult) for result in results)
    assert all(result.success is False for result in results)
    assert all("native Ipopt homogeneous reactive-speciation NLP route" in result.message for result in results)


def test_reactive_speciation_sweep_preserves_input_validation_failure_shape() -> None:
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )

    results = epcsaft.solve_reactive_speciation_sweep(
        species=["A", "B"],
        mixture_factory=lambda x, T, P: mix,
        points=[
            {"T": 298.15, "P": 1.0e5, "totals": {"missing": 1.0}, "initial_x": [0.5, 0.5]},
        ],
        balances={"total": {"A": 1.0, "B": 1.0}},
        reactions=[
            epcsaft.ReactionDefinition(
                {"A": -1.0, "B": 1.0},
                log_equilibrium_constant=math.log(3.0),
                standard_state="ideal_mole_fraction",
            )
        ],
        options=epcsaft.ReactiveSpeciationOptions(error_mode="result"),
    )

    assert len(results) == 1
    assert results[0].success is False
    assert "Missing total" in results[0].message
    assert results[0].diagnostics["structured_failure"] is True
