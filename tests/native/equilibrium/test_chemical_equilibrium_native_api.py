from __future__ import annotations

import math

import numpy as np
import pytest

import epcsaft
from epcsaft import _core
from tests.api.reactive.test_reactive_speciation_options import _assert_reactive_speciation_route_pending
from tests.equilibrium.core.test_stability import _assert_stability_route_pending
from tests.helpers.numeric import assert_allclose


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

def _mea_like_mixture() -> epcsaft.ePCSAFTMixture:
    params = {
        "m": np.asarray([1.2047, 2.5, 1.0, 1.0, 2.2, 2.7, 1.0]),
        "s": np.asarray([2.7927, 3.3, 3.0, 3.2, 3.5, 3.6, 2.0]),
        "e": np.asarray([353.95, 260.0, 190.0, 230.0, 250.0, 245.0, 120.0]),
        "z": np.asarray([0.0, 0.0, 0.0, -1.0, 1.0, -1.0, 1.0]),
        "dielc": np.asarray([78.09, 35.0, 12.0, 20.0, 25.0, 22.0, 8.0]),
        "d_born": np.asarray([0.0, 0.0, 0.0, 3.8, 3.4, 4.0, 2.0]),
        "MW": np.asarray([18.01528e-3, 61.08e-3, 44.01e-3, 61.02e-3, 62.09e-3, 104.1e-3, 1.008e-3]),
    }
    return epcsaft.ePCSAFTMixture.from_params(
        params,
        species=["H2O", "MEA", "CO2", "HCO3-", "MEAH+", "MEACOO-", "H+"],
    )

def _methanol_cyclohexane_mixture(kij: float = 0.051) -> epcsaft.ePCSAFTMixture:
    params = {
        "MW": np.asarray([32.042e-3, 84.147e-3]),
        "m": np.asarray([1.5255, 2.5303]),
        "s": np.asarray([3.2300, 3.8499]),
        "e": np.asarray([188.90, 278.11]),
        "e_assoc": np.asarray([2899.5, 0.0]),
        "vol_a": np.asarray([0.035176, 0.0]),
        "assoc_scheme": ["2B", None],
        "k_ij": np.asarray([[0.0, kij], [kij, 0.0]]),
        "z": np.asarray([0.0, 0.0]),
        "dielc": np.asarray([33.05, 2.02]),
    }
    return epcsaft.ePCSAFTMixture.from_params(params, species=["Methanol", "Cyclohexane"])

def _log_k_from_state(
    mix: epcsaft.ePCSAFTMixture,
    T: float,
    P: float,
    x: np.ndarray,
    stoichiometry: dict[str, float],
) -> float:
    state = mix.state(T=T, P=P, x=x, phase="liq")
    gamma = state.activity_coefficient(species=mix.species)
    return float(
        sum(
            nu * math.log(max(x[mix.species.index(label)] * gamma[label], 1.0e-300))
            for label, nu in stoichiometry.items()
        )
    )

def _neutral_log_k_from_fugacity_activity(
    mix: epcsaft.ePCSAFTMixture,
    T: float,
    P: float,
    x: np.ndarray,
    stoichiometry: dict[str, float],
) -> float:
    state = mix.state(T=T, P=P, x=x, phase="liq")
    ln_phi = state.fugacity_coefficient(natural_log=True)
    ln_gamma = []
    for idx in range(mix.ncomp):
        x_ref = np.full(mix.ncomp, 1.0e-14, dtype=float)
        x_ref[idx] = 1.0 - 1.0e-14 * float(mix.ncomp - 1)
        ref = mix.state(T=T, P=P, x=x_ref, phase="liq")
        ln_phi_ref = ref.fugacity_coefficient(natural_log=True)
        ln_gamma.append(float(ln_phi[idx] - ln_phi_ref[idx]))
    return float(
        sum(
            nu * (math.log(max(x[mix.species.index(label)], 1.0e-300)) + ln_gamma[mix.species.index(label)])
            for label, nu in stoichiometry.items()
        )
    )

def test_native_chemical_equilibrium_entrypoint_is_exposed() -> None:
    assert hasattr(_core, "_solve_chemical_equilibrium_native")
    assert hasattr(_core, "_evaluate_chemical_equilibrium_residual_native")

def test_native_chemical_equilibrium_residual_evaluator_uses_analytic_jacobian_by_default() -> None:
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )
    request = {
        "T": 298.15,
        "P": 1.0e5,
        "initial_x": [0.5, 0.5],
        "balance_matrix": [1.0, 1.0],
        "balance_rows": 1,
        "total_vector": [1.0],
        "reaction_stoichiometry": [-1.0, 1.0],
        "reaction_rows": 1,
        "log_equilibrium_constants": [math.log(3.0)],
        "reaction_standard_states": [1],
        "options": {"tolerance": 1.0e-10},
    }

    payload = _core._evaluate_chemical_equilibrium_residual_native(mix._native, request)

    assert payload["variable_model"] == "log_species_amounts"
    assert payload["jacobian_backend"] == "analytic"
    removed_hessian_key = "hessian" + "_backend"
    assert removed_hessian_key not in payload
    diagnostics = payload["diagnostics"]
    assert diagnostics["derivative_backend"] == "analytic"
    removed_derivative_label_key = "derivative" + "_status"
    assert removed_derivative_label_key not in diagnostics
    assert diagnostics["derivative_capability_path"] == "chemical_equilibrium:ideal_mole_fraction:log_amounts"
    assert diagnostics["derivative_available"] is True
    removed_reason_key = "not" + "_available_reason"
    assert removed_reason_key not in diagnostics
    assert removed_hessian_key not in diagnostics
    assert "hessian_kind" not in diagnostics
    assert diagnostics["hessian_available"] is False
    assert diagnostics["exact_hessian_available"] is False
    assert diagnostics["hessian_callback_available"] is False
    assert diagnostics["hessian_includes_second_residual_derivatives"] is False
    residual = np.asarray(payload["residual"], dtype=float)
    gradient = np.asarray(payload["gradient"], dtype=float)
    jacobian = np.asarray(payload["jacobian_row_major"], dtype=float).reshape(payload["jacobian_shape"])
    lower = np.asarray(payload["lower_bounds"], dtype=float)
    variables = np.asarray(payload["variables"], dtype=float)
    upper = np.asarray(payload["upper_bounds"], dtype=float)
    assert residual.shape == (3,)
    assert gradient.shape == (2,)
    assert jacobian.shape == (3, 2)
    assert np.isfinite(payload["objective"])
    assert np.all(np.isfinite(residual))
    assert np.all(np.isfinite(gradient))
    assert np.all(np.isfinite(jacobian))
    assert np.all(np.isfinite(lower))
    assert np.all(np.isfinite(variables))
    assert np.all(np.isfinite(upper))
    assert np.all(lower < upper)
    assert np.all(variables >= lower)
    assert np.all(variables <= upper)
    assert_allclose(gradient, jacobian.T @ residual, rtol=1.0e-10, atol=1.0e-10)
    assert payload["objective"] == pytest.approx(0.5 * float(residual @ residual))
    assert len(payload["lower_bounds"]) == len(payload["variables"]) == len(payload["upper_bounds"])

def test_mixture_equilibrium_routes_chemical_equilibrium_to_native_ipopt_gate() -> None:
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="chemical_equilibrium",
            T=298.15,
            P=1.0e5,
            z=[0.5, 0.5],
            balances={"total": {"A": 1.0, "B": 1.0}},
            totals={"total": 1.0},
            reactions=[
                epcsaft.ReactionDefinition(
                    {"A": -1.0, "B": 1.0},
                    math.log(3.0),
                    standard_state="ideal_mole_fraction",
                )
            ],
            options=epcsaft.ReactiveSpeciationOptions(tolerance=1.0e-10),
        )

    _assert_reactive_speciation_route_pending(excinfo)


def test_reactive_stability_requires_native_ipopt_stability_route_after_speciation(monkeypatch) -> None:
    mix = _methanol_cyclohexane_mixture()
    target_x = np.asarray([0.45, 0.55], dtype=float)

    def successful_speciation(*_args, **_kwargs):
        return epcsaft.ReactiveSpeciationResult(
            success=True,
            message="converged",
            x={"Methanol": float(target_x[0]), "Cyclohexane": float(target_x[1])},
            activity_coefficients={},
            mass_balance_residuals={"total": 0.0},
            charge_residual=0.0,
            reaction_residuals=[0.0],
            named_reaction_residuals={"methanol_to_cyclohexane": 0.0},
            state_failure_count=0,
            diagnostics={"phase_equilibrium_handoff": {}},
        )

    monkeypatch.setattr(mix, "chemical_equilibrium", successful_speciation)

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="reactive_stability",
            T=298.15,
            P=1.013e5,
            z=[0.5, 0.5],
            balances={"total": {"Methanol": 1.0, "Cyclohexane": 1.0}},
            totals={"total": 1.0},
            reactions=[
                epcsaft.ReactionDefinition(
                    {"Methanol": -1.0, "Cyclohexane": 1.0},
                    0.0,
                    name="methanol_to_cyclohexane",
                )
            ],
            parent_phase="liq",
            trial_phases=("liq",),
            options=epcsaft.ReactiveSpeciationOptions(tolerance=1.0e-10),
        )

    _assert_stability_route_pending(excinfo)
