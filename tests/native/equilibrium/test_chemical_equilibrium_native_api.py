from __future__ import annotations

import math

import numpy as np
import pytest

import epcsaft
from epcsaft import _core
from tests.equilibrium.core.test_stability import _assert_stability_native_ipopt_gate
from tests.helpers.numeric import assert_allclose


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
    diagnostics = payload["diagnostics"]
    assert diagnostics["derivative_backend"] == "analytic"
    assert diagnostics["derivative_capability_path"] == "chemical_equilibrium:ideal_mole_fraction:log_amounts"
    assert diagnostics["derivative_available"] is True
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


def test_native_chemical_equilibrium_residual_evaluator_uses_cppad_when_requested() -> None:
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
        "options": {"tolerance": 1.0e-10, "jacobian_backend": "cppad"},
    }

    cppad_payload = _core._evaluate_chemical_equilibrium_residual_native(mix._native, request)
    analytic_request = dict(request)
    analytic_request["options"] = {"tolerance": 1.0e-10, "jacobian_backend": "analytic"}
    analytic_payload = _core._evaluate_chemical_equilibrium_residual_native(mix._native, analytic_request)

    diagnostics = cppad_payload["diagnostics"]
    assert cppad_payload["jacobian_backend"] == "cppad"
    assert diagnostics["derivative_backend"] == "cppad"
    assert diagnostics["derivative_capability_path"] == "chemical_equilibrium:ideal_mole_fraction:cppad_log_amounts"
    assert diagnostics["derivative_available"] is True
    assert_allclose(
        np.asarray(cppad_payload["jacobian_row_major"], dtype=float),
        np.asarray(analytic_payload["jacobian_row_major"], dtype=float),
        rtol=1.0e-10,
        atol=1.0e-10,
    )

def test_mixture_equilibrium_auto_routes_ideal_chemical_equilibrium_to_native_ipopt_when_compiled() -> None:
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )
    kwargs = {
        "kind": "chemical_equilibrium",
        "T": 298.15,
        "P": 1.0e5,
        "z": [0.5, 0.5],
        "balances": {"total": {"A": 1.0, "B": 1.0}},
        "totals": {"total": 1.0},
        "reactions": [
            epcsaft.ReactionDefinition(
                {"A": -1.0, "B": 1.0},
                math.log(3.0),
                standard_state="ideal_mole_fraction",
            )
        ],
        "options": epcsaft.ReactiveSpeciationOptions(tolerance=1.0e-10),
    }

    if not _core._native_ipopt_smoke()["compiled"]:
        with pytest.raises(epcsaft.SolutionError, match=r"EPCSAFT_ENABLE_IPOPT=ON"):
            mix.equilibrium(**kwargs)
        return

    result = mix.equilibrium(**kwargs)

    assert result.success is True
    assert result.diagnostics["requested_solver_backend"] == "auto"
    assert result.diagnostics["selected_solver_backend"] == "native_ipopt"


def test_reactive_stability_requires_native_ipopt_stability_route_before_speciation(monkeypatch) -> None:
    mix = _methanol_cyclohexane_mixture()

    def successful_speciation(*_args, **_kwargs):
        pytest.fail("reactive_stability must not run a Python speciation handoff before the native stability gate")

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

    _assert_stability_native_ipopt_gate(excinfo, route="reactive_stability")
