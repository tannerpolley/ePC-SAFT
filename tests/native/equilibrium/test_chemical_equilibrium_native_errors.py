from __future__ import annotations

import math

import numpy as np
import pytest

import epcsaft
from epcsaft import _core
from epcsaft.epcsaft import create_struct


def _toy_mixture() -> epcsaft.ePCSAFTMixture:
    return epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )


def test_native_chemical_equilibrium_residual_evaluator_rejects_removed_backend() -> None:
    mix = _toy_mixture()
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
        "options": {"jacobian_backend": "finite" + "_difference", "finite" + "_difference_step": 1.0e-7},
    }

    with pytest.raises(_core.NativeValueError, match="jacobian_backend"):
        _core._evaluate_chemical_equilibrium_residual_native(mix._native, request)


@pytest.mark.parametrize(
    ("standard_state", "standard_state_code"),
    [("mole_fraction_activity", 0), ("concentration", 2)],
)
def test_native_chemical_equilibrium_residual_evaluator_uses_phase_state_cppad_derivatives(
    standard_state: str,
    standard_state_code: int,
) -> None:
    mix = _toy_mixture()
    request = {
        "T": 298.15,
        "P": 1.0e5,
        "initial_x": [0.4, 0.6],
        "balance_matrix": [1.0, 1.0],
        "balance_rows": 1,
        "total_vector": [1.0],
        "reaction_stoichiometry": [-1.0, 2.0],
        "reaction_rows": 1,
        "log_equilibrium_constants": [math.log(0.75)],
        "reaction_standard_states": [standard_state_code],
        "options": {"jacobian_backend": "auto", "activity_output": "always"},
    }

    payload = _core._evaluate_chemical_equilibrium_residual_native(mix._native, request)

    diagnostics = payload["diagnostics"]
    assert payload["jacobian_backend"] == "cppad_implicit"
    assert diagnostics["derivative_backend"] == "cppad_implicit"
    assert diagnostics["derivative_available"] is True
    assert diagnostics["activity_model"] == "eos_phase_state"
    assert diagnostics["activity_basis"] in {standard_state, "mole_fraction"}
    assert "phase_state_cppad_implicit" in diagnostics["derivative_capability_path"]
    assert diagnostics["activity_coefficients_evaluated"] is True

    x = np.asarray(payload["composition"], dtype=float)
    sensitivity = _core._native_phase_state_ln_fugacity_composition_sensitivity(
        request["T"],
        request["P"],
        x.tolist(),
        0,
        create_struct(mix.parameters),
    )
    assert sensitivity["supported"] is True
    ln_phi = np.asarray(sensitivity["ln_fugacity"], dtype=float)
    density = float(sensitivity["density"])
    stoich = np.asarray(request["reaction_stoichiometry"], dtype=float)
    if standard_state == "mole_fraction_activity":
        log_terms = np.log(x) + ln_phi
    else:
        log_terms = np.log(x * density)
    expected_residual = float(stoich @ log_terms - request["log_equilibrium_constants"][0])
    assert payload["reaction_residuals"][0] == pytest.approx(expected_residual, abs=1.0e-10)

    jacobian = np.asarray(payload["jacobian_row_major"], dtype=float).reshape(payload["jacobian_shape"])
    reaction_jacobian = jacobian[-1]
    logx_jacobian = np.asarray(
        [[(species == variable) - x[variable] for variable in range(x.size)] for species in range(x.size)],
        dtype=float,
    )
    if standard_state == "mole_fraction_activity":
        lnphi_jacobian = np.asarray(sensitivity["jacobian_row_major"], dtype=float).reshape((x.size, x.size))
        expected_species_jacobian = np.empty_like(logx_jacobian)
        for species in range(x.size):
            for variable in range(x.size):
                dlnphi = sum(
                    lnphi_jacobian[species, k] * x[k] * ((k == variable) - x[variable])
                    for k in range(x.size)
                )
                expected_species_jacobian[species, variable] = logx_jacobian[species, variable] + dlnphi
    else:
        drhodx = np.asarray(sensitivity["density_composition_derivative"], dtype=float)
        expected_species_jacobian = np.empty_like(logx_jacobian)
        for species in range(x.size):
            for variable in range(x.size):
                drho = sum(drhodx[k] * x[k] * ((k == variable) - x[variable]) for k in range(x.size))
                expected_species_jacobian[species, variable] = logx_jacobian[species, variable] + drho / density
    expected_reaction_jacobian = stoich @ expected_species_jacobian
    np.testing.assert_allclose(reaction_jacobian, expected_reaction_jacobian, rtol=1.0e-10, atol=1.0e-10)


def test_native_chemical_equilibrium_solve_accepts_ideal_cppad_derivative_request_when_compiled() -> None:
    mix = _toy_mixture()
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
        "options": {
            "solver_backend": "ipopt",
            "jacobian_backend": "cppad",
            "hessian_mode": "limited-memory",
            "iteration_history_limit": 2,
            "linear_solver": "mumps",
            "acceptable_tolerance": 9.0e-7,
            "constraint_violation_tolerance": 8.0e-8,
            "dual_infeasibility_tolerance": 7.0e-8,
            "complementarity_tolerance": 6.0e-8,
        },
    }

    if not _core._native_ipopt_smoke()["compiled"]:
        with pytest.raises(_core.NativeSolutionError, match=r"EPCSAFT_ENABLE_IPOPT=ON"):
            _core._solve_chemical_equilibrium_native(mix._native, request)
        return

    payload = _core._solve_chemical_equilibrium_native(mix._native, request)

    assert payload["diagnostics"]["requested_jacobian_backend"] == "cppad"
    assert payload["diagnostics"]["derivative_backend"] == "cppad"
    assert payload["diagnostics"]["implicit_sensitivity_backend"] == "cppad_implicit"
    assert payload["diagnostics"]["hessian_approximation"] == "limited-memory"
    assert payload["diagnostics"]["iteration_history_limit"] == 2
    assert payload["diagnostics"]["linear_solver_requested"] == "mumps"
    assert payload["diagnostics"]["acceptable_tolerance"] == pytest.approx(9.0e-7)
    assert payload["diagnostics"]["constraint_violation_tolerance"] == pytest.approx(8.0e-8)
    assert payload["diagnostics"]["dual_infeasibility_tolerance"] == pytest.approx(7.0e-8)
    assert payload["diagnostics"]["complementarity_tolerance"] == pytest.approx(6.0e-8)
    assert len(payload["diagnostics"]["iteration_history"]) <= 2
    assert len(payload["diagnostics"]["continuation_state"]["variables"]) == 2


def test_native_chemical_equilibrium_solve_routes_nonideal_speciation_to_ipopt() -> None:
    mix = _toy_mixture()
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
        "reaction_standard_states": [0],
        "options": {
            "solver_backend": "ipopt",
            "jacobian_backend": "cppad",
            "linear_solver": "mumps",
            "acceptable_tolerance": 9.0e-7,
            "constraint_violation_tolerance": 8.0e-8,
            "dual_infeasibility_tolerance": 7.0e-8,
            "complementarity_tolerance": 6.0e-8,
        },
    }

    if not _core._native_ipopt_smoke()["compiled"]:
        with pytest.raises(_core.NativeSolutionError, match=r"EPCSAFT_ENABLE_IPOPT=ON"):
            _core._solve_chemical_equilibrium_native(mix._native, request)
        return

    payload = _core._solve_chemical_equilibrium_native(mix._native, request)

    assert payload["success"] is True
    assert payload["composition"][1] / payload["composition"][0] == pytest.approx(3.0, rel=1.0e-7)
    diagnostics = payload["diagnostics"]
    assert diagnostics["problem_class"] == "homogeneous_nonideal_residual_speciation"
    assert diagnostics["derivative_backend"] == "cppad_implicit"
    assert diagnostics["jacobian_backend"] == "cppad_implicit"
    assert diagnostics["implicit_sensitivity_backend"] == "cppad_implicit"
    assert diagnostics["selected_solver_backend"] == "native_ipopt"
    assert diagnostics["linear_solver_requested"] == "mumps"
    assert diagnostics["acceptable_tolerance"] == pytest.approx(9.0e-7)
    assert diagnostics["constraint_violation_tolerance"] == pytest.approx(8.0e-8)
    assert diagnostics["dual_infeasibility_tolerance"] == pytest.approx(7.0e-8)
    assert diagnostics["complementarity_tolerance"] == pytest.approx(6.0e-8)


@pytest.mark.parametrize("standard_state_code", [0, 1])
def test_native_chemical_equilibrium_solve_rejects_exact_hessian_without_provider(
    standard_state_code: int,
) -> None:
    mix = _toy_mixture()
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
        "reaction_standard_states": [standard_state_code],
        "options": {"solver_backend": "ipopt", "hessian_mode": "exact"},
    }

    if not _core._native_ipopt_smoke()["compiled"]:
        with pytest.raises(_core.NativeSolutionError, match=r"EPCSAFT_ENABLE_IPOPT=ON"):
            _core._solve_chemical_equilibrium_native(mix._native, request)
        return

    with pytest.raises(_core.NativeValueError, match="exact Hessian mode requires"):
        _core._solve_chemical_equilibrium_native(mix._native, request)


def test_native_chemical_equilibrium_solve_accepts_mixed_nonideal_standard_states() -> None:
    mix = _toy_mixture()
    request = {
        "T": 298.15,
        "P": 1.0e5,
        "initial_x": [0.5, 0.5],
        "balance_matrix": [1.0, 1.0],
        "balance_rows": 1,
        "total_vector": [1.0],
        "reaction_stoichiometry": [-1.0, 1.0, 1.0, -1.0],
        "reaction_rows": 2,
        "log_equilibrium_constants": [math.log(3.0), math.log(1.0 / 3.0)],
        "reaction_standard_states": [0, 2],
        "options": {"solver_backend": "ipopt", "jacobian_backend": "cppad"},
    }

    if not _core._native_ipopt_smoke()["compiled"]:
        with pytest.raises(_core.NativeSolutionError, match=r"EPCSAFT_ENABLE_IPOPT=ON"):
            _core._solve_chemical_equilibrium_native(mix._native, request)
        return

    payload = _core._solve_chemical_equilibrium_native(mix._native, request)

    assert payload["success"] is True
    diagnostics = payload["diagnostics"]
    assert diagnostics["problem_class"] == "homogeneous_nonideal_residual_speciation"
    assert diagnostics["activity_basis"] == "mixed_standard_state"


def test_mixture_equilibrium_rejects_non_native_chemical_equilibrium_backend() -> None:
    mix = _toy_mixture()

    with pytest.raises(epcsaft.InputError, match="native"):
        mix.equilibrium(
            kind="chemical_equilibrium",
            T=298.15,
            P=1.0e5,
            z=[0.5, 0.5],
            balances={"total": {"A": 1.0, "B": 1.0}},
            totals={"total": 1.0},
            reactions=[epcsaft.ReactionDefinition({"A": -1.0, "B": 1.0}, math.log(3.0))],
            backend="python",
        )
