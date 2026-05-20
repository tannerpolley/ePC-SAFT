from __future__ import annotations

import math
from dataclasses import fields

import numpy as np
import pytest

import epcsaft


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


def _native_ipopt_compiled() -> bool:
    from epcsaft import _core

    return bool(_core._native_ipopt_smoke()["compiled"])


def _assert_reactive_speciation_native_ipopt_dependency_required(
    excinfo: pytest.ExceptionInfo[epcsaft.SolutionError],
) -> None:
    message = str(excinfo.value)
    assert "EPCSAFT_ENABLE_IPOPT=ON" in message


@pytest.mark.parametrize("standard_state", ["ideal_mole_fraction", "concentration", "mole_fraction_activity"])
def test_reaction_definition_accepts_supported_standard_states(standard_state: str) -> None:
    reaction = epcsaft.ReactionDefinition(
        stoichiometry={"NaCl": -1.0, "Na+": 1.0, "Cl-": 1.0},
        log_equilibrium_constant=0.0,
        standard_state=standard_state,
    )

    assert reaction.standard_state == standard_state

def test_reaction_definition_rejects_unknown_standard_state() -> None:
    with pytest.raises(epcsaft.InputError, match=r"ReactionDefinition\.standard_state"):
        epcsaft.ReactionDefinition(
            stoichiometry={"NaCl": -1.0, "Na+": 1.0, "Cl-": 1.0},
            log_equilibrium_constant=0.0,
            standard_state="unknown_basis",
        )

def test_reactive_speciation_options_public_surface_is_current_fields() -> None:
    assert {field.name for field in fields(epcsaft.ReactiveSpeciationOptions)} == {
        "max_iterations",
        "tolerance",
        "min_mole_fraction",
        "jacobian_backend",
        "solver_backend",
        "hessian_mode",
        "ipopt_iteration_history_limit",
        "ipopt_linear_solver",
        "ipopt_acceptable_tolerance",
        "ipopt_constraint_violation_tolerance",
        "ipopt_dual_infeasibility_tolerance",
        "ipopt_complementarity_tolerance",
        "continuation_state",
        "phase",
        "error_mode",
        "activity_output",
        "mass_tolerance",
        "charge_tolerance",
        "reaction_tolerance",
    }

@pytest.mark.parametrize(
    ("options", "message"),
    [
        (epcsaft.ReactiveSpeciationOptions(solver_backend="python_ipopt"), "solver_backend"),
        (epcsaft.ReactiveSpeciationOptions(jacobian_backend="autodiff"), "jacobian_backend"),
        (epcsaft.ReactiveSpeciationOptions(hessian_mode="unsupported-mode"), "hessian_mode"),
        (epcsaft.ReactiveSpeciationOptions(ipopt_iteration_history_limit=-1), "ipopt_iteration_history_limit"),
        (epcsaft.ReactiveSpeciationOptions(ipopt_iteration_history_limit=True), "ipopt_iteration_history_limit"),
        (epcsaft.ReactiveSpeciationOptions(ipopt_linear_solver=""), "ipopt_linear_solver"),
        (epcsaft.ReactiveSpeciationOptions(ipopt_acceptable_tolerance=0.0), "ipopt_acceptable_tolerance"),
        (
            epcsaft.ReactiveSpeciationOptions(ipopt_constraint_violation_tolerance=float("nan")),
            "ipopt_constraint_violation_tolerance",
        ),
        (
            epcsaft.ReactiveSpeciationOptions(ipopt_dual_infeasibility_tolerance=-1.0),
            "ipopt_dual_infeasibility_tolerance",
        ),
        (
            epcsaft.ReactiveSpeciationOptions(ipopt_complementarity_tolerance=True),
            "ipopt_complementarity_tolerance",
        ),
        (epcsaft.ReactiveSpeciationOptions(continuation_state=1), "continuation_state"),
    ],
)
def test_reactive_speciation_rejects_invalid_optimizer_options(options, message) -> None:
    species = ["H2O", "NaCl", "Na+", "Cl-"]
    mix = _salt_speciation_mixture()

    with pytest.raises(epcsaft.InputError, match=message):
        epcsaft.solve_reactive_speciation(
            species=species,
            mixture_factory=lambda x, T, P: mix,
            T=298.15,
            P=1.0e5,
            balances={
                "water_total": {"H2O": 1.0},
                "sodium_total": {"NaCl": 1.0, "Na+": 1.0},
                "chloride_total": {"NaCl": 1.0, "Cl-": 1.0},
            },
            totals={"water_total": 0.998, "sodium_total": 0.0015, "chloride_total": 0.0015},
            reactions=[
                epcsaft.ReactionDefinition(
                    stoichiometry={"NaCl": -1.0, "Na+": 1.0, "Cl-": 1.0},
                    log_equilibrium_constant=0.0,
                )
            ],
            initial_x=[0.998, 0.001, 0.0005, 0.0005],
            options=options,
        )


def test_reactive_speciation_builds_native_request_with_ipopt_tranche_options(monkeypatch) -> None:
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )
    recorded: dict[str, object] = {}

    def fake_native(_native, request):
        recorded["request"] = request
        return {
            "success": True,
            "message": "converged",
            "composition": [0.25, 0.75],
            "activity_coefficients": [1.0, 1.0],
            "mass_balance_residuals": [0.0],
            "charge_residual": 0.0,
            "reaction_residuals": [0.0],
            "diagnostics": {
                "derivative_backend": "analytic",
                "selected_solver_backend": "native_ipopt",
                "solver_selection_reason": "explicit_request",
                "hessian_approximation": "limited-memory",
                "hessian_backend": "limited-memory",
                "iteration_history_limit": 4,
                "linear_solver_requested": "mumps",
                "linear_solver_selected": "mumps",
                "acceptable_tolerance": 9.0e-7,
                "constraint_violation_tolerance": 8.0e-8,
                "dual_infeasibility_tolerance": 7.0e-8,
                "complementarity_tolerance": 6.0e-8,
                "iteration_history": [],
                "continuation_state": {
                    "variables": [0.25, 0.75],
                    "bound_lower_multipliers": [0.0, 0.0],
                    "bound_upper_multipliers": [0.0, 0.0],
                    "constraint_multipliers": [0.0],
                },
            },
        }

    monkeypatch.setattr(epcsaft._core, "_solve_chemical_equilibrium_native", fake_native)

    result = epcsaft.solve_reactive_speciation(
        species=["A", "B"],
        mixture_factory=lambda x, T, P: mix,
        T=298.15,
        P=1.0e5,
        balances={"total": {"A": 1.0, "B": 1.0}},
        totals={"total": 1.0},
        reactions=[
            epcsaft.ReactionDefinition(
                {"A": -1.0, "B": 1.0},
                log_equilibrium_constant=math.log(3.0),
                standard_state="ideal_mole_fraction",
            )
        ],
        initial_x=[0.5, 0.5],
        options=epcsaft.ReactiveSpeciationOptions(
            solver_backend="ipopt",
            hessian_mode="exact",
            ipopt_iteration_history_limit=4,
            ipopt_linear_solver="mumps",
            ipopt_acceptable_tolerance=9.0e-7,
            ipopt_constraint_violation_tolerance=8.0e-8,
            ipopt_dual_infeasibility_tolerance=7.0e-8,
            ipopt_complementarity_tolerance=6.0e-8,
            continuation_state={
                "variables": [0.5, 0.5],
                "bound_lower_multipliers": [0.0, 0.0],
                "bound_upper_multipliers": [0.0, 0.0],
                "constraint_multipliers": [0.0],
                "route_kind": "reactive_speciation",
                "species_order": ["A", "B"],
                "fixed_specs": {"fixed": ["T", "P", "totals"], "phase": "liq"},
            },
        ),
    )

    request = recorded["request"]
    assert request["options"]["hessian_mode"] == "exact"
    assert request["options"]["iteration_history_limit"] == 4
    assert request["options"]["linear_solver"] == "mumps"
    assert request["options"]["acceptable_tolerance"] == pytest.approx(9.0e-7)
    assert request["options"]["constraint_violation_tolerance"] == pytest.approx(8.0e-8)
    assert request["options"]["dual_infeasibility_tolerance"] == pytest.approx(7.0e-8)
    assert request["options"]["complementarity_tolerance"] == pytest.approx(6.0e-8)
    assert request["options"]["continuation_state"]["variables"] == pytest.approx([0.5, 0.5])
    assert result.success is True
    assert result.diagnostics["hessian_approximation"] == "limited-memory"
    assert result.diagnostics["linear_solver_requested"] == "mumps"
    assert result.diagnostics["linear_solver_selected"] == "mumps"
    assert result.diagnostics["acceptable_tolerance"] == pytest.approx(9.0e-7)
    assert result.diagnostics["constraint_violation_tolerance"] == pytest.approx(8.0e-8)
    assert result.diagnostics["dual_infeasibility_tolerance"] == pytest.approx(7.0e-8)
    assert result.diagnostics["complementarity_tolerance"] == pytest.approx(6.0e-8)
    assert result.diagnostics["continuation_state"]["route_kind"] == "reactive_speciation"
    assert result.diagnostics["continuation_state"]["species_order"] == ["A", "B"]


def test_reactive_speciation_rejects_incompatible_continuation_state_species_order() -> None:
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )

    with pytest.raises(epcsaft.InputError, match="species_order"):
        epcsaft.solve_reactive_speciation(
            species=["A", "B"],
            mixture_factory=lambda x, T, P: mix,
            T=298.15,
            P=1.0e5,
            balances={"total": {"A": 1.0, "B": 1.0}},
            totals={"total": 1.0},
            reactions=[
                epcsaft.ReactionDefinition(
                    {"A": -1.0, "B": 1.0},
                    log_equilibrium_constant=math.log(3.0),
                    standard_state="ideal_mole_fraction",
                )
            ],
            initial_x=[0.5, 0.5],
            options=epcsaft.ReactiveSpeciationOptions(
                solver_backend="ipopt",
                continuation_state={
                    "variables": [0.5, 0.5],
                    "species_order": ["B", "A"],
                },
            ),
        )

def test_reactive_speciation_requested_ipopt_routes_ideal_speciation_when_compiled() -> None:
    from epcsaft import _core

    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )
    kwargs = {
        "species": ["A", "B"],
        "mixture_factory": lambda x, T, P: mix,
        "T": 298.15,
        "P": 1.0e5,
        "balances": {"total": {"A": 1.0, "B": 1.0}},
        "totals": {"total": 1.0},
        "reactions": [
            epcsaft.ReactionDefinition(
                {"A": -1.0, "B": 1.0},
                log_equilibrium_constant=math.log(3.0),
                standard_state="ideal_mole_fraction",
            )
        ],
        "initial_x": [0.5, 0.5],
        "options": epcsaft.ReactiveSpeciationOptions(solver_backend="ipopt", tolerance=1.0e-9),
    }

    if not _core._native_ipopt_smoke()["compiled"]:
        with pytest.raises(epcsaft.SolutionError, match=r"EPCSAFT_ENABLE_IPOPT=ON"):
            epcsaft.solve_reactive_speciation(**kwargs)
        return

    result = epcsaft.solve_reactive_speciation(**kwargs)

    assert result.success is True
    assert result.x["B"] / result.x["A"] == pytest.approx(3.0, rel=1.0e-7)
    assert result.diagnostics["selected_solver_backend"] == "native_ipopt"
    assert result.diagnostics["problem_class"] == "homogeneous_ideal_gibbs_speciation"
    assert result.diagnostics["jacobian_backend"] == "analytic"

def test_reactive_speciation_requested_ipopt_handles_charged_ideal_constraint_when_compiled() -> None:
    from epcsaft import _core

    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0, 1.0]),
            "s": np.asarray([3.0, 3.0, 3.0]),
            "e": np.asarray([200.0, 200.0, 200.0]),
            "z": np.asarray([1.0, -1.0, 0.0]),
            "dielc": np.asarray([78.0, 78.0, 78.0]),
            "d_born": np.asarray([3.0, 3.0, 0.0]),
            "MW": np.asarray([20.0e-3, 20.0e-3, 40.0e-3]),
        },
        species=["C+", "A-", "N"],
    )
    kwargs = {
        "species": ["C+", "A-", "N"],
        "mixture_factory": lambda x, T, P: mix,
        "T": 298.15,
        "P": 1.0e5,
        "balances": {"formula_units": {"C+": 0.5, "A-": 0.5, "N": 1.0}},
        "totals": {"formula_units": 0.75},
        "reactions": [
            epcsaft.ReactionDefinition(
                {"C+": -1.0, "A-": -1.0, "N": 1.0},
                log_equilibrium_constant=math.log(8.0),
                standard_state="ideal_mole_fraction",
            )
        ],
        "initial_x": [0.3, 0.3, 0.4],
        "options": epcsaft.ReactiveSpeciationOptions(solver_backend="ipopt", tolerance=1.0e-9),
    }

    if not _core._native_ipopt_smoke()["compiled"]:
        with pytest.raises(epcsaft.SolutionError, match=r"EPCSAFT_ENABLE_IPOPT=ON"):
            epcsaft.solve_reactive_speciation(**kwargs)
        return

    result = epcsaft.solve_reactive_speciation(**kwargs)

    assert result.success is True
    assert result.x["C+"] == pytest.approx(0.25, rel=1.0e-7)
    assert result.x["A-"] == pytest.approx(0.25, rel=1.0e-7)
    assert result.x["N"] == pytest.approx(0.5, rel=1.0e-7)
    assert result.charge_residual == pytest.approx(0.0, abs=1.0e-10)
    assert result.diagnostics["charge_constraint_in_nlp"] is True

def test_reactive_speciation_auto_routes_ideal_speciation_to_native_ipopt_when_compiled() -> None:
    from epcsaft import _core

    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )
    kwargs = {
        "species": ["A", "B"],
        "mixture_factory": lambda x, T, P: mix,
        "T": 298.15,
        "P": 1.0e5,
        "balances": {"total": {"A": 1.0, "B": 1.0}},
        "totals": {"total": 1.0},
        "reactions": [
            epcsaft.ReactionDefinition(
                {"A": -1.0, "B": 1.0},
                log_equilibrium_constant=math.log(3.0),
                standard_state="ideal_mole_fraction",
            )
        ],
        "initial_x": [0.5, 0.5],
    }

    if not _core._native_ipopt_smoke()["compiled"]:
        with pytest.raises(epcsaft.SolutionError, match=r"EPCSAFT_ENABLE_IPOPT=ON"):
            epcsaft.solve_reactive_speciation(**kwargs)
        return

    result = epcsaft.solve_reactive_speciation(**kwargs)

    assert result.success is True
    assert result.diagnostics["requested_solver_backend"] == "auto"
    assert result.diagnostics["selected_solver_backend"] == "native_ipopt"
    assert result.diagnostics["solver_selection_reason"] == "auto_selected_native_ipopt"
    assert result.diagnostics["hessian_approximation"] == "exact"
    assert result.diagnostics["hessian_backend"] == "analytic"
    assert result.diagnostics["exact_hessian_available"] is True
    assert result.diagnostics["eval_h_calls"] > 0


def test_reactive_speciation_requested_ipopt_routes_nonideal_speciation_when_compiled() -> None:
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )
    kwargs = {
        "species": ["A", "B"],
        "mixture_factory": lambda x, T, P: mix,
        "T": 298.15,
        "P": 1.0e5,
        "balances": {"total": {"A": 1.0, "B": 1.0}},
        "totals": {"total": 1.0},
        "reactions": [
            epcsaft.ReactionDefinition(
                {"A": -1.0, "B": 1.0},
                log_equilibrium_constant=math.log(3.0),
                standard_state="mole_fraction_activity",
            )
        ],
        "initial_x": [0.5, 0.5],
        "options": epcsaft.ReactiveSpeciationOptions(solver_backend="ipopt", tolerance=1.0e-9),
    }

    if not _native_ipopt_compiled():
        with pytest.raises(epcsaft.SolutionError) as excinfo:
            epcsaft.solve_reactive_speciation(**kwargs)
        _assert_reactive_speciation_native_ipopt_dependency_required(excinfo)
        return

    result = epcsaft.solve_reactive_speciation(**kwargs)

    assert result.success is True
    assert result.x["B"] / result.x["A"] == pytest.approx(3.0, rel=1.0e-7)
    assert result.diagnostics["problem_class"] == "homogeneous_nonideal_residual_speciation"
    assert result.diagnostics["derivative_backend"] == "cppad_explicit_density"
    assert result.diagnostics["density_backend"] == "explicit_log_density_pressure_constraint"
    assert result.diagnostics["implicit_sensitivity_backend"] == "cppad_explicit_density_implicit"
    assert result.diagnostics["hessian_approximation"] == "exact"
    assert result.diagnostics["exact_hessian_available"] is True
    assert result.diagnostics["hessian_backend"] != "limited-memory"
    assert result.diagnostics["eval_h_calls"] > 0
