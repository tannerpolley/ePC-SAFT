from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from epcsaft import _core


def _toy_reactive_phase_case() -> tuple[
    epcsaft.ePCSAFTMixture,
    np.ndarray,
    epcsaft.ReactionDefinition,
]:
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "MW": np.asarray([32.042e-3, 84.147e-3]),
            "m": np.asarray([1.5255, 2.5303]),
            "s": np.asarray([3.2300, 3.8499]),
            "e": np.asarray([188.90, 278.11]),
            "e_assoc": np.asarray([2899.5, 0.0]),
            "vol_a": np.asarray([0.035176, 0.0]),
            "assoc_scheme": ["2B", None],
            "k_ij": np.asarray([[0.0, 0.051], [0.051, 0.0]]),
            "z": np.asarray([0.0, 0.0]),
            "dielc": np.asarray([33.05, 2.02]),
        },
        species=["Methanol", "Cyclohexane"],
    )
    beta2 = 0.48813098468607985
    liq1 = np.asarray([0.11757838279937723, 0.8824216172006228])
    liq2 = np.asarray([0.7985874308392054, 0.20141256916079467])
    feed = (1.0 - beta2) * liq1 + beta2 * liq2
    reaction = epcsaft.ReactionDefinition.from_literature_constant(
        {"Methanol": -1.0, "Cyclohexane": 1.0},
        log_equilibrium_constant=-0.079259405371,
        name="methanol_to_cyclohexane",
        standard_state="mole_fraction_activity",
        source="public route smoke fixture",
    )
    return mix, feed, reaction


def test_reactive_phase_equilibrium_problem_requires_native_ipopt_route(monkeypatch) -> None:
    mix, feed, reaction = _toy_reactive_phase_case()
    captured: dict[str, object] = {}

    def fail_if_staged(*_args, **_kwargs):
        raise AssertionError("ReactivePhaseEquilibriumProblem must not call the staged helper")

    def native_gate(*args, **kwargs):
        captured["args"] = args
        captured["kwargs"] = kwargs
        return {
            "backend": "ipopt",
            "compiled": False,
            "adapter_available": False,
            "problem_name": "reactive_two_phase_eos",
            "derivative_backend": "analytic_cppad",
            "accepted": False,
            "status": "ipopt_dependency_required",
            "postsolve": {},
        }

    monkeypatch.setattr(mix, "reactive_staged_equilibrium", fail_if_staged)
    monkeypatch.setattr(_core, "_native_reactive_lle_eos_route_result", native_gate)
    problem = epcsaft.ReactivePhaseEquilibriumProblem(
        T=298.15,
        P=1.013e5,
        z=feed,
        balances={"total": {"Methanol": 1.0, "Cyclohexane": 1.0}},
        totals={"total": 1.0},
        reactions=[reaction],
        phase_kind="lle_flash",
        phase_options=epcsaft.EquilibriumOptions(
            max_iterations=80,
            tolerance=1.0e-8,
            min_composition=1.0e-12,
            ipopt_linear_solver="mumps",
            ipopt_acceptable_tolerance=9.0e-7,
        ),
    )

    with pytest.raises(epcsaft.InputError, match="native Ipopt reactive phase-equilibrium NLP route"):
        mix.solve_equilibrium(problem)

    args = captured["args"]
    assert args[0] is mix._native
    assert args[1] == pytest.approx(298.15)
    assert args[2] == pytest.approx(1.013e5)
    assert args[3] == pytest.approx(feed.tolist())
    assert args[4] == 1
    assert args[5] == pytest.approx([1.0, 1.0])
    assert args[6] == pytest.approx([1.0])
    assert args[7] == 1
    assert args[8] == pytest.approx([-1.0, 1.0])
    assert args[9] == pytest.approx([-0.079259405371])
    assert captured["kwargs"]["linear_solver"] == "mumps"
    assert captured["kwargs"]["acceptable_tolerance"] == pytest.approx(9.0e-7)
    assert captured["kwargs"]["constraint_violation_tolerance"] == pytest.approx(1.0e-8)
    assert captured["kwargs"]["dual_infeasibility_tolerance"] == pytest.approx(1.0e-8)
    assert captured["kwargs"]["complementarity_tolerance"] == pytest.approx(1.0e-8)


def test_reactive_phase_equilibrium_problem_preserves_native_route_metadata_in_result_diagnostics(monkeypatch) -> None:
    mix, feed, reaction = _toy_reactive_phase_case()

    def accepted_native_route(*_args, **_kwargs):
        return {
            "backend": "ipopt",
            "compiled": True,
            "adapter_available": True,
            "ran": True,
            "accepted": True,
            "status": "accepted",
            "solver_accepted": True,
            "solver_status": "success",
            "application_status": "solve_succeeded",
            "problem_name": "reactive_liquid_root_eos",
            "variable_model": "log_phase_species_amounts_plus_log_density",
            "density_backend": "explicit_log_density_pressure_constraint",
            "residual_families": [
                "conserved_balance",
                "reaction_stationarity",
                "phase_equilibrium",
                "phase_charge",
            ],
            "constraint_families": [
                "conserved_balance",
                "phase_pressure_consistency",
                "phase_distance",
            ],
            "postsolve": {
                "accepted": True,
                "density_backend": "explicit_log_density_pressure_constraint",
                "phase_compositions": [[0.12, 0.88], [0.80, 0.20]],
                "phase_amount_totals": [0.5, 0.5],
                "phase_volumes": [0.01, 0.02],
                "phase_ln_fugacity_coefficients": [[0.0, 0.0], [0.0, 0.0]],
                "phase_distance": 0.68,
            },
        }

    monkeypatch.setattr(_core, "_native_reactive_lle_eos_route_result", accepted_native_route)
    problem = epcsaft.ReactivePhaseEquilibriumProblem(
        T=298.15,
        P=1.013e5,
        z=feed,
        balances={"total": {"Methanol": 1.0, "Cyclohexane": 1.0}},
        totals={"total": 1.0},
        reactions=[reaction],
        phase_kind="lle_flash",
    )

    result = mix.solve_equilibrium(problem)

    assert result.diagnostics["reactive_route_status"] == "accepted"
    assert result.diagnostics["variable_model"] == "log_phase_species_amounts_plus_log_density"
    assert result.diagnostics["residual_families"] == [
        "conserved_balance",
        "reaction_stationarity",
        "phase_equilibrium",
        "phase_charge",
    ]
    assert result.diagnostics["constraint_families"] == [
        "conserved_balance",
        "phase_pressure_consistency",
        "phase_distance",
    ]


def test_reactive_phase_equilibrium_problem_rejects_non_lle_production_kind() -> None:
    mix, feed, reaction = _toy_reactive_phase_case()
    problem = epcsaft.ReactivePhaseEquilibriumProblem(
        T=298.15,
        P=1.013e5,
        z=feed,
        balances={"total": {"Methanol": 1.0, "Cyclohexane": 1.0}},
        totals={"total": 1.0},
        reactions=[reaction],
        phase_kind="tp_flash",
    )

    with pytest.raises(epcsaft.InputError, match="reactive_staged_equilibrium"):
        mix.solve_equilibrium(problem)
