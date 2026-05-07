from __future__ import annotations

import pytest

import epcsaft
from epcsaft import reactive_electrolyte
from epcsaft.electrolyte_bubble import ElectrolyteBubbleResult
from epcsaft.reactive_speciation import ReactiveSpeciationResult


def _salt_mixture(x, T, P):
    _ = P
    return epcsaft.ePCSAFTMixture.from_dataset("2026_Khudaida", ["H2O", "Na+", "Cl-"], x, T)


def test_reactive_electrolyte_bubble_runs_native_speciation_then_bubble() -> None:
    result = epcsaft.solve_reactive_electrolyte_bubble(
        species=["H2O", "Na+", "Cl-"],
        mixture_factory=_salt_mixture,
        T=298.15,
        P_seed=101325.0,
        balances={
            "water": {"H2O": 1.0},
            "sodium": {"Na+": 1.0},
            "chloride": {"Cl-": 1.0},
        },
        totals={"water": 0.98, "sodium": 0.01, "chloride": 0.01},
        reactions=[],
        initial_x=[0.98, 0.01, 0.01],
        vapor_species=["H2O"],
    )

    assert result.success
    assert result.P_total > 0.0
    assert result.y_vap == pytest.approx({"H2O": 1.0})
    assert result.named_reaction_residuals == {}
    assert (
        result.diagnostics["native_entrypoint"]
        == "_solve_chemical_equilibrium_native_then__solve_electrolyte_bubble_native"
    )


def test_reactive_electrolyte_bubble_sweep_uses_continuation() -> None:
    results = epcsaft.solve_reactive_electrolyte_bubble_sweep(
        species=["H2O", "Na+", "Cl-"],
        mixture_factory=_salt_mixture,
        points=[
            {"T": 298.15, "totals": {"water": 0.98, "sodium": 0.01, "chloride": 0.01}, "initial_x": [0.98, 0.01, 0.01]},
            {
                "T": 298.15,
                "totals": {"water": 0.982, "sodium": 0.009, "chloride": 0.009},
                "initial_x": [0.982, 0.009, 0.009],
            },
        ],
        balances={
            "water": {"H2O": 1.0},
            "sodium": {"Na+": 1.0},
            "chloride": {"Cl-": 1.0},
        },
        reactions=[],
        vapor_species=["H2O"],
    )

    assert len(results) == 2
    assert all(result.success for result in results)


def test_reactive_electrolyte_bubble_accepts_phase_handoff_speciation_residuals(monkeypatch) -> None:
    chemical = ReactiveSpeciationResult(
        success=False,
        message="reactive speciation residual family tolerances were not met",
        x={"H2O": 0.98, "Na+": 0.01, "Cl-": 0.01},
        activity_coefficients={"H2O": 1.0, "Na+": 1.0, "Cl-": 1.0},
        mass_balance_residuals={"water": 0.0, "sodium": 0.0, "chloride": 0.0},
        charge_residual=1.0e-14,
        reaction_residuals=[3.7e-6],
        named_reaction_residuals={"toy": 3.7e-6},
        state_failure_count=0,
        diagnostics={
            "native_success": False,
            "mass_residual_norm": 0.0,
            "charge_residual_abs": 1.0e-14,
            "reaction_residual_norm": 3.7e-6,
        },
    )
    bubble = ElectrolyteBubbleResult(
        success=True,
        message="converged",
        P=101325.0,
        y_vap={"H2O": 1.0},
        x_liq=[0.98, 0.01, 0.01],
        ln_phi_liq={"H2O": 0.0},
        ln_phi_vap={"H2O": 0.0},
        fugacity_residual={"H2O": 1.0e-7},
        fugacity_residual_norm=1.0e-7,
        charge_residual=1.0e-14,
        partial_pressures={"H2O": 101325.0},
        diagnostics={"state_failure_count": 0},
    )

    monkeypatch.setattr("epcsaft.reactive_electrolyte.solve_reactive_speciation", lambda **kwargs: chemical)
    monkeypatch.setattr("epcsaft.reactive_electrolyte.electrolyte_bubble_pressure", lambda *args, **kwargs: bubble)

    result = epcsaft.solve_reactive_electrolyte_bubble(
        species=["H2O", "Na+", "Cl-"],
        mixture_factory=_salt_mixture,
        T=298.15,
        P_seed=101325.0,
        balances={"water": {"H2O": 1.0}},
        totals={"water": 0.98},
        reactions=[],
        initial_x=[0.98, 0.01, 0.01],
        vapor_species=["H2O"],
    )

    assert result.success
    assert result.message == "converged"
    assert not result.diagnostics["speciation_strict_success"]
    assert result.diagnostics["speciation_phase_handoff_success"]
    assert not result.diagnostics["speciation_phase_handoff"]["native_success"]
    assert result.diagnostics["speciation_phase_handoff"]["reason"] == "residuals_within_phase_handoff_tolerances"


def test_reactive_electrolyte_bubble_respects_configured_phase_handoff_tolerances(monkeypatch) -> None:
    chemical = ReactiveSpeciationResult(
        success=False,
        message="reactive speciation residual family tolerances were not met",
        x={"H2O": 0.98, "Na+": 0.01, "Cl-": 0.01},
        activity_coefficients={"H2O": 1.0, "Na+": 1.0, "Cl-": 1.0},
        mass_balance_residuals={"water": 0.0, "sodium": 0.0, "chloride": 0.0},
        charge_residual=1.0e-14,
        reaction_residuals=[3.7e-6],
        named_reaction_residuals={"toy": 3.7e-6},
        state_failure_count=0,
        diagnostics={
            "native_success": False,
            "mass_residual_norm": 0.0,
            "charge_residual_abs": 1.0e-14,
            "reaction_residual_norm": 3.7e-6,
        },
    )
    bubble = ElectrolyteBubbleResult(
        success=True,
        message="converged",
        P=101325.0,
        y_vap={"H2O": 1.0},
        x_liq=[0.98, 0.01, 0.01],
        ln_phi_liq={"H2O": 0.0},
        ln_phi_vap={"H2O": 0.0},
        fugacity_residual={"H2O": 1.0e-7},
        fugacity_residual_norm=1.0e-7,
        charge_residual=1.0e-14,
        partial_pressures={"H2O": 101325.0},
        diagnostics={"state_failure_count": 0},
    )

    monkeypatch.setattr("epcsaft.reactive_electrolyte.solve_reactive_speciation", lambda **kwargs: chemical)
    monkeypatch.setattr("epcsaft.reactive_electrolyte.electrolyte_bubble_pressure", lambda *args, **kwargs: bubble)

    result = epcsaft.solve_reactive_electrolyte_bubble(
        species=["H2O", "Na+", "Cl-"],
        mixture_factory=_salt_mixture,
        T=298.15,
        P_seed=101325.0,
        balances={"water": {"H2O": 1.0}},
        totals={"water": 0.98},
        reactions=[],
        initial_x=[0.98, 0.01, 0.01],
        vapor_species=["H2O"],
        options=epcsaft.ReactiveElectrolyteBubbleOptions(
            phase_handoff_reaction_tolerance=1.0e-7,
        ),
    )

    assert not result.success
    assert result.message == "reactive electrolyte speciation did not meet phase-handoff tolerances"
    assert not result.diagnostics["speciation_phase_handoff_success"]
    handoff = result.diagnostics["speciation_phase_handoff"]
    assert handoff["reason"] == "residuals_exceed_phase_handoff_tolerances"
    assert handoff["reaction_tolerance"] == pytest.approx(1.0e-7)


def test_reactive_electrolyte_bubble_sweep_preserves_phase_handoff_tolerances(monkeypatch) -> None:
    calls = []
    chemical = ReactiveSpeciationResult(
        success=False,
        message="reactive speciation residual family tolerances were not met",
        x={"H2O": 0.98, "Na+": 0.01, "Cl-": 0.01},
        activity_coefficients={"H2O": 1.0, "Na+": 1.0, "Cl-": 1.0},
        mass_balance_residuals={"water": 0.0, "sodium": 0.0, "chloride": 0.0},
        charge_residual=1.0e-14,
        reaction_residuals=[3.7e-6],
        named_reaction_residuals={"toy": 3.7e-6},
        state_failure_count=0,
        diagnostics={
            "native_success": False,
            "mass_residual_norm": 0.0,
            "charge_residual_abs": 1.0e-14,
            "reaction_residual_norm": 3.7e-6,
        },
    )
    bubble = ElectrolyteBubbleResult(
        success=True,
        message="converged",
        P=101325.0,
        y_vap={"H2O": 1.0},
        x_liq=[0.98, 0.01, 0.01],
        ln_phi_liq={"H2O": 0.0},
        ln_phi_vap={"H2O": 0.0},
        fugacity_residual={"H2O": 1.0e-7},
        fugacity_residual_norm=1.0e-7,
        charge_residual=1.0e-14,
        partial_pressures={"H2O": 101325.0},
        diagnostics={"state_failure_count": 0},
    )

    monkeypatch.setattr("epcsaft.reactive_electrolyte.solve_reactive_speciation", lambda **kwargs: chemical)
    monkeypatch.setattr("epcsaft.reactive_electrolyte.electrolyte_bubble_pressure", lambda *args, **kwargs: bubble)

    original = reactive_electrolyte.solve_reactive_electrolyte_bubble

    def wrapped_solve(**kwargs):
        calls.append(kwargs["options"])
        return original(**kwargs)

    monkeypatch.setattr("epcsaft.reactive_electrolyte.solve_reactive_electrolyte_bubble", wrapped_solve)

    results = epcsaft.solve_reactive_electrolyte_bubble_sweep(
        species=["H2O", "Na+", "Cl-"],
        mixture_factory=_salt_mixture,
        points=[
            {"T": 298.15, "totals": {"water": 0.98}, "initial_x": [0.98, 0.01, 0.01]},
            {"T": 298.15, "totals": {"water": 0.98}, "initial_x": [0.98, 0.01, 0.01]},
        ],
        balances={"water": {"H2O": 1.0}},
        reactions=[],
        vapor_species=["H2O"],
        options=epcsaft.ReactiveElectrolyteBubbleOptions(
            bubble_options=epcsaft.ElectrolyteBubbleOptions(initial_pressure=101325.0),
            phase_handoff_reaction_tolerance=1.0e-7,
        ),
    )

    assert [result.success for result in results] == [False, False]
    assert len(calls) == 2
    assert all(call.phase_handoff_reaction_tolerance == pytest.approx(1.0e-7) for call in calls)


def test_reactive_electrolyte_bubble_strict_bubble_failure_raises(monkeypatch) -> None:
    chemical = _successful_chemical_result()
    monkeypatch.setattr("epcsaft.reactive_electrolyte.solve_reactive_speciation", lambda **kwargs: chemical)

    def fail_bubble(*args, **kwargs):
        _ = args, kwargs
        raise epcsaft.SolutionError(
            "bubble failed",
            {
                "best_P": 88000.0,
                "best_y_vap": [1.0],
                "best_partial_pressures": [88000.0],
                "best_fugacity_residual_norm": 2.0e-4,
            },
        )

    monkeypatch.setattr("epcsaft.reactive_electrolyte.electrolyte_bubble_pressure", fail_bubble)

    with pytest.raises(epcsaft.SolutionError, match="bubble failed"):
        epcsaft.solve_reactive_electrolyte_bubble(
            species=["H2O", "Na+", "Cl-"],
            mixture_factory=_salt_mixture,
            T=298.15,
            P_seed=101325.0,
            balances={"water": {"H2O": 1.0}},
            totals={"water": 0.98},
            reactions=[],
            initial_x=[0.98, 0.01, 0.01],
            vapor_species=["H2O"],
        )


def test_reactive_electrolyte_bubble_result_mode_returns_bubble_failure(monkeypatch) -> None:
    chemical = _successful_chemical_result()
    monkeypatch.setattr("epcsaft.reactive_electrolyte.solve_reactive_speciation", lambda **kwargs: chemical)

    def fail_bubble(*args, **kwargs):
        _ = args, kwargs
        raise epcsaft.SolutionError(
            "bubble failed",
            {
                "best_P": 88000.0,
                "best_y_vap": [1.0],
                "best_partial_pressures": [88000.0],
                "best_fugacity_residual_norm": 2.0e-4,
                "state_failure_count": 3,
            },
        )

    monkeypatch.setattr("epcsaft.reactive_electrolyte.electrolyte_bubble_pressure", fail_bubble)

    result = epcsaft.solve_reactive_electrolyte_bubble(
        species=["H2O", "Na+", "Cl-"],
        mixture_factory=_salt_mixture,
        T=298.15,
        P_seed=101325.0,
        balances={"water": {"H2O": 1.0}},
        totals={"water": 0.98},
        reactions=[],
        initial_x=[0.98, 0.01, 0.01],
        vapor_species=["H2O"],
        options=epcsaft.ReactiveElectrolyteBubbleOptions(error_mode="result"),
    )

    assert result.success is False
    assert result.message == "bubble failed"
    assert result.P_total == pytest.approx(88000.0)
    assert result.y_vap == pytest.approx({"H2O": 1.0})
    assert result.partial_pressures == pytest.approx({"H2O": 88000.0})
    assert result.fugacity_residual_norm == pytest.approx(2.0e-4)
    assert result.diagnostics["speciation_strict_success"] is True
    assert result.diagnostics["speciation_phase_handoff_success"] is True
    assert result.diagnostics["bubble_success"] is False
    assert result.diagnostics["bubble_failure"]["message"] == "bubble failed"
    assert result.diagnostics["bubble"]["diagnostics"]["state_failure_count"] == 3


def test_reactive_electrolyte_bubble_sweep_continuation_survives_result_mode_failure(monkeypatch) -> None:
    chemical = _successful_chemical_result()
    calls = []
    monkeypatch.setattr("epcsaft.reactive_electrolyte.solve_reactive_speciation", lambda **kwargs: chemical)

    def bubble_sequence(*args, **kwargs):
        _ = args
        calls.append(kwargs["options"])
        if len(calls) == 1:
            return ElectrolyteBubbleResult(
                success=True,
                message="converged",
                P=112000.0,
                y_vap={"H2O": 1.0},
                x_liq=[0.98, 0.01, 0.01],
                ln_phi_liq={"H2O": 0.0},
                ln_phi_vap={"H2O": 0.0},
                fugacity_residual={"H2O": 1.0e-8},
                fugacity_residual_norm=1.0e-8,
                charge_residual=0.0,
                partial_pressures={"H2O": 112000.0},
                diagnostics={"state_failure_count": 0},
            )
        raise epcsaft.SolutionError(
            "bubble failed",
            {
                "best_P": 111000.0,
                "best_y_vap": [1.0],
                "best_partial_pressures": [111000.0],
                "best_fugacity_residual_norm": 1.0e-3,
            },
        )

    monkeypatch.setattr("epcsaft.reactive_electrolyte.electrolyte_bubble_pressure", bubble_sequence)

    results = epcsaft.solve_reactive_electrolyte_bubble_sweep(
        species=["H2O", "Na+", "Cl-"],
        mixture_factory=_salt_mixture,
        points=[
            {"T": 298.15, "totals": {"water": 0.98}, "initial_x": [0.98, 0.01, 0.01]},
            {"T": 298.15, "totals": {"water": 0.98}, "initial_x": [0.98, 0.01, 0.01]},
        ],
        balances={"water": {"H2O": 1.0}},
        reactions=[],
        vapor_species=["H2O"],
        options=epcsaft.ReactiveElectrolyteBubbleOptions(
            bubble_options=epcsaft.ElectrolyteBubbleOptions(initial_pressure=101325.0),
            error_mode="result",
        ),
    )

    assert [result.success for result in results] == [True, False]
    assert len(calls) == 2
    assert calls[1].initial_pressure == pytest.approx(112000.0)
    assert calls[1].initial_y_vap == pytest.approx({"H2O": 1.0})
    assert results[1].P_total == pytest.approx(111000.0)


def _successful_chemical_result() -> ReactiveSpeciationResult:
    return ReactiveSpeciationResult(
        success=True,
        message="converged",
        x={"H2O": 0.98, "Na+": 0.01, "Cl-": 0.01},
        activity_coefficients={"H2O": 1.0, "Na+": 1.0, "Cl-": 1.0},
        mass_balance_residuals={"water": 0.0, "sodium": 0.0, "chloride": 0.0},
        charge_residual=0.0,
        reaction_residuals=[],
        named_reaction_residuals={},
        state_failure_count=0,
        diagnostics={
            "native_success": True,
            "mass_residual_norm": 0.0,
            "charge_residual_abs": 0.0,
            "reaction_residual_norm": 0.0,
            "state_failure_count": 0,
        },
    )
