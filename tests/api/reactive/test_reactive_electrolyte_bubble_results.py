from __future__ import annotations

import inspect

import pytest

import epcsaft
from epcsaft.electrolyte_bubble import ElectrolyteBubbleResult
from epcsaft.reactive_speciation import ReactiveSpeciationResult


def _salt_mixture(x, T, P):
    _ = P
    return epcsaft.ePCSAFTMixture.from_dataset("2026_Khudaida", ["H2O", "Na+", "Cl-"], x, T)

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
            P=101325.0,
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
        P=101325.0,
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

def test_reactive_electrolyte_bubble_sweep_does_not_continue_bubble_seed_controls(monkeypatch) -> None:
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
            bubble_options=epcsaft.ElectrolyteBubbleOptions(),
            error_mode="result",
        ),
    )

    assert [result.success for result in results] == [True, False]
    assert len(calls) == 2
    assert results[1].P_total == pytest.approx(111000.0)

def test_reactive_electrolyte_bubble_sweeps_expose_no_continuation_flag() -> None:
    assert "continuation" not in inspect.signature(epcsaft.solve_reactive_electrolyte_bubble_sweep).parameters
    assert "continuation" not in inspect.signature(epcsaft.ePCSAFTMixture.equilibrium_sweep).parameters


def test_reactive_electrolyte_bubble_public_api_rejects_pressure_seed_key() -> None:
    removed = "P" + "_seed"

    assert removed not in inspect.signature(epcsaft.solve_reactive_electrolyte_bubble).parameters
    assert removed not in inspect.signature(epcsaft.ePCSAFTMixture.reactive_electrolyte_bubble_p).parameters
    assert removed not in epcsaft.ReactiveElectrolyteBubbleProblem.__dataclass_fields__
    with pytest.raises(epcsaft.InputError, match=removed):
        epcsaft.solve_reactive_electrolyte_bubble_sweep(
            species=["H2O", "Na+", "Cl-"],
            mixture_factory=_salt_mixture,
            points=[{"T": 298.15, removed: 101325.0, "totals": {"water": 0.98}}],
            balances={"water": {"H2O": 1.0}},
            reactions=[],
            vapor_species=["H2O"],
        )


def test_reactive_electrolyte_bubble_sweep_honors_point_options(monkeypatch) -> None:
    chemical = _successful_chemical_result()
    calls = []
    monkeypatch.setattr("epcsaft.reactive_electrolyte.solve_reactive_speciation", lambda **kwargs: chemical)

    def bubble_sequence(*args, **kwargs):
        _ = args
        calls.append(kwargs["options"])
        return ElectrolyteBubbleResult(
            success=True,
            message="converged",
            P=101325.0,
            y_vap={"H2O": 1.0},
            x_liq=[0.98, 0.01, 0.01],
            ln_phi_liq={"H2O": 0.0},
            ln_phi_vap={"H2O": 0.0},
            fugacity_residual={"H2O": 1.0e-8},
            fugacity_residual_norm=1.0e-8,
            charge_residual=0.0,
            partial_pressures={"H2O": 101325.0},
            diagnostics={"state_failure_count": 0},
        )

    monkeypatch.setattr("epcsaft.reactive_electrolyte.electrolyte_bubble_pressure", bubble_sequence)

    point_options = epcsaft.ReactiveElectrolyteBubbleOptions(
        bubble_options=epcsaft.ElectrolyteBubbleOptions(max_iterations=7),
        error_mode="result",
    )
    epcsaft.solve_reactive_electrolyte_bubble_sweep(
        species=["H2O", "Na+", "Cl-"],
        mixture_factory=_salt_mixture,
        points=[
            {
                "T": 298.15,
                "totals": {"water": 0.98},
                "initial_x": [0.98, 0.01, 0.01],
                "options": point_options,
            }
        ],
        balances={"water": {"H2O": 1.0}},
        reactions=[],
        vapor_species=["H2O"],
        options=epcsaft.ReactiveElectrolyteBubbleOptions(
            bubble_options=epcsaft.ElectrolyteBubbleOptions(),
            error_mode="result",
        ),
    )

    assert calls[0].max_iterations == 7
