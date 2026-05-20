from __future__ import annotations

from dataclasses import fields

import numpy as np
import pytest

import epcsaft
from epcsaft import InputError, _core
from epcsaft.electrolyte_bubble import electrolyte_bubble_pressure


def _salt_mixture() -> epcsaft.ePCSAFTMixture:
    x = np.asarray([0.98, 0.01, 0.01], dtype=float)
    return epcsaft.ePCSAFTMixture.from_dataset("2026_Khudaida", ["H2O", "Na+", "Cl-"], x, 298.15)


class _FakeState:
    def fugacity_coefficient(self):
        return [0.0, 0.0, 0.0]


class _FakeElectrolyteBubbleMixture:
    species = ["H2O", "Na+", "Cl-"]
    _native = object()

    def state(self, *args, **kwargs):
        return _FakeState()


def test_electrolyte_bubble_pressure_builds_native_route_before_ipopt_gate(monkeypatch) -> None:
    assert {field.name for field in fields(epcsaft.ElectrolyteBubbleOptions)} == {
        "max_iterations",
        "tolerance",
        "min_composition",
        "charge_tolerance",
        "hessian_mode",
        "ipopt_iteration_history_limit",
        "ipopt_linear_solver",
        "ipopt_acceptable_tolerance",
        "ipopt_constraint_violation_tolerance",
        "ipopt_dual_infeasibility_tolerance",
        "ipopt_complementarity_tolerance",
        "continuation_state",
    }
    mix = _salt_mixture()
    calls: list[dict[str, object]] = []

    def fake_route(
        _native,
        temperature,
        liquid_composition,
        max_iterations,
        tolerance,
        hessian_mode,
        iteration_history_limit,
        phase_total_tolerance,
        pressure_tolerance,
        charge_tolerance,
        chemical_potential_tolerance,
        phase_distance_tolerance,
        continuation_state,
        **ipopt_controls,
    ):
        calls.append(
            {
                "temperature": temperature,
                "liquid_composition": liquid_composition,
                "max_iterations": max_iterations,
                "tolerance": tolerance,
                "continuation_state": continuation_state,
                "phase_total_tolerance": phase_total_tolerance,
                "pressure_tolerance": pressure_tolerance,
                "charge_tolerance": charge_tolerance,
                "chemical_potential_tolerance": chemical_potential_tolerance,
                "phase_distance_tolerance": phase_distance_tolerance,
                "ipopt_controls": ipopt_controls,
            }
        )
        return {
            "backend": "ipopt",
            "compiled": False,
            "ran": False,
            "accepted": False,
            "status": "ipopt_dependency_required",
            "postsolve": {"accepted": False},
        }

    monkeypatch.setattr(_core, "_native_electrolyte_bubble_p_eos_route_result", fake_route)

    with pytest.raises(InputError, match="native Ipopt equilibrium route builder"):
        mix.equilibrium(
            kind="electrolyte_bubble_pressure",
            T=298.15,
            x_liq=[0.98, 0.01, 0.01],
            vapor_species=["H2O"],
            volatile_species=["H2O"],
            nonvolatile_species=["Na+", "Cl-"],
            backend="native",
            options=epcsaft.ElectrolyteBubbleOptions(
                ipopt_linear_solver="mumps",
                ipopt_acceptable_tolerance=9.0e-7,
                ipopt_constraint_violation_tolerance=8.0e-8,
                ipopt_dual_infeasibility_tolerance=7.0e-8,
                ipopt_complementarity_tolerance=6.0e-8,
            ),
        )

    assert len(calls) == 1
    call = calls[0]
    assert call["temperature"] == pytest.approx(298.15)
    assert call["liquid_composition"] == pytest.approx([0.98, 0.01, 0.01])
    assert call["max_iterations"] == 80
    assert call["tolerance"] == pytest.approx(1.0e-6)
    assert call["phase_total_tolerance"] == pytest.approx(1.0e-6)
    assert call["pressure_tolerance"] == pytest.approx(0.1)
    assert call["charge_tolerance"] == pytest.approx(1.0e-8)
    assert call["chemical_potential_tolerance"] == pytest.approx(1.0e-6)
    assert call["phase_distance_tolerance"] > 0.0
    assert call["ipopt_controls"]["linear_solver"] == "mumps"
    assert call["ipopt_controls"]["acceptable_tolerance"] == pytest.approx(9.0e-7)
    assert call["ipopt_controls"]["constraint_violation_tolerance"] == pytest.approx(8.0e-8)
    assert call["ipopt_controls"]["dual_infeasibility_tolerance"] == pytest.approx(7.0e-8)
    assert call["ipopt_controls"]["complementarity_tolerance"] == pytest.approx(6.0e-8)


def test_electrolyte_bubble_pressure_reports_phase_eligibility_from_vapor_controls(monkeypatch) -> None:
    mix = _FakeElectrolyteBubbleMixture()

    def fake_route(*_args, **_kwargs):
        return {
            "backend": "ipopt",
            "compiled": True,
            "ran": True,
            "accepted": True,
            "status": "accepted",
            "variables": [1.0e5],
            "phase_amounts": [[0.98, 0.01, 0.01], [1.0, 0.0, 0.0]],
            "phase_volumes": [1.0e-4, 1.0],
            "postsolve": {"accepted": True, "charge_balance_norm": 0.0},
        }

    monkeypatch.setattr(_core, "_native_electrolyte_bubble_p_eos_route_result", fake_route)

    result = electrolyte_bubble_pressure(
        mix,
        T=298.15,
        x_liq=[0.98, 0.01, 0.01],
        vapor_species=["H2O"],
        nonvolatile_species=["Na+", "Cl-"],
    )

    diagnostics = result.diagnostics
    assert diagnostics["phase_eligibility_mask_available"] is True
    assert diagnostics["phase_eligibility_rows"] == 2
    assert diagnostics["phase_eligibility_cols"] == 3
    assert diagnostics["phase_eligibility_shape"] == [2, 3]
    assert diagnostics["phase_eligibility_phase_labels"] == ["liq", "vap"]
    assert diagnostics["phase_eligibility_species_labels"] == ["H2O", "Na+", "Cl-"]
    assert diagnostics["phase_eligibility_mask_source"] == "public_vapor_species"
    assert diagnostics["phase_eligibility_mask"] == pytest.approx([1.0, 1.0, 1.0, 1.0, 0.0, 0.0])


def test_electrolyte_bubble_pressure_rejects_unknown_nonvolatile_species() -> None:
    mix = _salt_mixture()

    with pytest.raises(InputError, match="Unknown nonvolatile species: K\\+"):
        mix.equilibrium(
            kind="electrolyte_bubble_pressure",
            T=298.15,
            x_liq=[0.98, 0.01, 0.01],
            vapor_species=["H2O"],
            nonvolatile_species=["K+"],
            backend="native",
        )


def test_electrolyte_bubble_pressure_rejects_python_backend_alias() -> None:
    mix = _salt_mixture()

    with pytest.raises(InputError, match="backend must be None or 'native'"):
        mix.equilibrium(
            kind="electrolyte_bubble_pressure",
            T=298.15,
            x_liq=[0.98, 0.01, 0.01],
            vapor_species=["H2O"],
            backend="python",
        )
