from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from epcsaft import InputError, SolutionError


def _salt_mixture() -> epcsaft.ePCSAFTMixture:
    x = np.asarray([0.98, 0.01, 0.01], dtype=float)
    return epcsaft.ePCSAFTMixture.from_dataset("2026_Khudaida", ["H2O", "Na+", "Cl-"], x, 298.15)


def test_electrolyte_bubble_pressure_converges_for_water_salt() -> None:
    mix = _salt_mixture()

    result = mix.equilibrium(
        kind="electrolyte_bubble_pressure",
        T=298.15,
        x_liq=[0.98, 0.01, 0.01],
        vapor_species=["H2O"],
        volatile_species=["H2O"],
        nonvolatile_species=["Na+", "Cl-"],
        backend="native",
    )

    assert result.success
    assert np.isfinite(result.P)
    assert result.P > 0.0
    assert result.y_vap == pytest.approx({"H2O": 1.0})
    assert result.charge_residual == pytest.approx(0.0, abs=1e-8)
    assert result.fugacity_residual_norm <= result.diagnostics["acceptance_tolerance"]
    assert result.diagnostics["accepted_by_diagnostic_envelope"]
    assert result.diagnostics["native_entrypoint"] == "_solve_electrolyte_bubble_native"


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


def test_electrolyte_bubble_pressure_rejects_ionic_vapor_species() -> None:
    mix = _salt_mixture()

    with pytest.raises(InputError, match="vapor_species must be neutral"):
        mix.equilibrium(
            kind="electrolyte_bubble_pressure",
            T=298.15,
            x_liq=[0.98, 0.01, 0.01],
            vapor_species=["Na+"],
            backend="native",
        )


def test_electrolyte_bubble_pressure_best_effort_returns_diagnostics() -> None:
    mix = _salt_mixture()
    options = epcsaft.ElectrolyteBubbleOptions(
        min_pressure=9.0e7,
        max_pressure=1.0e8,
        initial_pressure=1.0e8,
        max_bracket_expansions=0,
        max_iterations=1,
        return_best_effort=True,
    )

    result = mix.equilibrium(
        kind="electrolyte_bubble_pressure",
        T=298.15,
        x_liq=[0.98, 0.01, 0.01],
        vapor_species=["H2O"],
        options=options,
    )

    assert not result.success
    assert np.isfinite(result.P)
    assert result.diagnostics["best_P"] == pytest.approx(result.P)
    assert "best_objective" in result.diagnostics
    assert "state_failure_count" in result.diagnostics


def test_electrolyte_bubble_pressure_nonconvergence_raises_by_default() -> None:
    mix = _salt_mixture()
    options = epcsaft.ElectrolyteBubbleOptions(
        min_pressure=9.0e7,
        max_pressure=1.0e8,
        initial_pressure=1.0e8,
        max_bracket_expansions=0,
        max_iterations=1,
    )

    with pytest.raises(SolutionError) as excinfo:
        mix.equilibrium(
            kind="electrolyte_bubble_pressure",
            T=298.15,
            x_liq=[0.98, 0.01, 0.01],
            vapor_species=["H2O"],
            options=options,
        )
    assert "best_P" in excinfo.value.diagnostics
