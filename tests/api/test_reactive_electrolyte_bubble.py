from __future__ import annotations

import pytest

import epcsaft


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
