from __future__ import annotations

import pytest

import epcsaft


def _salt_mixture(x, T, P):
    _ = P
    return epcsaft.ePCSAFTMixture.from_dataset("2026_Khudaida", ["H2O", "Na+", "Cl-"], x, T)


def test_reactive_electrolyte_bubble_requires_native_speciation_request_before_bubble_route() -> None:
    with pytest.raises(epcsaft.InputError, match="requires at least one reaction"):
        epcsaft.solve_reactive_electrolyte_bubble(
            species=["H2O", "Na+", "Cl-"],
            mixture_factory=_salt_mixture,
            T=298.15,
            P=101325.0,
            balances={
                "water": {"H2O": 1.0},
                "sodium": {"Na+": 1.0},
                "chloride": {"Cl-": 1.0},
            },
            totals={"water": 0.98, "sodium": 0.01, "chloride": 0.01},
            reactions=[],
            initial_x=[0.98, 0.01, 0.01],
            vapor_species=["H2O"],
            volatile_species=["H2O"],
            nonvolatile_species=["Na+", "Cl-"],
        )
