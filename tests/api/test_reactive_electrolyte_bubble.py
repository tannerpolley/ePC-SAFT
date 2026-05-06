from __future__ import annotations

import pytest

import epcsaft
from epcsaft import InputError


def test_reactive_electrolyte_bubble_requires_native_backend() -> None:
    with pytest.raises(InputError, match="native C\\+\\+ backend"):
        epcsaft.solve_reactive_electrolyte_bubble(
            species=["H2O", "Na+", "Cl-"],
            mixture_factory=lambda x, T, P: None,
            T=298.15,
            P_seed=101325.0,
            balances={"water": {"H2O": 1.0}},
            totals={"water": 1.0},
            reactions=[],
            initial_x=[0.98, 0.01, 0.01],
            vapor_species=["H2O"],
        )


def test_reactive_electrolyte_bubble_sweep_requires_native_backend() -> None:
    with pytest.raises(InputError, match="native C\\+\\+ backend"):
        epcsaft.solve_reactive_electrolyte_bubble_sweep(
            species=["H2O", "Na+", "Cl-"],
            mixture_factory=lambda x, T, P: None,
            points=[{"T": 298.15, "totals": {"water": 1.0}, "initial_x": [0.98, 0.01, 0.01]}],
            balances={"water": {"H2O": 1.0}},
            reactions=[],
            vapor_species=["H2O"],
        )
