from __future__ import annotations

import json

import numpy as np
import pytest

import epcsaft


def _hydrocarbon_mixture() -> epcsaft.ePCSAFTMixture:
    params = {
        "m": np.asarray([1.0, 1.6069, 2.0020]),
        "s": np.asarray([3.7039, 3.5206, 3.6184]),
        "e": np.asarray([150.03, 191.42, 208.11]),
        "k_ij": np.asarray(
            [
                [0.0, 3.0e-4, 1.15e-2],
                [3.0e-4, 0.0, 5.10e-3],
                [1.15e-2, 5.10e-3, 0.0],
            ]
        ),
    }
    return epcsaft.ePCSAFTMixture.from_params(params, species=["Methane", "Ethane", "Propane"])


def test_solve_reactive_staged_equilibrium_returns_chemical_and_phase_results() -> None:
    mix = _hydrocarbon_mixture()

    result = epcsaft.solve_reactive_staged_equilibrium(
        species=mix.species,
        mixture_factory=lambda x, T, P: mix,
        T=220.0,
        P=1.0e5,
        balances={
            "methane": {"Methane": 1.0},
            "ethane": {"Ethane": 1.0},
            "propane": {"Propane": 1.0},
        },
        totals={"methane": 0.1, "ethane": 0.3, "propane": 0.6},
        reactions=[],
        initial_x=[0.1, 0.3, 0.6],
        phase_kind="tp_flash",
    )

    assert isinstance(result, epcsaft.ReactiveStagedEquilibriumResult)
    assert result.success is True
    assert result.chemical.success is True
    assert isinstance(result.phase, epcsaft.EquilibriumResult)
    assert result.phase.problem_kind == "tp_flash"
    assert result.z == pytest.approx({"Methane": 0.1, "Ethane": 0.3, "Propane": 0.6})
    assert result.diagnostics["workflow"] == "chemical_equilibrium_then_phase_equilibrium"
    assert result.diagnostics["phase_kind"] == "tp_flash"
    assert "reactive_flash_tp" not in dir(mix)
    json.dumps(result.to_dict(), allow_nan=False)
