from __future__ import annotations

import json
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


def test_solve_reactive_speciation_returns_balanced_activity_coupled_state() -> None:
    species = ["H2O", "NaCl", "Na+", "Cl-"]
    mix = _salt_speciation_mixture()
    initial_x = np.asarray([0.998, 0.001, 0.0005, 0.0005], dtype=float)
    state = mix.state(T=298.15, P=1.0e5, x=initial_x, phase="liq")
    gamma = state.activity_coefficient(species=species)
    log_k = math.log(initial_x[2] * gamma["Na+"]) + math.log(initial_x[3] * gamma["Cl-"])
    log_k -= math.log(initial_x[1] * gamma["NaCl"])

    result = epcsaft.solve_reactive_speciation(
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
                log_equilibrium_constant=log_k,
            )
        ],
        initial_x=initial_x,
        options=epcsaft.ReactiveSpeciationOptions(max_iterations=8, tolerance=1.0e-8),
    )

    assert isinstance(result, epcsaft.ReactiveSpeciationResult)
    assert result.success is True
    assert sum(result.x.values()) == pytest.approx(1.0, abs=1.0e-10)
    assert abs(result.charge_residual) <= 1.0e-10
    assert max(abs(value) for value in result.mass_balance_residuals.values()) <= 1.0e-8
    assert max(abs(value) for value in result.reaction_residuals) <= 1.0e-8
    assert set(result.activity_coefficients) == set(species)
    assert result.state_failure_count == 0
    assert result.diagnostics["solver_language"] == "c++"
    assert result.diagnostics["backend"] == "native"
    assert result.diagnostics["native_entrypoint"] == "_solve_chemical_equilibrium_native"
    json.dumps(result.to_dict(), allow_nan=False)


def test_reactive_speciation_options_have_no_backend_selector() -> None:
    assert "backend" not in {field.name for field in fields(epcsaft.ReactiveSpeciationOptions)}
