from __future__ import annotations

import json
import math

import numpy as np
import pytest

import epcsaft


def _co2_water_salt_reaction_fixture() -> tuple[epcsaft.ePCSAFTMixture, list[str], float, np.ndarray]:
    species = ["CO2", "H2O", "NaCl", "Na+", "Cl-"]
    params = {
        "m": np.asarray([2.0729, 1.2047, 1.0, 1.0, 1.0]),
        "s": np.asarray([2.7852, 2.7927, 3.0, 2.8232, 2.7560]),
        "e": np.asarray([169.21, 353.95, 200.0, 230.0, 170.0]),
        "z": np.asarray([0.0, 0.0, 0.0, 1.0, -1.0]),
        "dielc": np.asarray([1.6, 78.09, 8.0, 8.0, 8.0]),
        "d_born": np.asarray([0.0, 0.0, 0.0, 3.445, 4.1]),
        "MW": np.asarray([44.0095e-3, 18.01528e-3, 58.44e-3, 22.989e-3, 35.45e-3]),
    }
    mix = epcsaft.ePCSAFTMixture.from_params(params, species=species)
    initial_x = np.asarray([0.02, 0.978, 0.001, 0.0005, 0.0005], dtype=float)
    state = mix.state(T=313.15, P=1.0e5, x=initial_x, phase="liq")
    gamma = state.activity_coefficient(species=species)
    log_k = math.log(initial_x[3] * gamma["Na+"]) + math.log(initial_x[4] * gamma["Cl-"])
    log_k -= math.log(initial_x[2] * gamma["NaCl"])
    return mix, species, log_k, initial_x


def _balances() -> dict[str, dict[str, float]]:
    return {
        "co2_total": {"CO2": 1.0},
        "water_total": {"H2O": 1.0},
        "sodium_total": {"NaCl": 1.0, "Na+": 1.0},
        "chloride_total": {"NaCl": 1.0, "Cl-": 1.0},
    }


def _totals(x: np.ndarray) -> dict[str, float]:
    return {
        "co2_total": float(x[0]),
        "water_total": float(x[1]),
        "sodium_total": float(x[2] + x[3]),
        "chloride_total": float(x[2] + x[4]),
    }


def test_solve_reactive_electrolyte_bubble_returns_structured_result() -> None:
    mix, species, log_k, initial_x = _co2_water_salt_reaction_fixture()

    result = epcsaft.solve_reactive_electrolyte_bubble(
        species=species,
        mixture_factory=lambda x, T, P: mix,
        T=313.15,
        P_seed=1.0e5,
        balances=_balances(),
        totals=_totals(initial_x),
        reactions=[
            epcsaft.ReactionDefinition(
                {"NaCl": -1.0, "Na+": 1.0, "Cl-": 1.0},
                log_equilibrium_constant=log_k,
                name="salt_dissociation",
            )
        ],
        initial_x=initial_x,
        vapor_species=["CO2", "H2O"],
        nonvolatile_species=["NaCl", "Na+", "Cl-"],
        options=epcsaft.ReactiveElectrolyteBubbleOptions(
            speciation_options=epcsaft.ReactiveSpeciationOptions(max_iterations=8, tolerance=1.0e-8),
            bubble_options=epcsaft.ElectrolyteBubbleOptions(initial_pressure=1.0e5, max_iterations=80),
        ),
    )

    assert isinstance(result, epcsaft.ReactiveElectrolyteBubbleResult)
    assert result.success is True
    assert sum(result.x_liq.values()) == pytest.approx(1.0, abs=1.0e-10)
    assert result.P_total > 0.0
    assert result.partial_pressures["CO2"] > 0.0
    assert set(result.y_vap) == {"CO2", "H2O"}
    assert set(result.named_reaction_residuals) == {"salt_dissociation"}
    assert result.penalty_residuals == []
    json.dumps(result.to_dict(), allow_nan=False)


def test_reactive_electrolyte_bubble_result_mode_keeps_fixed_shape_on_failure() -> None:
    mix, species, log_k, initial_x = _co2_water_salt_reaction_fixture()

    result = epcsaft.solve_reactive_electrolyte_bubble(
        species=species,
        mixture_factory=lambda x, T, P: mix,
        T=313.15,
        P_seed=1.0e5,
        balances=_balances(),
        totals=_totals(initial_x),
        reactions=[
            epcsaft.ReactionDefinition(
                {"NaCl": -1.0, "Na+": 1.0, "Cl-": 1.0},
                log_equilibrium_constant=log_k,
                name="salt_dissociation",
            )
        ],
        initial_x=initial_x,
        vapor_species=["CO2", "H2O"],
        nonvolatile_species=["NaCl", "Na+", "Cl-"],
        options=epcsaft.ReactiveElectrolyteBubbleOptions(
            bubble_options=epcsaft.ElectrolyteBubbleOptions(
                initial_pressure=1.0e5,
                max_iterations=1,
                tolerance=1.0e-20,
            ),
            error_mode="result",
        ),
    )

    assert result.success is False
    assert result.P_total > 0.0
    assert set(result.x_liq) == set(species)
    assert len(result.penalty_residuals) == 1
    assert result.diagnostics["failure_stage"] in {"nonconverged_subsolve", "electrolyte_bubble_pressure"}
    json.dumps(result.to_dict(), allow_nan=False)


def test_reactive_electrolyte_bubble_sweep_reuses_continuation_seeds() -> None:
    mix, species, log_k, initial_x = _co2_water_salt_reaction_fixture()
    next_x = np.asarray([0.021, 0.977, 0.001, 0.0005, 0.0005], dtype=float)

    results = mix.equilibrium_sweep(
        kind="reactive_electrolyte_bubble_pressure",
        points=[
            {"T": 313.15, "P_seed": 1.0e5, "totals": _totals(initial_x), "initial_x": initial_x},
            {"T": 313.15, "totals": _totals(next_x), "initial_x": next_x},
        ],
        balances=_balances(),
        reactions=[
            epcsaft.ReactionDefinition(
                {"NaCl": -1.0, "Na+": 1.0, "Cl-": 1.0},
                log_equilibrium_constant=log_k,
                name="salt_dissociation",
            )
        ],
        vapor_species=["CO2", "H2O"],
        nonvolatile_species=["NaCl", "Na+", "Cl-"],
        options=epcsaft.ReactiveElectrolyteBubbleOptions(
            speciation_options=epcsaft.ReactiveSpeciationOptions(max_iterations=8, tolerance=1.0e-8),
            bubble_options=epcsaft.ElectrolyteBubbleOptions(initial_pressure=1.0e5, max_iterations=80),
            error_mode="result",
        ),
    )

    assert len(results) == 2
    assert all(isinstance(result, epcsaft.ReactiveElectrolyteBubbleResult) for result in results)
    assert results[0].success is True
    assert results[1].diagnostics["continuation_used"] is True
    assert results[1].P_total > 0.0
