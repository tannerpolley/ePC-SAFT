from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

import epcsaft

FIXTURE_ROOT = Path(__file__).resolve().parents[2] / "fixtures" / "literature" / "pure_neutral"


def _salt_mixture(x, T, P):
    _ = P
    return epcsaft.ePCSAFTMixture.from_dataset("2026_Khudaida", ["H2O", "Na+", "Cl-"], x, T)


def _generic_amine_co2_water_mixture() -> epcsaft.ePCSAFTMixture:
    fixture = json.loads((FIXTURE_ROOT / "mea_co2_h2o_benchmark.json").read_text(encoding="utf-8"))
    amine_key, amine_h_key, carbamate_key, bicarbonate_key = fixture["components"]
    reference = fixture["reference_values"]
    params = {
        "m": np.asarray(
            [
                1.2047,
                reference[amine_key]["m"],
                2.0729,
                1.0,
                1.0,
                1.0,
                1.0,
            ]
        ),
        "s": np.asarray(
            [
                2.7927,
                reference[amine_key]["s"],
                2.7852,
                reference[bicarbonate_key]["s"],
                reference[amine_h_key]["s"],
                reference[carbamate_key]["s"],
                2.0,
            ]
        ),
        "e": np.asarray(
            [
                353.95,
                reference[amine_key]["e"],
                169.21,
                reference[bicarbonate_key]["e"],
                reference[amine_h_key]["e"],
                reference[carbamate_key]["e"],
                120.0,
            ]
        ),
        "e_assoc": np.asarray(
            [
                2425.67,
                reference[amine_key]["e_assoc"],
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ]
        ),
        "vol_a": np.asarray(
            [
                0.0451,
                reference[amine_key]["vol_a"],
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ]
        ),
        "assoc_scheme": ["4C", "2B", None, None, None, None, None],
        "z": np.asarray([0.0, 0.0, 0.0, -1.0, 1.0, -1.0, 1.0]),
        "dielc": np.asarray([78.09, 35.0, 12.0, 20.0, 25.0, 22.0, 8.0]),
        "d_born": np.asarray(
            [
                0.0,
                0.0,
                0.0,
                reference[bicarbonate_key]["d_born"],
                reference[amine_h_key]["d_born"],
                reference[carbamate_key]["d_born"],
                2.0,
            ]
        ),
        "MW": np.asarray([18.01528e-3, 61.08e-3, 44.01e-3, 61.02e-3, 62.09e-3, 104.1e-3, 1.008e-3]),
    }
    return epcsaft.ePCSAFTMixture.from_params(
        params,
        species=["H2O", "Amine", "CO2", "HCO3-", "AmineH+", "AmineCOO-", "H+"],
    )


def _log_k_from_state(
    mix: epcsaft.ePCSAFTMixture,
    *,
    T: float,
    P: float,
    x: np.ndarray,
    stoichiometry: dict[str, float],
) -> float:
    state = mix.state(T=T, P=P, x=x, phase="liq")
    gamma = state.activity_coefficient(species=mix.species)
    return float(
        sum(
            nu * math.log(max(x[mix.species.index(label)] * gamma[label], 1.0e-300))
            for label, nu in stoichiometry.items()
        )
    )


def test_reactive_electrolyte_bubble_reports_volatile_neutral_partial_pressure() -> None:
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
        volatile_species=["H2O"],
        nonvolatile_species=["Na+", "Cl-"],
    )

    assert result.success is True
    assert result.P_total > 0.0
    assert result.y_vap == pytest.approx({"H2O": 1.0})
    assert result.partial_pressures["H2O"] == pytest.approx(result.P_total)
    assert "Na+" not in result.y_vap
    assert "Cl-" not in result.y_vap
    assert result.fugacity_residual_norm <= result.diagnostics["bubble"]["diagnostics"]["acceptance_tolerance"]
    assert result.diagnostics["speciation"]["diagnostics"]["solver_language"] == "c++"
    assert result.diagnostics["bubble"]["diagnostics"]["native_entrypoint"] == "_solve_electrolyte_bubble_native"


def test_generic_co2_amine_water_benchmark_reports_native_speciation_pressure_path() -> None:
    mix = _generic_amine_co2_water_mixture()
    target_x = np.asarray([0.7899, 0.19, 0.0001, 0.005, 0.005, 0.005, 0.005], dtype=float)
    initial_x = np.asarray([0.8, 0.18, 0.001, 0.004, 0.006, 0.004, 0.005], dtype=float)
    reactions = [
        {"CO2": -1.0, "H2O": -1.0, "HCO3-": 1.0, "H+": 1.0},
        {"Amine": -1.0, "H+": -1.0, "AmineH+": 1.0},
        {"Amine": -1.0, "CO2": -1.0, "AmineCOO-": 1.0, "H+": 1.0},
    ]
    reaction_definitions = [
        epcsaft.ReactionDefinition(
            stoichiometry,
            _log_k_from_state(mix, T=313.15, P=1.0e5, x=target_x, stoichiometry=stoichiometry),
            name=f"r{index}",
        )
        for index, stoichiometry in enumerate(reactions, start=1)
    ]

    result = mix.equilibrium(
        kind="reactive_electrolyte_bubble_pressure",
        T=313.15,
        P=1.0e5,
        balances={
            "water_total": {"H2O": 1.0},
            "amine_total": {"Amine": 1.0, "AmineH+": 1.0, "AmineCOO-": 1.0},
            "carbon_total": {"CO2": 1.0, "HCO3-": 1.0, "AmineCOO-": 1.0},
        },
        totals={
            "water_total": float(target_x[0]),
            "amine_total": float(target_x[1] + target_x[4] + target_x[5]),
            "carbon_total": float(target_x[2] + target_x[3] + target_x[5]),
        },
        reactions=reaction_definitions,
        initial_x=initial_x,
        vapor_species=["CO2", "H2O"],
        volatile_species=["CO2", "H2O"],
        nonvolatile_species=["Amine", "HCO3-", "AmineH+", "AmineCOO-", "H+"],
        options=epcsaft.ReactiveElectrolyteBubbleOptions(
            speciation_options=epcsaft.ReactiveSpeciationOptions(tolerance=1.0e-7, max_iterations=40, damping=0.7),
            bubble_options=epcsaft.ElectrolyteBubbleOptions(
                initial_pressure=1.0e5,
                tolerance=1.0e-6,
                return_best_effort=True,
            ),
            error_mode="result",
        ),
    )

    speciation_diagnostics = result.diagnostics["speciation"]["diagnostics"]
    bubble_diagnostics = result.diagnostics["bubble"]["diagnostics"]
    benchmark_report = {
        "reaction_residual_norm": max(abs(value) for value in result.reaction_residuals),
        "charge_residual_norm": abs(result.charge_residual),
        "material_residual_norm": max(abs(value) for value in result.mass_balance_residuals.values()),
        "co2_partial_pressure": result.partial_pressures["CO2"],
        "liquid_speciation": dict(result.x_liq),
        "activity_convention": tuple(speciation_diagnostics["reaction_standard_states"]),
        "derivative_backend": speciation_diagnostics["derivative_backend_by_block"]["reactive_speciation_variables"],
        "native_solver_iterations": {
            "speciation": speciation_diagnostics["iterations"],
            "bubble": bubble_diagnostics["iterations"],
        },
    }

    assert result.success is True
    assert benchmark_report["reaction_residual_norm"] <= 1.0e-7
    assert benchmark_report["charge_residual_norm"] <= 1.0e-8
    assert benchmark_report["material_residual_norm"] <= 1.0e-8
    assert benchmark_report["co2_partial_pressure"] > 0.0
    assert benchmark_report["liquid_speciation"]["CO2"] == pytest.approx(target_x[2], abs=5.0e-7)
    assert set(benchmark_report["activity_convention"]) == {"mole_fraction_activity"}
    assert benchmark_report["derivative_backend"] == "analytic_implicit"
    assert benchmark_report["native_solver_iterations"]["speciation"] > 0
    assert benchmark_report["native_solver_iterations"]["bubble"] > 0
    assert result.y_vap["CO2"] > 0.0
    assert result.y_vap["H2O"] > 0.0
    assert "HCO3-" not in result.y_vap
    assert "AmineH+" not in result.y_vap
    assert result.fugacity_residual_norm <= 1.0e-5
    assert speciation_diagnostics["solver_language"] == "c++"
    assert speciation_diagnostics["activity_model"] == "epcsaft_component_activity"
    assert speciation_diagnostics["activity_or_fugacity_terms_in_residual"] is True
    assert speciation_diagnostics["implicit_solve_results"]["reactive_speciation_variables"]["sensitivity"]
    assert bubble_diagnostics["native_entrypoint"] == "_solve_electrolyte_bubble_native"
    assert bubble_diagnostics["solver_language"] == "c++"
