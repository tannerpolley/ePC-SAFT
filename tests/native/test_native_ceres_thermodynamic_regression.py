from __future__ import annotations

import math

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


def _salt_speciation_ssmds_mixture() -> epcsaft.ePCSAFTMixture:
    params = {
        "m": np.asarray([1.2047, 1.0, 1.0, 1.0]),
        "s": np.asarray([2.7927, 3.0, 2.8232, 2.7560]),
        "e": np.asarray([353.95, 200.0, 230.0, 170.0]),
        "z": np.asarray([0.0, 0.0, 1.0, -1.0]),
        "dielc": np.asarray([78.09, 8.0, 8.0, 8.0]),
        "d_born": np.asarray([0.0, 0.0, 3.445, 4.1]),
        "f_solv": np.asarray([1.5, 1.0, 1.0, 1.0]),
        "MW": np.asarray([18.01528e-3, 58.44e-3, 22.989e-3, 35.45e-3]),
        "elec_model": {
            "include_born_model": True,
            "born_model": {
                "d_Born_mode": 3,
                "solvation_shell_model": True,
                "dielectric_saturation": True,
                "bulk_mode": "solvent",
                "mu_born_model": {
                    "comp_dep_rel_perm": True,
                    "include_sum_term": True,
                    "comp_dep_delta_d": True,
                },
            },
        },
    }
    return epcsaft.ePCSAFTMixture.from_params(params, species=["H2O", "NaCl", "Na+", "Cl-"])


def _salt_speciation_row(log_k: float, observed_na: float) -> dict[str, object]:
    initial_x = [0.998, 0.001, 0.0005, 0.0005]
    return {
        "row_id": "speciation_1",
        "row_mode": "reactive_speciation",
        "T": 298.15,
        "P": 1.0e5,
        "initial_x": initial_x,
        "balance_matrix": [
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            1.0,
            0.0,
            0.0,
            1.0,
            0.0,
            1.0,
        ],
        "balance_rows": 3,
        "total_vector": [0.998, 0.0015, 0.0015],
        "reaction_stoichiometry": [0.0, -1.0, 1.0, 1.0],
        "reaction_rows": 1,
        "log_equilibrium_constants": [log_k],
        "reaction_standard_states": [1],
        "options": {"jacobian_backend": "auto", "max_iterations": 50, "tolerance": 1.0e-8},
        "targets": [{"family": "speciation", "target": "Na+", "index": 2, "observed": observed_na, "scale": 1.0}],
    }


def test_native_thermo_regression_evaluates_reactive_speciation_row_in_cpp() -> None:
    species = ["H2O", "NaCl", "Na+", "Cl-"]
    mix = _salt_speciation_mixture()
    initial_x = [0.998, 0.001, 0.0005, 0.0005]
    log_k = math.log(initial_x[2]) + math.log(initial_x[3]) - math.log(initial_x[1])

    result = epcsaft.evaluate_native_thermo_regression_rows(
        mix,
        {
            "species": species,
            "rows": [
                {
                    "row_id": "speciation_1",
                    "row_mode": "reactive_speciation",
                    "T": 298.15,
                    "P": 1.0e5,
                    "initial_x": initial_x,
                    "balance_matrix": [
                        1.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        1.0,
                        1.0,
                        0.0,
                        0.0,
                        1.0,
                        0.0,
                        1.0,
                    ],
                    "balance_rows": 3,
                    "total_vector": [0.998, 0.0015, 0.0015],
                    "reaction_stoichiometry": [0.0, -1.0, 1.0, 1.0],
                    "reaction_rows": 1,
                    "log_equilibrium_constants": [log_k],
                    "reaction_standard_states": [1],
                    "options": {"jacobian_backend": "auto", "max_iterations": 50, "tolerance": 1.0e-8},
                    "targets": [
                        {"family": "speciation", "target": "Na+", "index": 2, "observed": initial_x[2], "scale": 1.0},
                        {
                            "family": "reaction",
                            "target": "salt_dissociation",
                            "index": 0,
                            "observed": 0.0,
                            "scale": 1.0,
                        },
                    ],
                }
            ],
        },
    )

    assert result["success_count"] == 1
    assert result["failure_count"] == 0
    assert result["fixed_shape_residuals"] is True
    assert result["row_diagnostics"][0]["solve_backend"] == "native_chemical_equilibrium"
    assert result["row_diagnostics"][0]["derivative_backend"] in {"analytic", "autodiff"}
    assert [entry["family"] for entry in result["residual_schema"]] == ["speciation", "reaction"]
    assert max(abs(value) for value in result["residuals"]) <= 1.0e-8


def test_native_thermo_regression_fit_reports_fixed_shape_result() -> None:
    species = ["H2O", "NaCl", "Na+", "Cl-"]
    mix = _salt_speciation_mixture()
    log_k = math.log(0.0005) + math.log(0.0005) - math.log(0.001)

    result = epcsaft.fit_native_thermo_regression(
        mix,
        {
            "species": species,
            "rows": [_salt_speciation_row(log_k, 0.00065)],
            "parameters": [
                {
                    "name": "salt.logK",
                    "kind": "reaction_equilibrium_constant",
                    "initial": log_k,
                    "lower": log_k - 5.0,
                    "upper": log_k + 5.0,
                    "metadata": {"row_id": "speciation_1", "reaction_index": "0"},
                }
            ],
            "options": {"max_iterations": 3, "derivative_backend": "implicit"},
        },
    )

    assert result["status"] in {"backend_unavailable", "converged", "max_iterations"}
    assert result["optimizer_backend"] in {"backend_unavailable", "ceres"}
    assert result["derivative_backend"] in {"implicit", "analytic_implicit", "cppad_implicit"}
    assert result["initial_cost"] >= 0.0
    assert result["objective_result"]["fixed_shape_residuals"] is True


def test_native_thermo_regression_reports_ssmds_born_derivatives_unavailable() -> None:
    species = ["H2O", "NaCl", "Na+", "Cl-"]
    mix = _salt_speciation_mixture()
    log_k = math.log(0.0005) + math.log(0.0005) - math.log(0.001)

    result = epcsaft.fit_native_thermo_regression(
        mix,
        {
            "species": species,
            "rows": [_salt_speciation_row(log_k, 0.00065)],
            "parameters": [
                {
                    "name": "Na+.d_born",
                    "kind": "born_radius",
                    "initial": 3.445,
                    "lower": 2.0,
                    "upper": 5.0,
                    "metadata": {"component_index": "2"},
                },
                {
                    "name": "H2O.f_solv",
                    "kind": "f_solv",
                    "initial": 1.5,
                    "lower": 0.5,
                    "upper": 3.0,
                    "metadata": {"component_index": "0"},
                },
            ],
            "options": {"max_iterations": 3, "derivative_backend": "implicit"},
        },
    )

    assert result["status"] == "backend_unavailable"
    assert result["optimizer_backend"] == "backend_unavailable"
    assert "activity-standard-state reactive_speciation rows only" in result["message"]
    assert result["objective_result"]["fixed_shape_residuals"] is True


def test_native_thermo_regression_supports_ssmds_born_radius_parameter_on_activity_rows() -> None:
    species = ["H2O", "NaCl", "Na+", "Cl-"]
    mix = _salt_speciation_ssmds_mixture()
    initial_x = np.asarray([0.998, 0.001, 0.0005, 0.0005], dtype=float)
    state = mix.state(T=298.15, P=1.0e5, x=initial_x, phase="liq")
    gamma = state.activity_coefficient(species=species)
    log_k = math.log(initial_x[2] * gamma["Na+"]) + math.log(initial_x[3] * gamma["Cl-"])
    log_k -= math.log(initial_x[1] * gamma["NaCl"])

    row = _salt_speciation_row(log_k, 0.00065)
    row["reaction_standard_states"] = [0]

    result = epcsaft.fit_native_thermo_regression(
        mix,
        {
            "species": species,
            "rows": [row],
            "parameters": [
                {
                    "name": "Na+.d_born",
                    "kind": "born_radius",
                    "initial": 3.30,
                    "lower": 2.0,
                    "upper": 5.0,
                    "metadata": {"component_index": "2"},
                }
            ],
            "options": {"max_iterations": 3, "derivative_backend": "implicit"},
        },
    )

    cppad_enabled = bool(epcsaft.runtime_build_info()["native_dependencies"]["cppad"]["enabled"])
    if not cppad_enabled:
        assert result["status"] == "backend_unavailable"
        assert "CppAD-enabled build" in result["message"]
        return

    assert result["status"] in {"converged", "max_iterations"}
    assert result["optimizer_backend"] == "ceres"
    assert result["derivative_backend"] in {"analytic_implicit", "cppad_implicit", "mixed_implicit"}


def test_native_thermo_regression_supports_ssmds_solvation_factor_parameter_on_activity_rows() -> None:
    species = ["H2O", "NaCl", "Na+", "Cl-"]
    mix = _salt_speciation_ssmds_mixture()
    initial_x = np.asarray([0.998, 0.001, 0.0005, 0.0005], dtype=float)
    state = mix.state(T=298.15, P=1.0e5, x=initial_x, phase="liq")
    gamma = state.activity_coefficient(species=species)
    log_k = math.log(initial_x[2] * gamma["Na+"]) + math.log(initial_x[3] * gamma["Cl-"])
    log_k -= math.log(initial_x[1] * gamma["NaCl"])

    row = _salt_speciation_row(log_k, 0.00065)
    row["reaction_standard_states"] = [0]

    result = epcsaft.fit_native_thermo_regression(
        mix,
        {
            "species": species,
            "rows": [row],
            "parameters": [
                {
                    "name": "H2O.f_solv",
                    "kind": "f_solv",
                    "initial": 1.35,
                    "lower": 0.5,
                    "upper": 3.0,
                    "metadata": {"component_index": "0"},
                }
            ],
            "options": {"max_iterations": 3, "derivative_backend": "implicit"},
        },
    )

    cppad_enabled = bool(epcsaft.runtime_build_info()["native_dependencies"]["cppad"]["enabled"])
    if not cppad_enabled:
        assert result["status"] == "backend_unavailable"
        assert "CppAD-enabled build" in result["message"]
        return

    assert result["status"] in {"converged", "max_iterations"}
    assert result["optimizer_backend"] == "ceres"
    assert result["derivative_backend"] in {"analytic_implicit", "cppad_implicit", "mixed_implicit"}


def test_native_thermo_regression_ideal_speciation_targets_are_invariant_to_born_inputs() -> None:
    species = ["H2O", "NaCl", "Na+", "Cl-"]
    initial_x = [0.998, 0.001, 0.0005, 0.0005]
    log_k = math.log(initial_x[2]) + math.log(initial_x[3]) - math.log(initial_x[1])

    baseline = epcsaft.evaluate_native_thermo_regression_rows(
        _salt_speciation_mixture(),
        {"species": species, "rows": [_salt_speciation_row(log_k, 0.00065)]},
    )

    shifted = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.2047, 1.0, 1.0, 1.0]),
            "s": np.asarray([2.7927, 3.0, 2.8232, 2.7560]),
            "e": np.asarray([353.95, 200.0, 230.0, 170.0]),
            "z": np.asarray([0.0, 0.0, 1.0, -1.0]),
            "dielc": np.asarray([78.09, 8.0, 8.0, 8.0]),
            "d_born": np.asarray([0.0, 0.0, 4.2, 4.8], dtype=float),
            "f_solv": np.asarray([2.25, 1.0, 1.0, 1.0], dtype=float),
            "MW": np.asarray([18.01528e-3, 58.44e-3, 22.989e-3, 35.45e-3]),
        },
        species=species,
    )
    moved = epcsaft.evaluate_native_thermo_regression_rows(
        shifted,
        {"species": species, "rows": [_salt_speciation_row(log_k, 0.00065)]},
    )

    assert baseline["row_diagnostics"][0]["derivative_backend"] in {"analytic", "autodiff"}
    assert moved["row_diagnostics"][0]["derivative_backend"] in {"analytic", "autodiff"}
    assert baseline["residuals"] == pytest.approx(moved["residuals"], abs=1.0e-12)
    assert baseline["cost"] == pytest.approx(moved["cost"], abs=1.0e-12)


def test_native_thermo_regression_row_can_request_autodiff_speciation_jacobian() -> None:
    species = ["H2O", "NaCl", "Na+", "Cl-"]
    mix = _salt_speciation_mixture()
    log_k = math.log(0.0005) + math.log(0.0005) - math.log(0.001)
    row = _salt_speciation_row(log_k, 0.00065)
    row["options"] = {"jacobian_backend": "autodiff", "max_iterations": 50, "tolerance": 1.0e-8}

    cppad_enabled = bool(epcsaft.runtime_build_info()["native_dependencies"]["cppad"]["enabled"])
    if not cppad_enabled:
        with pytest.raises(ValueError, match="CppAD-enabled build"):
            epcsaft.evaluate_native_thermo_regression_rows(mix, {"species": species, "rows": [row]})
        return

    result = epcsaft.evaluate_native_thermo_regression_rows(mix, {"species": species, "rows": [row]})
    assert result["row_diagnostics"][0]["derivative_backend"] == "autodiff"


def test_native_thermo_regression_penalizes_unsupported_row_mode() -> None:
    mix = _salt_speciation_mixture()

    result = epcsaft.evaluate_native_thermo_regression_rows(
        mix,
        {
            "species": ["H2O", "NaCl", "Na+", "Cl-"],
            "penalty_residual": 123.0,
            "rows": [
                {
                    "row_id": "bad",
                    "row_mode": "python_objective",
                    "T": 298.15,
                    "targets": [{"family": "speciation", "target": "Na+", "observed": 0.0}],
                }
            ],
        },
    )

    assert result["success_count"] == 0
    assert result["failure_count"] == 1
    assert result["residuals"] == pytest.approx([123.0])
    assert result["row_diagnostics"][0]["penalty_applied"] is True


def test_native_thermo_regression_supports_concentration_standard_state() -> None:
    species = ["H2O", "NaCl", "Na+", "Cl-"]
    mix = _salt_speciation_mixture()
    initial_x = np.asarray([0.998, 0.001, 0.0005, 0.0005], dtype=float)
    density = mix.state(T=298.15, P=1.0e5, x=initial_x, phase="liq").molar_density()
    log_k = math.log(density * initial_x[2]) + math.log(density * initial_x[3]) - math.log(density * initial_x[1])

    row = _salt_speciation_row(log_k, 0.00065)
    row["reaction_standard_states"] = [2]

    result = epcsaft.fit_native_thermo_regression(
        mix,
        {
            "species": species,
            "rows": [row],
            "parameters": [
                {
                    "name": "salt.logK",
                    "kind": "reaction_equilibrium_constant",
                    "initial": log_k,
                    "lower": log_k - 5.0,
                    "upper": log_k + 5.0,
                    "metadata": {"row_id": "speciation_1", "reaction_index": "0"},
                }
            ],
            "options": {"max_iterations": 3, "derivative_backend": "implicit"},
        },
    )

    cppad_enabled = bool(epcsaft.runtime_build_info()["native_dependencies"]["cppad"]["enabled"])
    if not cppad_enabled:
        assert result["status"] == "backend_unavailable"
        assert "CppAD-enabled build" in result["message"]
        return

    assert result["status"] in {"converged", "max_iterations"}
    assert result["optimizer_backend"] == "ceres"
    assert result["derivative_backend"] == "cppad_implicit"


def test_native_thermo_regression_supports_activity_standard_state() -> None:
    species = ["H2O", "NaCl", "Na+", "Cl-"]
    mix = _salt_speciation_mixture()
    initial_x = np.asarray([0.998, 0.001, 0.0005, 0.0005], dtype=float)
    state = mix.state(T=298.15, P=1.0e5, x=initial_x, phase="liq")
    gamma = state.activity_coefficient(species=species)
    log_k = math.log(initial_x[2] * gamma["Na+"]) + math.log(initial_x[3] * gamma["Cl-"])
    log_k -= math.log(initial_x[1] * gamma["NaCl"])

    row = _salt_speciation_row(log_k, 0.00065)
    row["reaction_standard_states"] = [0]

    result = epcsaft.fit_native_thermo_regression(
        mix,
        {
            "species": species,
            "rows": [row],
            "parameters": [
                {
                    "name": "salt.logK",
                    "kind": "reaction_equilibrium_constant",
                    "initial": log_k,
                    "lower": log_k - 5.0,
                    "upper": log_k + 5.0,
                    "metadata": {"row_id": "speciation_1", "reaction_index": "0"},
                }
            ],
            "options": {"max_iterations": 3, "derivative_backend": "implicit"},
        },
    )

    cppad_enabled = bool(epcsaft.runtime_build_info()["native_dependencies"]["cppad"]["enabled"])
    if not cppad_enabled:
        assert result["status"] == "backend_unavailable"
        assert "CppAD-enabled build" in result["message"]
        return

    assert result["status"] in {"converged", "max_iterations"}
    assert result["optimizer_backend"] == "ceres"
    assert result["derivative_backend"] == "cppad_implicit"


def test_native_thermo_regression_supports_single_vapor_bubble_pressure_slice() -> None:
    species = ["H2O", "NaCl", "Na+", "Cl-"]
    mix = _salt_speciation_ssmds_mixture()
    x_liq = [0.998, 0.001, 0.0005, 0.0005]
    bubble = mix.electrolyte_bubble_p(
        298.15,
        x_liq=x_liq,
        vapor_species=["H2O"],
    )

    result = epcsaft.fit_native_thermo_regression(
        mix,
        {
            "species": species,
            "rows": [
                {
                    "row_id": "bubble_1",
                    "row_mode": "reactive_electrolyte_bubble",
                    "T": 298.15,
                    "initial_x": x_liq,
                    "x_liq": x_liq,
                    "vapor_species": ["H2O"],
                    "initial_pressure": bubble.P,
                    "min_pressure": bubble.P * 0.3,
                    "max_pressure": bubble.P * 3.0,
                    "targets": [
                        {
                            "family": "pressure",
                            "target": "bubble_pressure",
                            "index": 0,
                            "observed": bubble.P,
                            "scale": 1.0e-5,
                        }
                    ],
                }
            ],
            "parameters": [
                {
                    "name": "H2O.f_solv",
                    "kind": "f_solv",
                    "initial": 1.5,
                    "lower": 0.5,
                    "upper": 3.0,
                    "metadata": {"component_index": "0"},
                }
            ],
            "options": {"max_iterations": 2, "derivative_backend": "implicit"},
        },
    )

    assert result["status"] == "converged"
    assert result["optimizer_backend"] == "ceres"
    assert result["derivative_backend"] == "cppad_implicit"
    assert result["objective_result"]["row_diagnostics"][0]["solve_backend"] == "native_electrolyte_bubble"
    assert result["objective_result"]["row_diagnostics"][0]["derivative_backend"] == "autodiff"
