from __future__ import annotations

import json

import pytest

import epcsaft
from tests.helpers.regression_cases import _methane_like_records, _minimal_neutral_metadata


def test_ceres_pure_neutral_regression_owns_optimizer_loop() -> None:
    ceres = epcsaft.runtime_build_info()["optional_dependencies"]["ceres"]
    if not ceres["compiled"]:
        pytest.skip("Ceres support is not enabled in this native build.")
    initial_guess = {"m": 1.08, "s": 3.55, "e": 155.0}

    result = epcsaft.fit_pure_neutral(
        _methane_like_records(),
        "Methane",
        assoc_scheme="",
        fixed_parameters=_minimal_neutral_metadata(16.043e-3),
        initial_guess=initial_guess,
        bounds={
            "m": (0.5, 3.5),
            "s": (2.0, 5.0),
            "e": (50.0, 400.0),
        },
        optimizer_backend="ceres",
    )

    assert result.success, result.message
    assert result.optimizer_backend == "ceres"
    assert result.backend == "ceres"
    assert result.derivative_backend == "cppad_implicit"
    assert result.jacobian_backend == "cppad_implicit"
    assert result.python_objective_used is False
    assert result.objective_final < result.objective_initial
    assert any(abs(result.fitted_values[name] - initial_guess[name]) > 1.0e-8 for name in ("m", "s", "e"))
    assert result.parameter_map == pytest.approx(result.fitted_values)
    assert result.row_diagnostics
    assert result.n_residual_evaluations >= 1
    assert result.n_jacobian_evaluations >= 1


def test_ceres_pure_ion_regression_uses_cppad_implicit_for_density_osmotic_miac(tmp_path) -> None:
    build = epcsaft.runtime_build_info()["optional_dependencies"]
    if not build["ceres"]["compiled"] or not build["cppad"]["compiled"]:
        pytest.skip("Ceres and CppAD support are not enabled in this native build.")

    dataset = tmp_path / "synthetic_ion_dataset"
    pure_dir = dataset / "pure"
    pure_dir.mkdir(parents=True)
    (pure_dir / "any_solvent.csv").write_text(
        "\n".join(
            [
                "component,MW,m,s,e,e_assoc,vol_a,assoc_scheme,z,dielc,d_born,f_solv",
                "Solv,0.018,1.0,3.7,150,0,0,,0,78,0,1",
                "Cat+,0.023,1,2.8,100,0,0,,1,8,3.4,1",
                "An-,0.035,1,2.7,100,0,0,,-1,8,4.1,1",
            ]
        ),
        encoding="utf-8",
    )
    (dataset / "user_options.json").write_text(
        json.dumps(
            {
                "elec_model": {
                    "rel_perm": {"rule": "constant"},
                    "include_born_model": True,
                    "born_model": {
                        "d_Born_mode": 3,
                        "solvation_shell_model": False,
                        "dielectric_saturation": False,
                        "mu_born_model": {"differential_mode": "auto", "comp_dep_delta_d": False},
                    },
                }
            }
        ),
        encoding="utf-8",
    )

    result = epcsaft.fit_pure_ion(
        [
            {
                "T": 298.15,
                "P": 101325.0,
                "x_Solv": 0.999998,
                "x_Cat+": 1.0e-6,
                "x_An-": 1.0e-6,
                "rho": 40.95,
                "osmotic_coefficient": 1.0,
                "mean_ionic_activity": 1.0,
            }
        ],
        "Cat+",
        dataset=dataset,
        species=["Solv", "Cat+", "An-"],
        solvent="Solv",
        fit_targets=["s", "e", "d_born"],
        initial_guess={"s": 2.8, "e": 100.0, "d_born": 3.4},
        bounds={"s": (2.4, 3.2), "e": (50.0, 300.0), "d_born": (2.0, 5.0)},
    )

    assert result.success, result.message
    assert result.backend == "ceres"
    assert result.jacobian_backend == "cppad_implicit"
    assert result.python_objective_used is False
    assert result.problem.fit_targets == ("s", "e", "d_born")
    assert set(result.metrics_by_term) == {"density", "osmotic_coefficient", "mean_ionic_activity"}
    assert result.provenance_report["parameter_movement"].keys() == {"s", "e", "d_born"}
