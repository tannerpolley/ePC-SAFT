from __future__ import annotations

import json

import pytest

import epcsaft


def _write_ssmds_dataset(tmp_path):
    dataset = tmp_path / "synthetic_liquid_electrolyte"
    pure_dir = dataset / "pure"
    pure_dir.mkdir(parents=True)
    (pure_dir / "components.csv").write_text(
        "\n".join(
            [
                "component,MW,m,s,e,e_assoc,vol_a,assoc_scheme,z,dielc,d_born,f_solv",
                "Solv,0.018,1.0,3.7,150,0,0,,0,78,0,1.4",
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
                    "rel_perm": {"rule": "linear", "differential_mode": "auto"},
                    "include_born_model": True,
                    "born_model": {
                        "d_Born_mode": 3,
                        "solvation_shell_model": True,
                        "dielectric_saturation": True,
                        "bulk_mode": "mix",
                        "mu_born_model": {"differential_mode": "auto", "comp_dep_delta_d": True},
                    },
                }
            }
        ),
        encoding="utf-8",
    )
    return dataset


def test_ceres_liquid_electrolyte_regression_uses_native_ssmds_derivatives(tmp_path) -> None:
    build = epcsaft.runtime_build_info()["optional_dependencies"]
    if not build["ceres"]["compiled"] or not build["cppad"]["compiled"]:
        pytest.skip("Ceres and CppAD support are not enabled in this native build.")

    dataset = _write_ssmds_dataset(tmp_path)
    species = ("Solv", "Cat+", "An-")
    x = [0.98, 0.01, 0.01]
    temperature = 298.15
    pressure = 101325.0
    mixture = epcsaft.ePCSAFTMixture.from_dataset(dataset, species, x, temperature)
    state = mixture.state(T=temperature, x=x, P=pressure, phase="liq")
    epsilon_r, _ = state.relative_permittivity()
    miac = state.activity_coefficient(mean_ionic_form=True, basis="mole")["Cat+An-"]

    result = epcsaft.fit_liquid_electrolyte_parameters(
        species=species,
        data_rows=[
            {
                "T": temperature,
                "P": pressure,
                "x_Solv": x[0],
                "x_Cat+": x[1],
                "x_An-": x[2],
                "rho": state.molar_density(),
                "epsilon_r_exp": epsilon_r,
                "osmotic_coefficient": float(state.osmotic_coefficient()[0]),
                "mean_ionic_activity": miac,
            }
        ],
        parameter_set=dataset,
        parameters_to_fit=("d_born", "f_solv"),
        initial_guess={"d_born": 2.0, "f_solv": 0.5},
        bounds={"d_born": (1.0, 8.0), "f_solv": (0.1, 3.0)},
        solver_options={"optimizer_backend": "ceres"},
    )

    assert result.success, result.message
    assert result.backend == "ceres"
    assert result.optimizer_backend == "ceres"
    assert result.derivative_backend == "cppad_implicit"
    assert result.jacobian_backend == "cppad_implicit"
    assert result.jacobian_fallback_used is False
    assert result.python_objective_used is False
    assert result.objective_final < result.objective_initial
    assert result.problem.mode == "liquid_electrolyte"
    assert result.problem.fit_targets == ("d_born", "f_solv")
    assert result.problem.solver_options["target_components"] == {"d_born": "Cat+", "f_solv": "Solv"}
    assert result.row_diagnostics
    assert {row["row_family"] for row in result.row_diagnostics} == {
        "density",
        "relative_permittivity",
        "osmotic_coefficient",
        "mean_ionic_activity",
    }
    assert result.parameter_movement.keys() == {"d_born", "f_solv"}
    assert any(abs(delta) > 1.0e-8 for delta in result.parameter_movement.values())
