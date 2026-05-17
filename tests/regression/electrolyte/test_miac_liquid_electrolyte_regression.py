from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

import epcsaft
from analyses.data_validation.miac_fits.scripts import validate_miac_fits as vmf
from scripts._epcsaft_oop import as_mixture

REPO_ROOT = Path(__file__).resolve().parents[3]
FIXTURE = REPO_ROOT / "tests" / "fixtures" / "literature" / "figiel_2025" / "miac_liquid_electrolyte.json"


def _load_fixture() -> dict:
    with FIXTURE.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def test_miac_liquid_electrolyte_fixture_records_figiel_provenance() -> None:
    fixture = _load_fixture()

    assert fixture["family"] == "liquid electrolyte d_born/f_solv checks"
    assert fixture["source"]["paper"] == "Figiel, Yu, Held 2025"
    for path in fixture["source"]["local_assets"]:
        assert (REPO_ROOT / path).exists()
    assert fixture["phase_scope"] == "liquid-electrolyte only"
    assert fixture["excluded_scope"] == "vapor Born derivatives"


def test_miac_liquid_electrolyte_regression_uses_native_ceres_without_nonexact_derivatives() -> None:
    fixture = _load_fixture()
    result = epcsaft.fit_liquid_electrolyte_parameters(
        species=tuple(fixture["regression_probe"]["species"]),
        data_rows=fixture["regression_probe"]["data_rows"],
        parameter_set=fixture["regression_probe"]["parameter_set"],
        parameters_to_fit=tuple(fixture["regression_probe"]["parameters_to_fit"]),
        initial_guess=fixture["regression_probe"]["initial_guess"],
        bounds={key: tuple(value) for key, value in fixture["regression_probe"]["bounds"].items()},
        solver_options={"optimizer_backend": "ceres"},
    )

    assert result.success is True
    assert result.backend == "ceres"
    assert result.optimizer_backend == "ceres"
    assert result.derivative_backend == "cppad_implicit"
    assert result.jacobian_backend == "cppad_implicit"
    assert result.python_objective_used is False
    assert result.objective_final <= result.objective_initial
    assert result.parameter_movement.keys() == {"d_born", "f_solv"}
    assert any(abs(value) > 1.0e-8 for value in result.parameter_movement.values())
    assert {row["row_family"] for row in result.row_diagnostics} == {
        "density",
        "relative_permittivity",
        "osmotic_coefficient",
        "mean_ionic_activity",
    }


def test_miac_liquid_electrolyte_figiel_outputs_and_coverage_matrix_are_finite() -> None:
    fixture = _load_fixture()
    case = fixture["miac_cases"][0]
    combo = next(
        combo
        for combo in vmf.discover_combos(solvent_scope=case["solvent_system"], salt_scope=case["salt"])
        if combo.get("comp_signature") == tuple(tuple(item) for item in case["comp_signature"])
    )
    species = vmf._species_for_combo(case["salt"], case["solvent_system"])
    params = vmf.build_params_for_variant("2025_Figiel", combo, user_options={})
    mixture = as_mixture(params, species=species)
    x = vmf._molality_to_molefraction_combo(
        case["molality"],
        case["salt"],
        case["solvent_system"],
        dict(combo.get("comp", {})),
    )
    state = mixture.state(T=vmf.T_REF, x=x, P=vmf.P_REF, phase="liq")

    mean_molality = state.activity_coefficient(species=species, mean_ionic_form=True, basis="molality")
    derivatives = state.born_ssmds_liquid_derivatives()
    coverage_rows = state.derivative_coverage_matrix()
    backend_labels = {str(row["backend"]).lower() for row in coverage_rows}

    assert case["pair_key"] in mean_molality
    assert mean_molality[case["pair_key"]] == pytest.approx(float(case["expected_miac"]), rel=1.0e-11, abs=1.0e-12)
    assert derivatives["supported"] is True
    assert derivatives["backend"] in {"analytic", "cppad"}
    assert np.all(np.isfinite(derivatives["a_born_d_d_born"]))
    assert np.all(np.isfinite(derivatives["a_born_d_f_solv"]))
    assert backend_labels <= {"analytic", "cppad", "analytic_implicit", "cppad_implicit", "out_of_scope"}
