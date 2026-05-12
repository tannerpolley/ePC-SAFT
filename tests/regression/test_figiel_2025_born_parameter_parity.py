from __future__ import annotations

import numpy as np
import pytest

import epcsaft
import epcsaft.regression as regression_module
from analyses.miac_fits.scripts import validate_miac_fits as vmf


def _figiel_nabr_water_combo() -> dict[str, object]:
    return next(
        combo
        for combo in vmf.discover_combos(solvent_scope="water", salt_scope="NaBr")
        if not combo.get("comp_signature")
    )


def _figiel_user_options() -> dict[str, object]:
    return dict(vmf.DATASET_VARIANTS["2025_Figiel"].get("user_options", {}))


def _figiel_species() -> list[str]:
    return vmf._species_for_combo("NaBr", "water")


def _figiel_params_and_mixture() -> tuple[list[str], dict[str, object], epcsaft.ePCSAFTMixture]:
    combo = _figiel_nabr_water_combo()
    species = _figiel_species()
    user_options = _figiel_user_options()
    params = vmf.build_params_for_variant("2025_Figiel", combo, user_options=user_options)
    return species, user_options, epcsaft.ePCSAFTMixture.from_params(params, species=species)


def _synthetic_figiel_miac_records() -> tuple[list[str], dict[str, object], list[dict[str, float]]]:
    combo = _figiel_nabr_water_combo()
    species, user_options, mixture = _figiel_params_and_mixture()
    records: list[dict[str, float]] = []
    for molality in (0.05, 0.1, 0.2, 0.4, 0.8):
        x = vmf._molality_to_molefraction_combo(molality, "NaBr", "water", dict(combo.get("comp", {})))
        state = mixture.state(T=vmf.T_REF, x=x, P=vmf.P_REF, phase="liq")
        mean_ionic_activity = state.activity_coefficient(
            species=species,
            mean_ionic_form=True,
            basis="molality",
        )["Na+Br-"]
        record: dict[str, float] = {
            "T": vmf.T_REF,
            "P": vmf.P_REF,
            "mean_ionic_activity": float(mean_ionic_activity),
            "pair_label": "Na+Br-",
        }
        for label, value in zip(species, x):
            record[f"x_{label}"] = float(value)
        records.append(record)
    return species, user_options, records


def _native_generic_fsolv_fit() -> dict[str, object]:
    species, user_options, records = _synthetic_figiel_miac_records()
    solvent = species[2]
    terms = regression_module._build_pure_ion_terms(regression_module._normalize_records(records))
    fixed_payloads: list[dict[str, object]] = []
    native_records: list[dict[str, object]] = []

    for term in terms:
        scale = regression_module._family_scale(term)
        for record in term.records:
            temperature = regression_module._float_from_record(record, "T", required=True)
            x = regression_module._ion_composition_from_record(record, species, solvent)
            miac_i, miac_j = regression_module._native_miac_pair_indices(record, species)
            native_records.append(
                {
                    "term_name": term.term_type,
                    "term": regression_module.NATIVE_TERM_KINDS[term.term_type],
                    "T": temperature,
                    "P": regression_module._float_from_record(record, "P", required=True),
                    "phase": regression_module.phase_to_int("liq"),
                    "x": x.tolist(),
                    "target": regression_module._float_from_record(
                        record,
                        "mean_ionic_activity",
                        "mean_ionic_activity_coefficient",
                        "miac",
                        required=True,
                    ),
                    "target_index": miac_i,
                    "target_index_2": miac_j,
                    "activity_basis": 1,
                    "solvent_index": species.index(solvent),
                    "scale": scale,
                }
            )
            fixed_payloads.append(
                regression_module._params_for_native_record(
                    "2025_Figiel",
                    species,
                    x,
                    temperature,
                    user_options,
                )
            )

    return regression_module._run_native_generic_least_squares(
        fixed_payloads,
        native_records,
        ["solvation_factor"],
        species,
        np.asarray([1.2], dtype=float),
        np.asarray([0.5], dtype=float),
        np.asarray([3.0], dtype=float),
        component=solvent,
        max_nfev=200,
    )


def test_fit_pure_ion_recovers_figiel_2025_na_born_radius_from_synthetic_miac() -> None:
    species, user_options, records = _synthetic_figiel_miac_records()
    result = epcsaft.fit_pure_ion(
        records,
        "Na+",
        dataset="2025_Figiel",
        species=species,
        solvent=species[2],
        fit_targets=("d_born",),
        initial_guess={"d_born": 3.1},
        bounds={"d_born": (2.0, 5.0)},
        user_options=user_options,
    )

    assert result.success, result.message
    assert result.backend == "least_squares_native"
    assert result.fitted_values["d_born"] == pytest.approx(3.445, abs=1.0e-6)
    assert result.metrics_by_term["mean_ionic_activity"] <= 1.0e-8


def test_native_generic_liquid_electrolyte_fit_recovers_figiel_2025_water_f_solv() -> None:
    result = _native_generic_fsolv_fit()

    assert result["success"] is True
    assert result["backend"] == "least_squares_native"
    assert float(result["x"][0]) == pytest.approx(1.5, abs=1.0e-6)
    assert result["metrics_by_term"]["mean_ionic_activity"] <= 1.0e-8
