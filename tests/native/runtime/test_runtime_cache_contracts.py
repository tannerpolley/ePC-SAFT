"""Native runtime contracts for pressure-vs-density closure and contribution terms."""

from __future__ import annotations

import json

import numpy as np
import pytest
from epcsaft._core import NativeValueError

from epcsaft import InputError, SolutionError, ePCSAFTMixture
from tests.helpers.native_cases import _ionic_state, _neutral_state

def _assert_close_terms(observed: dict[str, float], expected: dict[str, float]) -> None:
    assert set(observed) == set(expected)
    for key, value in expected.items():
        assert observed[key] == pytest.approx(value, rel=1e-10, abs=1e-12)

def _assert_finite_mapping_values(values: dict[str, float]) -> None:
    assert values
    assert all(np.isfinite(float(value)) for value in values.values())

def test_runtime_cache_stats_track_density_and_reference_state_reuse() -> None:
    mix, species, pressure, _, temperature, composition = _ionic_state()
    mix.clear_runtime_caches()
    mix.reset_runtime_cache_stats()

    first = mix.state(T=temperature, x=composition, P=pressure, phase="liq")
    second = mix.state(T=temperature, x=composition, P=pressure, phase="liq")

    assert first.density() == pytest.approx(second.density())
    stats = mix.runtime_cache_stats()
    assert stats["density_warm_start_hits"] >= 1
    assert stats["density_warm_start_fallbacks"] == 0

    for _ in range(3):
        second.activity_coefficient(species=species, mean_ionic_form=True, basis="molality")

    stats = mix.runtime_cache_stats()
    assert stats["reference_state_cache_misses"] == 1
    assert stats["reference_state_cache_hits"] >= 2

def test_activity_coefficient_cache_behavior_distinguishes_aux_cache_from_solvent_override() -> None:
    mix, species, pressure, _, temperature, composition = _ionic_state()
    state = mix.state(T=temperature, x=composition, P=pressure, phase="liq")
    mix.clear_runtime_caches()
    mix.reset_runtime_cache_stats()

    first_solvation = state.solvation_free_energy(species=species)
    after_first_solvation = mix.runtime_cache_stats()
    second_solvation = state.solvation_free_energy(species=species)
    after_second_solvation = mix.runtime_cache_stats()

    assert second_solvation == pytest.approx(first_solvation)
    assert after_second_solvation == after_first_solvation

    first_override = state.activity_coefficient(species=species, solvent="water")
    after_first_override = mix.runtime_cache_stats()
    second_override = state.activity_coefficient(species=species, solvent="water")
    after_second_override = mix.runtime_cache_stats()

    assert second_override == pytest.approx(first_override)
    assert after_first_override["reference_state_cache_hits"] > after_second_solvation["reference_state_cache_hits"]
    assert after_second_override["reference_state_cache_hits"] > after_first_override["reference_state_cache_hits"]
    assert after_second_override["reference_state_cache_misses"] == after_first_override["reference_state_cache_misses"]

def test_runtime_cache_stats_track_warm_start_fallbacks_without_hiding_failures() -> None:
    mix, _, pressure, _, temperature, composition = _ionic_state()
    mix.clear_runtime_caches()
    mix.reset_runtime_cache_stats()

    first = mix.state(T=temperature, x=composition, P=pressure, phase="liq")
    before_failure = mix.runtime_cache_stats()
    with pytest.raises(SolutionError):
        mix.state(T=temperature, x=composition, P=1.0e12, phase="liq")
    after_failure = mix.runtime_cache_stats()

    assert first.density() > 0.0
    assert before_failure == {
        "reference_state_cache_hits": 0,
        "reference_state_cache_misses": 0,
        "density_warm_start_hits": 0,
        "density_warm_start_fallbacks": 0,
    }
    assert after_failure["density_warm_start_fallbacks"] == 1

def test_miac_electrolyte_variants_cover_water_nonaqueous_and_mixed_solvents() -> None:
    from analyses.miac_fits.scripts import validate_miac_fits as vmf
    from scripts._epcsaft_oop import as_mixture

    def find_combo(solvent_system: str, salt: str, comp_signature=None) -> dict[str, object]:
        for combo in vmf.discover_combos(solvent_scope=solvent_system, salt_scope=salt):
            if comp_signature is None:
                if not combo.get("comp_signature"):
                    return combo
            elif combo.get("comp_signature") == comp_signature:
                return combo
        raise AssertionError(f"Missing MIAC combo for {solvent_system}/{salt} with {comp_signature!r}.")

    cases = [
        ("water", "NaCl", None, "Na+Cl-"),
        ("methanol", "KCl", None, "K+Cl-"),
        ("water-ethanol", "NaBr", (("water", 0.0), ("ethanol", 1.0)), "Na+Br-"),
    ]

    for solvent_system, salt, comp_signature, pair_key in cases:
        combo = find_combo(solvent_system, salt, comp_signature)
        species = vmf._species_for_combo(salt, solvent_system)
        params = vmf.build_params_for_variant(
            "2025_Figiel",
            combo,
            user_options=dict(vmf.DATASET_VARIANTS["2025_Figiel"].get("user_options", {})),
        )
        mix = as_mixture(params, species=species)
        composition = vmf._molality_to_molefraction_combo(0.5, salt, solvent_system, dict(combo.get("comp", {})))
        state = mix.state(T=vmf.T_REF, x=composition, P=vmf.P_REF, phase="liq")

        component_activity = state.activity_coefficient(species=species)
        mean_molality = state.activity_coefficient(species=species, mean_ionic_form=True, basis="molality")
        mean_mole = state.activity_coefficient(species=species, mean_ionic_form=True, basis="mole")
        solvation_free_energy = state.solvation_free_energy(species=species)
        osmotic_coefficient = state.osmotic_coefficient()
        dadt = state.temperature_derivative_residual_helmholtz(return_contribution_terms=True)
        dadx = state.composition_derivative_residual_helmholtz()

        _assert_finite_mapping_values(component_activity)
        _assert_finite_mapping_values(mean_molality)
        _assert_finite_mapping_values(mean_mole)
        _assert_finite_mapping_values(solvation_free_energy)
        assert np.all(np.isfinite(np.asarray(osmotic_coefficient, dtype=float)))
        assert set(dadt["terms"]) == {"hc", "disp", "assoc", "ion", "born"}
        assert set(dadx["terms"]) == {"hc", "disp", "assoc", "ion", "born"}
        assert pair_key in mean_molality
        assert pair_key in mean_mole
        assert state.pressure() == pytest.approx(vmf.P_REF)

def test_miac_activity_cache_reuse_keeps_results_stable() -> None:
    from analyses.miac_fits.scripts import validate_miac_fits as vmf
    from scripts._epcsaft_oop import as_mixture

    combo = next(
        combo
        for combo in vmf.discover_combos(solvent_scope="water-ethanol", salt_scope="NaBr")
        if combo.get("comp_signature") == (("water", 0.0), ("ethanol", 1.0))
    )
    species = vmf._species_for_combo("NaBr", "water-ethanol")
    params = vmf.build_params_for_variant(
        "2025_Figiel",
        combo,
        user_options=dict(vmf.DATASET_VARIANTS["2025_Figiel"].get("user_options", {})),
    )
    mix = as_mixture(params, species=species)
    composition = vmf._molality_to_molefraction_combo(0.5, "NaBr", "water-ethanol", dict(combo.get("comp", {})))

    mix.clear_runtime_caches()
    mix.reset_runtime_cache_stats()
    first = mix.state(T=vmf.T_REF, x=composition, P=vmf.P_REF, phase="liq")
    first_gamma = first.activity_coefficient(species=species, mean_ionic_form=True, basis="molality")
    second = mix.state(T=vmf.T_REF, x=composition, P=vmf.P_REF, phase="liq")
    second_gamma = second.activity_coefficient(species=species, mean_ionic_form=True, basis="molality")
    stats = mix.runtime_cache_stats()

    assert second.density() == pytest.approx(first.density())
    assert second_gamma == pytest.approx(first_gamma)
    assert stats["density_warm_start_hits"] >= 1
