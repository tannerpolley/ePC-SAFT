# -*- coding: utf-8 -*-
"""Native runtime contracts for pressure-vs-density closure and contribution terms."""

from __future__ import annotations

import numpy as np
import pytest

from epcsaft import InputError
from epcsaft import SolutionError
from epcsaft import ePCSAFTMixture
from epcsaft._core import NativeValueError
from tests.helpers.native_cases import _ionic_state
from tests.helpers.native_cases import _neutral_state


def _assert_close_terms(observed: dict[str, float], expected: dict[str, float]) -> None:
    assert set(observed) == set(expected)
    for key, value in expected.items():
        assert observed[key] == pytest.approx(value, rel=1e-10, abs=1e-12)


def _assert_finite_mapping_values(values: dict[str, float]) -> None:
    assert values
    assert all(np.isfinite(float(value)) for value in values.values())


def test_pressure_based_and_density_based_states_match_for_neutral_system() -> None:
    mix, _, pressure, density, temperature, composition = _neutral_state()
    from_rho = mix.state(T=temperature, x=composition, rho=density)
    from_p = mix.state(T=temperature, x=composition, P=pressure, phase="liq")

    assert from_p.density() == pytest.approx(from_rho.density())
    assert from_p.pressure() == pytest.approx(from_rho.pressure())
    assert from_p.z() == pytest.approx(from_rho.z())
    assert from_p.ares() == pytest.approx(from_rho.ares())


def test_pressure_based_and_density_based_states_match_for_ionic_system() -> None:
    mix, _, pressure, density, temperature, composition = _ionic_state()
    from_rho = mix.state(T=temperature, x=composition, rho=density)
    from_p = mix.state(T=temperature, x=composition, P=pressure, phase="liq")

    assert from_p.density() == pytest.approx(from_rho.density())
    assert from_p.pressure() == pytest.approx(from_rho.pressure())
    assert from_p.z() == pytest.approx(from_rho.z())
    assert from_p.ares() == pytest.approx(from_rho.ares())


def test_pressure_density_edge_cases_cover_vapor_and_liquid_extremes() -> None:
    mix, _, _, _, _, composition = _neutral_state()

    vapor = mix.state(T=600.0, x=composition, P=1.0, phase="vap")
    liquid = mix.state(T=220.0, x=composition, P=5.0e7, phase="liq")

    assert vapor.phase == 1
    assert liquid.phase == 0
    assert vapor.density() == pytest.approx(2.0045400150430712e-4, rel=1e-10)
    assert liquid.density() == pytest.approx(16076.977238412512, rel=1e-10)
    assert vapor.pressure() == pytest.approx(1.0)
    assert liquid.pressure() == pytest.approx(5.0e7)
    assert np.isfinite(vapor.z())
    assert np.isfinite(liquid.z())


def test_pressure_density_phase_branches_do_not_cross_at_two_phase_like_state() -> None:
    mix, _, _, _, _, composition = _neutral_state()
    temperature = 300.0
    pressure = 1.0e3

    vapor = mix.state(T=temperature, x=composition, P=pressure, phase="vap")
    liquid = mix.state(T=temperature, x=composition, P=pressure, phase="liq")

    assert vapor.phase == 1
    assert liquid.phase == 0
    assert vapor.pressure() == pytest.approx(pressure)
    assert liquid.pressure() == pytest.approx(pressure)
    assert vapor.density() == pytest.approx(0.4009505832238275, rel=1e-10)
    assert liquid.density() == pytest.approx(10700.137898056397, rel=1e-10)
    assert liquid.density() / vapor.density() > 1.0e4


def test_ionic_high_pressure_liquid_density_branch_remains_stable() -> None:
    mix, species, _, _, _, composition = _ionic_state()

    state = mix.state(T=320.0, x=composition, P=5.0e7, phase="liq")
    gamma = state.activity_coefficient(species=species, mean_ionic_form=True, basis="molality")

    assert state.phase == 0
    assert state.pressure() == pytest.approx(5.0e7)
    assert state.density() == pytest.approx(55177.833730914346, rel=1e-10)
    assert state.z() == pytest.approx(0.3405817367443399, rel=1e-10)
    assert gamma == {"Na+Cl-": pytest.approx(0.9295936238832602, rel=1e-10)}


def test_pressure_density_failure_reports_state_context_and_native_outcome() -> None:
    mix, _, _, _, temperature, composition = _ionic_state()

    with pytest.raises(SolutionError) as excinfo:
        mix.state(T=temperature, x=composition, P=1.0e-12, phase="liq")

    message = str(excinfo.value)
    assert "pressure-based state solve failed" in message
    assert "T=298.15" in message
    assert "P=1e-12" in message
    assert "phase=liq" in message
    assert "ncomp=3" in message
    assert "No valid density root found for liquid phase" in message


def test_pressure_density_invalid_inputs_fail_before_native_density_search() -> None:
    mix, _, pressure, _, temperature, composition = _neutral_state()

    with pytest.raises(InputError, match="State composition length"):
        mix.state(T=temperature, x=composition[:-1], P=pressure, phase="liq")
    with pytest.raises(InputError, match="phase must be"):
        mix.state(T=temperature, x=composition, P=pressure, phase="solid")
    with pytest.raises(InputError, match="Provide exactly one of P or rho"):
        mix.state(T=temperature, x=composition, P=pressure, rho=1.0, phase="liq")


def test_nonionic_state_rejects_electrolyte_only_activity_methods() -> None:
    mix, species, _, density, temperature, composition = _neutral_state()
    state = mix.state(T=temperature, x=composition, rho=density)

    with pytest.raises(NativeValueError, match="requires ionic species"):
        state.activity_coefficient(species=species)
    with pytest.raises(NativeValueError, match="requires ionic species"):
        state.solvation_free_energy(species=species)
    with pytest.raises(NativeValueError, match="requires ionic species"):
        state.osmotic_coefficient()


def test_native_residual_helmholtz_and_compressibility_contributions_match_neutral_contract() -> None:
    mix, _, _, density, temperature, composition = _neutral_state()
    state = mix.state(T=temperature, x=composition, rho=density)
    ares = state.ares(return_contribution_terms=True)
    z = state.z(return_contribution_terms=True)

    assert ares["total"] == pytest.approx(-3.54988545131505, rel=0.0, abs=1e-12)
    assert z["total"] == pytest.approx(0.04594621208078564, rel=0.0, abs=1e-12)
    _assert_close_terms(
        ares["terms"],
        {
            "hc": 3.774229851214634,
            "disp": -7.324115302529684,
            "assoc": 0.0,
            "ion": 0.0,
            "born": 0.0,
        },
    )
    _assert_close_terms(
        z["terms"],
        {
            "hc": 7.122473867439451,
            "disp": -8.076527655358666,
            "assoc": 0.0,
            "ion": 0.0,
            "born": 0.0,
            "ideal": 1.0,
        },
    )


def test_native_residual_helmholtz_and_compressibility_contributions_match_ionic_contract() -> None:
    mix, _, _, density, temperature, composition = _ionic_state()
    state = mix.state(T=temperature, x=composition, rho=density)
    ares = state.ares(return_contribution_terms=True)
    z = state.z(return_contribution_terms=True)

    assert ares["total"] == pytest.approx(-9.7214027218058, rel=0.0, abs=1e-12)
    assert z["total"] == pytest.approx(0.000728884077611683, rel=0.0, abs=1e-12)
    _assert_close_terms(
        ares["terms"],
        {
            "hc": 4.5498342977047095,
            "disp": -8.862194941025747,
            "assoc": -5.369357675632981,
            "ion": -1.1229434731248254e-05,
            "born": -0.03967317341704953,
        },
    )
    _assert_close_terms(
        z["terms"],
        {
            "hc": 10.033753448769597,
            "disp": -7.956283347485374,
            "assoc": -3.0767358436684233,
            "ion": -5.373538203676085e-06,
            "born": 0.0,
            "ideal": 1.0,
        },
    )


def test_temperature_derivative_matches_neutral_finite_difference() -> None:
    mix, _, _, density, temperature, composition = _neutral_state()
    state = mix.state(T=temperature, x=composition, rho=density)
    delta_t = 1.0e-3

    plus = mix.state(T=temperature + delta_t, x=composition, rho=density)
    minus = mix.state(T=temperature - delta_t, x=composition, rho=density)
    finite_difference = (plus.ares() - minus.ares()) / (2.0 * delta_t)
    derivative = state.temperature_derivative_residual_helmholtz(return_contribution_terms=True)

    assert derivative["total"] == pytest.approx(finite_difference, rel=1e-9, abs=1e-11)
    assert sum(derivative["terms"].values()) == pytest.approx(derivative["total"])


def test_temperature_derivative_matches_neutral_finite_difference_across_density_branches() -> None:
    mix, _, _, _, _, composition = _neutral_state()
    states = [
        mix.state(T=300.0, x=composition, P=1.0e3, phase="vap"),
        mix.state(T=300.0, x=composition, P=1.0e3, phase="liq"),
    ]
    delta_t = 1.0e-3

    for state in states:
        density = state.density()
        plus = mix.state(T=state.T + delta_t, x=composition, rho=density, phase="liq")
        minus = mix.state(T=state.T - delta_t, x=composition, rho=density, phase="liq")
        finite_difference = (plus.ares() - minus.ares()) / (2.0 * delta_t)
        derivative = state.temperature_derivative_residual_helmholtz(return_contribution_terms=True)

        assert derivative["total"] == pytest.approx(finite_difference, rel=1e-8, abs=1e-10)
        assert sum(derivative["terms"].values()) == pytest.approx(derivative["total"])


def test_composition_derivative_contribution_terms_are_accounted_for() -> None:
    for state_factory in (_neutral_state, _ionic_state):
        mix, _, _, density, temperature, composition = state_factory()
        state = mix.state(T=temperature, x=composition, rho=density)
        derivative = state.composition_derivative_residual_helmholtz()

        total_from_terms = sum(derivative["terms"].values())
        np.testing.assert_allclose(total_from_terms, derivative["total"])
        assert derivative["z_total"] == pytest.approx(1.0 + sum(derivative["z_terms"].values()))
        assert set(derivative["terms"]) == {"hc", "disp", "assoc", "ion", "born"}


def test_composition_derivative_matches_constrained_composition_finite_difference() -> None:
    for state_factory in (_neutral_state, _ionic_state):
        mix, _, _, density, temperature, composition = state_factory()
        state = mix.state(T=temperature, x=composition, rho=density)
        derivative = np.asarray(state.composition_derivative_residual_helmholtz()["total"], dtype=float)

        for i, j in ((0, 1), (1, 2), (0, 2)):
            delta_x = min(1.0e-6, 0.25 * float(composition[i]), 0.25 * float(composition[j]))
            plus = composition.copy()
            minus = composition.copy()
            plus[i] += delta_x
            plus[j] -= delta_x
            minus[i] -= delta_x
            minus[j] += delta_x
            finite_difference = (
                mix.state(T=temperature, x=plus, rho=density).ares()
                - mix.state(T=temperature, x=minus, rho=density).ares()
            ) / (2.0 * delta_x)

            assert derivative[i] - derivative[j] == pytest.approx(finite_difference, rel=1e-7, abs=1e-8)


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


def test_public_methods_expose_eqid_owned_contribution_groups() -> None:
    mix, _, _, density, temperature, composition = _ionic_state()
    state = mix.state(T=temperature, x=composition, rho=density)
    contribution_families = {"hc", "disp", "assoc", "ion", "born"}

    ares = state.ares(return_contribution_terms=True)
    mures = state.residual_chemical_potential(return_contribution_terms=True)
    fugcoef = state.fugacity_coefficient(return_contribution_terms=True)

    assert set(ares["terms"]) == contribution_families
    assert set(mures["terms"]) == contribution_families
    assert set(fugcoef["terms"]) == contribution_families
    np.testing.assert_allclose(sum(mures["terms"].values()), mures["total"])
    np.testing.assert_allclose(fugcoef["terms_total_natural_log"], sum(fugcoef["terms"].values()))


def test_miac_electrolyte_variants_cover_water_nonaqueous_and_mixed_solvents() -> None:
    from scripts.fits import validate_miac_fits as vmf
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
    from scripts.fits import validate_miac_fits as vmf
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
