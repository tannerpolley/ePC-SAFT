# -*- coding: utf-8 -*-
"""Native runtime contracts for pressure-vs-density closure and contribution terms."""

from __future__ import annotations

import numpy as np
import pytest

from epcsaft import SolutionError
from epcsaft import ePCSAFTMixture


def _neutral_state() -> tuple[object, list[str], float, float, float, np.ndarray]:
    species = ["A", "B", "C"]
    params = {
        "m": np.asarray([1.0000, 1.6069, 2.0020]),
        "s": np.asarray([3.7039, 3.5206, 3.6184]),
        "e": np.asarray([150.03, 191.42, 208.11]),
        "k_ij": np.asarray([
            [0.0, 3.0e-4, 1.15e-2],
            [3.0e-4, 0.0, 5.10e-3],
            [1.15e-2, 5.10e-3, 0.0],
        ]),
    }
    mix = ePCSAFTMixture.from_params(params, species=species)
    composition = np.array([0.1, 0.3, 0.6])
    temperature = 233.15
    density = 14330.417110
    pressure = 1276374.1152948933
    return mix, species, pressure, density, temperature, composition


def _ionic_state() -> tuple[object, list[str], float, float, float, np.ndarray]:
    temperature = 298.15
    species = ["water", "Na+", "Cl-"]
    s_water = 2.7927 + 10.11 * np.exp(-0.01775 * temperature) - 1.417 * np.exp(-0.01146 * temperature)
    params = {
        "MW": np.asarray([18.01528e-3, 22.98e-3, 35.45e-3]),
        "m": np.asarray([1.2047, 1.0, 1.0]),
        "s": np.asarray([s_water, 2.8232, 2.7560]),
        "e": np.asarray([353.95, 230.0, 170.0]),
        "e_assoc": np.asarray([2425.7, 0.0, 0.0]),
        "vol_a": np.asarray([0.04509, 0.0, 0.0]),
        "assoc_scheme": ["2B", None, None],
        "z": np.asarray([0.0, 1.0, -1.0]),
        "dielc": np.asarray([78.09, 8.0, 8.0]),
        "d_born": np.asarray([0.0, 3.445, 4.1]),
        "f_solv": np.asarray([1.5, 1.0, 1.0]),
        "k_ij": np.asarray([
            [0.0, 0.0045, -0.25],
            [0.0045, 0.0, 0.317],
            [-0.25, 0.317, 0.0],
        ]),
        "l_ij": np.zeros((3, 3)),
        "k_hb": np.zeros((3, 3)),
    }
    mix = ePCSAFTMixture.from_params(params, species=species)
    composition = np.array([0.9998, 1.0e-4, 1.0e-4])
    density = 55344.274540081075
    pressure = 1.0e5
    return mix, species, pressure, density, temperature, composition


def _assert_close_terms(observed: dict[str, float], expected: dict[str, float]) -> None:
    assert set(observed) == set(expected)
    for key, value in expected.items():
        assert observed[key] == pytest.approx(value, rel=0.0, abs=1e-12)


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


def test_native_residual_helmholtz_and_compressibility_contributions_match_neutral_contract() -> None:
    mix, _, _, density, temperature, composition = _neutral_state()
    state = mix.state(T=temperature, x=composition, rho=density)
    ares = state.ares(return_contribution_terms=True)
    z = state.z(return_contribution_terms=True)

    assert ares["total"] == pytest.approx(-3.54988545131505, rel=0.0, abs=1e-12)
    assert z["total"] == pytest.approx(0.04594621208078564, rel=0.0, abs=1e-12)
    _assert_close_terms(ares["terms"], {
        "hc": 3.774229851214634,
        "disp": -7.324115302529684,
        "assoc": 0.0,
        "ion": 0.0,
        "born": 0.0,
    })
    _assert_close_terms(z["terms"], {
        "hc": 7.122473867439451,
        "disp": -8.076527655358666,
        "assoc": 0.0,
        "ion": 0.0,
        "born": 0.0,
        "ideal": 1.0,
    })


def test_native_residual_helmholtz_and_compressibility_contributions_match_ionic_contract() -> None:
    mix, _, _, density, temperature, composition = _ionic_state()
    state = mix.state(T=temperature, x=composition, rho=density)
    ares = state.ares(return_contribution_terms=True)
    z = state.z(return_contribution_terms=True)

    assert ares["total"] == pytest.approx(-9.7214027218058, rel=0.0, abs=1e-12)
    assert z["total"] == pytest.approx(0.000728884077611683, rel=0.0, abs=1e-12)
    _assert_close_terms(ares["terms"], {
        "hc": 4.5498342977047095,
        "disp": -8.862194941025747,
        "assoc": -5.369357675632981,
        "ion": -1.1229434731248254e-05,
        "born": -0.03967317341704953,
    })
    _assert_close_terms(z["terms"], {
        "hc": 10.033753448769597,
        "disp": -7.956283347485374,
        "assoc": -3.0767358436684233,
        "ion": -5.373538203676085e-06,
        "born": 0.0,
        "ideal": 1.0,
    })


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


def test_composition_derivative_contribution_terms_are_accounted_for() -> None:
    for state_factory in (_neutral_state, _ionic_state):
        mix, _, _, density, temperature, composition = state_factory()
        state = mix.state(T=temperature, x=composition, rho=density)
        derivative = state.composition_derivative_residual_helmholtz()

        total_from_terms = sum(derivative["terms"].values())
        np.testing.assert_allclose(total_from_terms, derivative["total"])
        assert derivative["z_total"] == pytest.approx(1.0 + sum(derivative["z_terms"].values()))
        assert set(derivative["terms"]) == {"hc", "disp", "assoc", "ion", "born"}


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
