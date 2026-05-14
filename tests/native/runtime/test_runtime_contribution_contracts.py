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

def test_temperature_derivative_reports_finite_accounted_terms() -> None:
    mix, _, _, density, temperature, composition = _neutral_state()
    state = mix.state(T=temperature, x=composition, rho=density)

    derivative = state.temperature_derivative_residual_helmholtz(return_contribution_terms=True)

    assert np.isfinite(derivative["total"])
    assert sum(derivative["terms"].values()) == pytest.approx(derivative["total"])

def test_temperature_derivative_is_available_across_density_branches() -> None:
    mix, _, _, _, _, composition = _neutral_state()
    states = [
        mix.state(T=300.0, x=composition, P=1.0e3, phase="vap"),
        mix.state(T=300.0, x=composition, P=1.0e3, phase="liq"),
    ]

    for state in states:
        derivative = state.temperature_derivative_residual_helmholtz(return_contribution_terms=True)

        assert np.isfinite(derivative["total"])
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

def test_composition_derivative_reports_finite_accounted_terms() -> None:
    for state_factory in (_neutral_state, _ionic_state):
        mix, _, _, density, temperature, composition = state_factory()
        state = mix.state(T=temperature, x=composition, rho=density)
        derivative = np.asarray(state.composition_derivative_residual_helmholtz()["total"], dtype=float)

        assert derivative.shape == composition.shape
        assert np.all(np.isfinite(derivative))

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
