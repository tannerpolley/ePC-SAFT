from __future__ import annotations

import numpy as np

from epcsaft import ePCSAFTMixture


def _neutral_state():
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
    state = mix.state(T=233.15, x=np.array([0.1, 0.3, 0.6]), rho=14330.417110)
    return state, species


def _ionic_params():
    t = 298.15
    s_water = 2.7927 + 10.11 * np.exp(-0.01775 * t) - 1.417 * np.exp(-0.01146 * t)
    return {
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


def _ionic_state():
    t = 298.15
    species = ["water", "Na+", "Cl-"]
    params = _ionic_params()
    mix = ePCSAFTMixture.from_params(params, species=species)
    state = mix.state(T=t, x=np.array([0.9998, 1.0e-4, 1.0e-4]), P=1.0e5)
    return state, species


def _ionic_state_with_elec_model(elec_model):
    t = 298.15
    species = ["water", "Na+", "Cl-"]
    params = _ionic_params()
    params["elec_model"] = elec_model
    mix = ePCSAFTMixture.from_params(params, species=species)
    state = mix.state(T=t, x=np.array([0.9998, 1.0e-4, 1.0e-4]), P=1.0e5)
    return state, species


def _assert_array(actual, expected, rtol=1e-8, atol=1e-10):
    np.testing.assert_allclose(np.asarray(actual, dtype=float), np.asarray(expected, dtype=float), rtol=rtol, atol=atol)


def _sum_term_arrays(terms):
    total = None
    for value in terms.values():
        arr = np.asarray(value, dtype=float)
        total = arr.copy() if total is None else total + arr
    return total
