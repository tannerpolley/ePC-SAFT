from __future__ import annotations

import numpy as np

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
