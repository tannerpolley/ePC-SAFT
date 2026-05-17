from __future__ import annotations

import numpy as np

from epcsaft.epcsaft import _fit_generic_native_ceres


def test_ceres_binary_kij_regression_uses_native_cppad_implicit_jacobian() -> None:
    payload = {
        "MW": np.asarray([0.016043, 0.03007], dtype=float),
        "m": np.asarray([1.0, 1.6069], dtype=float),
        "s": np.asarray([3.7039, 3.5206], dtype=float),
        "e": np.asarray([150.03, 191.42], dtype=float),
        "e_assoc": np.asarray([0.0, 0.0], dtype=float),
        "vol_a": np.asarray([0.0, 0.0], dtype=float),
        "assoc_scheme": [None, None],
        "k_ij": np.asarray([[0.0, 0.01], [0.01, 0.0]], dtype=float),
        "l_ij": np.zeros((2, 2), dtype=float),
        "k_hb": np.zeros((2, 2), dtype=float),
    }
    records = [
        {
            "term_name": "binary_vle_fugacity_balance",
            "term": 5,
            "T": 180.0,
            "P": 1.0e6,
            "x": [0.35, 0.65],
            "y": [0.74, 0.26],
            "target_index": 0,
            "target_index_2": 1,
            "scale": 1.0,
        },
        {
            "term_name": "binary_vle_fugacity_balance",
            "term": 5,
            "T": 190.0,
            "P": 1.2e6,
            "x": [0.42, 0.58],
            "y": [0.79, 0.21],
            "target_index": 0,
            "target_index_2": 1,
            "scale": 1.0,
        },
    ]

    result = _fit_generic_native_ceres(
        [payload, payload],
        records,
        np.asarray([6], dtype=int),
        np.asarray([0], dtype=int),
        np.asarray([1], dtype=int),
        np.asarray([0.01], dtype=float),
        np.asarray([-0.15], dtype=float),
        np.asarray([0.10], dtype=float),
    )

    assert result["success"] is True
    assert result["backend"] == "ceres"
    assert result["optimizer_backend"] == "ceres"
    assert result["derivative_backend"] == "cppad_implicit"
    assert result["jacobian_backend"] == "cppad_implicit"
    assert result["jacobian_available"] is True
    assert result["cost"] <= result["initial_cost"] + 1.0e-12
    assert -0.15 <= float(result["x"][0]) <= 0.10
    assert "binary_vle_fugacity_balance" in result["metrics_by_term"]


def test_ceres_binary_kij_regression_accepts_associating_neutral_rows() -> None:
    payload = {
        "MW": np.asarray([0.01801528, 0.03007], dtype=float),
        "m": np.asarray([1.2047, 1.6069], dtype=float),
        "s": np.asarray([3.0, 3.5206], dtype=float),
        "e": np.asarray([353.95, 191.42], dtype=float),
        "e_assoc": np.asarray([2425.7, 0.0], dtype=float),
        "vol_a": np.asarray([0.04509, 0.0], dtype=float),
        "assoc_scheme": ["2B", None],
        "k_ij": np.asarray([[0.0, 0.01], [0.01, 0.0]], dtype=float),
    }
    records = [
        {
            "term_name": "binary_vle_fugacity_balance",
            "term": 5,
            "T": 298.15,
            "P": 100000.0,
            "x": [0.55, 0.45],
            "y": [0.2, 0.8],
            "target_index": 0,
            "target_index_2": 1,
            "scale": 1.0,
        }
    ]

    result = _fit_generic_native_ceres(
        [payload],
        records,
        np.asarray([6], dtype=int),
        np.asarray([0], dtype=int),
        np.asarray([1], dtype=int),
        np.asarray([0.01], dtype=float),
        np.asarray([-0.15], dtype=float),
        np.asarray([0.10], dtype=float),
        max_nfev=1,
    )

    assert result["success"] is True
    assert result["backend"] == "ceres"
    assert result["optimizer_backend"] == "ceres"
    assert result["derivative_backend"] == "cppad_implicit"
    assert result["jacobian_backend"] == "cppad_implicit"
    assert result["jacobian_available"] is True
    assert "binary_vle_fugacity_balance" in result["metrics_by_term"]
