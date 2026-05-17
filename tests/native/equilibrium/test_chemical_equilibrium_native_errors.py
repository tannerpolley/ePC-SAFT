from __future__ import annotations

import math

import numpy as np
import pytest

import epcsaft
from epcsaft import _core


def _toy_mixture() -> epcsaft.ePCSAFTMixture:
    return epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )


def test_native_chemical_equilibrium_residual_evaluator_rejects_removed_backend() -> None:
    mix = _toy_mixture()
    request = {
        "T": 298.15,
        "P": 1.0e5,
        "initial_x": [0.5, 0.5],
        "balance_matrix": [1.0, 1.0],
        "balance_rows": 1,
        "total_vector": [1.0],
        "reaction_stoichiometry": [-1.0, 1.0],
        "reaction_rows": 1,
        "log_equilibrium_constants": [math.log(3.0)],
        "reaction_standard_states": [1],
        "options": {"jacobian_backend": "finite" + "_difference", "finite" + "_difference_step": 1.0e-7},
    }

    with pytest.raises(_core.NativeValueError, match="jacobian_backend"):
        _core._evaluate_chemical_equilibrium_residual_native(mix._native, request)


def test_native_chemical_equilibrium_solve_rejects_unimplemented_cppad_nlp_derivatives() -> None:
    mix = _toy_mixture()
    request = {
        "T": 298.15,
        "P": 1.0e5,
        "initial_x": [0.5, 0.5],
        "balance_matrix": [1.0, 1.0],
        "balance_rows": 1,
        "total_vector": [1.0],
        "reaction_stoichiometry": [-1.0, 1.0],
        "reaction_rows": 1,
        "log_equilibrium_constants": [math.log(3.0)],
        "reaction_standard_states": [1],
        "options": {"solver_backend": "ipopt", "jacobian_backend": "cppad"},
    }

    with pytest.raises(_core.NativeValueError, match="CppAD chemical-equilibrium NLP derivatives"):
        _core._solve_chemical_equilibrium_native(mix._native, request)


def test_mixture_equilibrium_rejects_non_native_chemical_equilibrium_backend() -> None:
    mix = _toy_mixture()

    with pytest.raises(epcsaft.InputError, match="native"):
        mix.equilibrium(
            kind="chemical_equilibrium",
            T=298.15,
            P=1.0e5,
            z=[0.5, 0.5],
            balances={"total": {"A": 1.0, "B": 1.0}},
            totals={"total": 1.0},
            reactions=[epcsaft.ReactionDefinition({"A": -1.0, "B": 1.0}, math.log(3.0))],
            backend="legacy",
        )
