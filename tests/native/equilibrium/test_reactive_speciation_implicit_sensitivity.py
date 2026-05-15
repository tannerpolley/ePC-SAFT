from __future__ import annotations

import math

import numpy as np
import pytest

import epcsaft


def test_reactive_speciation_reports_solved_state_derivative_boundaries() -> None:
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )

    result = epcsaft.solve_reactive_speciation(
        species=["A", "B"],
        mixture_factory=lambda x, T, P: mix,
        T=298.15,
        P=1.0e5,
        balances={"total": {"A": 1.0, "B": 1.0}},
        totals={"total": 1.0},
        reactions=[
            epcsaft.ReactionDefinition.from_literature_constant(
                {"A": -1.0, "B": 1.0},
                log_equilibrium_constant=math.log(3.0),
                standard_state="ideal_mole_fraction",
            )
        ],
        initial_x=[0.5, 0.5],
    )

    diagnostics = result.diagnostics
    assert diagnostics["derivative_policy"]["numerical_derivative_backend_available"] is False
    assert diagnostics["derivative_policy"]["unsupported_derivative_status"] == "not_available"
    assert "reactive_speciation_variables" in diagnostics["solved_state_derivative_blocks"]
    assert "association_site_fractions" in diagnostics["solved_state_derivative_blocks"]
    assert diagnostics["derivative_backend_by_block"]["reactive_speciation_variables"] == "analytic_implicit"
    assert diagnostics["derivative_backend_by_block"]["association_site_fractions"] in {
        "analytic_implicit",
        "cppad_implicit",
        "not_available",
    }
    implicit_results = diagnostics["implicit_solve_results"]
    assert implicit_results["reactive_speciation_variables"]["backend"] == "analytic_implicit"
    assert implicit_results["reactive_speciation_variables"]["status"] == "residual_jacobian_available"
    backend_values = {str(value).lower() for value in diagnostics["derivative_backend_by_block"].values()}
    assert "numerical_derivative" not in backend_values


def test_explicit_cppad_request_reports_not_available_without_numerical_derivative() -> None:
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )

    with pytest.raises(epcsaft.InputError, match="not_available"):
        epcsaft.solve_reactive_speciation(
            species=["A", "B"],
            mixture_factory=lambda x, T, P: mix,
            T=298.15,
            P=1.0e5,
            balances={"total": {"A": 1.0, "B": 1.0}},
            totals={"total": 1.0},
            reactions=[
                epcsaft.ReactionDefinition.from_literature_constant(
                    {"A": -1.0, "B": 1.0},
                    log_equilibrium_constant=math.log(3.0),
                    standard_state="ideal_mole_fraction",
                )
            ],
            initial_x=[0.5, 0.5],
            options=epcsaft.ReactiveSpeciationOptions(jacobian_backend="cppad"),
        )
