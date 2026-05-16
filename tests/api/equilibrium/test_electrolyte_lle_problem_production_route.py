from __future__ import annotations

import numpy as np

import epcsaft
from tests.equilibrium.electrolyte.test_salting_out_lle_benchmark import _salting_out_fixture
from tests.helpers.numeric import assert_allclose


def test_electrolyte_lle_problem_routes_to_native_ceres_production_solver() -> None:
    mix, feed, initial_phases = _salting_out_fixture()
    problem = epcsaft.ElectrolyteLLEProblem(
        T=298.15,
        P=1.013e5,
        z=feed,
        initial_phases=initial_phases,
        options=epcsaft.EquilibriumOptions(max_iterations=80, tolerance=1.0e-8, min_composition=1.0e-12),
    )

    result = mix.solve_equilibrium(problem)
    diagnostics = result.diagnostics

    assert result.backend == "electrolyte_lle"
    assert result.split_detected is True
    assert diagnostics["equilibrium_route"] == "electrolyte_lle"
    assert diagnostics["solver_backend"] == "ceres"
    assert diagnostics["solver_method"] == "ceres_trust_region_residual_solve"
    assert diagnostics["jacobian_backend"] == "cppad_implicit"
    assert diagnostics["derivative_backend"] == "cppad_implicit"
    assert diagnostics["jacobian_available"] is True
    assert diagnostics["derivative_available"] is True
    assert diagnostics["neutral_fugacity_residual_norm"] <= 1.0e-8
    assert diagnostics["ionic_equilibrium_residual_norm"] <= 1.0e-8

    reconstructed = np.zeros_like(feed)
    for phase in result.phases:
        reconstructed += phase.phase_fraction * phase.composition
    assert_allclose(reconstructed, feed, atol=1.0e-10)
