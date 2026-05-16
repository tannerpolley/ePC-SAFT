from __future__ import annotations

import pytest

import epcsaft
from tests.equilibrium.electrolyte.test_electrolyte_lle_smokes import _assert_electrolyte_lle_route_pending
from tests.equilibrium.electrolyte.test_salting_out_lle_benchmark import _salting_out_fixture


def test_electrolyte_lle_problem_requires_native_ipopt_route() -> None:
    mix, feed, initial_phases = _salting_out_fixture()
    problem = epcsaft.ElectrolyteLLEProblem(
        T=298.15,
        P=1.013e5,
        z=feed,
        initial_phases=initial_phases,
        options=epcsaft.EquilibriumOptions(max_iterations=80, tolerance=1.0e-8, min_composition=1.0e-12),
    )

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.solve_equilibrium(problem)

    _assert_electrolyte_lle_route_pending(excinfo)
