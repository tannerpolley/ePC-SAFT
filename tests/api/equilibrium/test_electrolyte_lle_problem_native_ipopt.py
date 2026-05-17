from __future__ import annotations

import numpy as np

import epcsaft
from epcsaft import ePCSAFTMixture


def _typed_problem_fixture() -> tuple[ePCSAFTMixture, np.ndarray]:
    species = ["H2O", "Butanol", "Na+", "Cl-"]
    aq = np.asarray([0.798324680201737, 0.016320352824141723, 0.09267748348706063, 0.09267748348706063])
    org = np.asarray([0.37006036048879404, 0.6214918588210971, 0.004223890345054407, 0.004223890345054407])
    beta_org = 0.613766575013417
    feed = (1.0 - beta_org) * aq + beta_org * org
    mix = ePCSAFTMixture.from_dataset("2022_Ascani", species, feed, 298.15)
    return mix, feed


def test_electrolyte_lle_problem_executes_native_ipopt_route() -> None:
    mix, feed = _typed_problem_fixture()
    problem = epcsaft.ElectrolyteLLEProblem(
        T=298.15,
        P=1.013e5,
        z=feed,
        options=epcsaft.EquilibriumOptions(max_iterations=500, tolerance=1.0e-8, min_composition=1.0e-12),
    )

    result = mix.solve_equilibrium(problem)

    assert result.problem_kind == "electrolyte_lle"
    assert result.diagnostics["route_status"] == "accepted"
    assert result.diagnostics["solver_status"] == "success"
    assert result.diagnostics["solver_backend"] == "ipopt"
