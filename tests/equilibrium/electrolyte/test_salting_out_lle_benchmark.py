from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from epcsaft import ePCSAFTMixture
from tests.equilibrium.electrolyte.test_electrolyte_lle_smokes import _assert_electrolyte_lle_route_pending


def _salting_out_fixture() -> tuple[ePCSAFTMixture, np.ndarray, dict[str, object]]:
    species = ["H2O", "Butanol", "Na+", "Cl-"]
    aq = np.asarray([0.798324680201737, 0.016320352824141723, 0.09267748348706063, 0.09267748348706063])
    org = np.asarray([0.37006036048879404, 0.6214918588210971, 0.004223890345054407, 0.004223890345054407])
    beta_org = 0.613766575013417
    feed = (1.0 - beta_org) * aq + beta_org * org
    mix = ePCSAFTMixture.from_dataset("2022_Ascani", species, feed, 298.15)
    initial_phases = {"aq": aq, "org": org, "phase_fraction": beta_org}
    return mix, feed, initial_phases


def test_quaternary_salting_out_lle_benchmark_requires_native_ipopt_route() -> None:
    mix, feed, initial_phases = _salting_out_fixture()

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="electrolyte_lle",
            T=298.15,
            P=1.013e5,
            z=feed,
            initial_phases=initial_phases,
            options=epcsaft.EquilibriumOptions(max_iterations=80, tolerance=1.0e-8, min_composition=1.0e-12),
        )

    _assert_electrolyte_lle_route_pending(excinfo)
