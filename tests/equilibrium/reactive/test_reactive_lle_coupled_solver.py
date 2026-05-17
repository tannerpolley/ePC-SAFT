from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from epcsaft import ePCSAFTMixture


def _neutral_reactive_lle_fixture() -> tuple[ePCSAFTMixture, np.ndarray, epcsaft.ReactionDefinition]:
    params = {
        "MW": np.asarray([32.042e-3, 84.147e-3]),
        "m": np.asarray([1.5255, 2.5303]),
        "s": np.asarray([3.2300, 3.8499]),
        "e": np.asarray([188.90, 278.11]),
        "e_assoc": np.asarray([2899.5, 0.0]),
        "vol_a": np.asarray([0.035176, 0.0]),
        "assoc_scheme": ["2B", None],
        "k_ij": np.asarray([[0.0, 0.051], [0.051, 0.0]]),
        "z": np.asarray([0.0, 0.0]),
        "dielc": np.asarray([33.05, 2.02]),
    }
    mix = ePCSAFTMixture.from_params(params, species=["Methanol", "Cyclohexane"])
    beta2 = 0.48813098468607985
    liq1 = np.asarray([0.11757838279937723, 0.8824216172006228])
    liq2 = np.asarray([0.7985874308392054, 0.20141256916079467])
    feed = (1.0 - beta2) * liq1 + beta2 * liq2
    reaction = epcsaft.ReactionDefinition.from_literature_constant(
        {"Methanol": -1.0, "Cyclohexane": 1.0},
        log_equilibrium_constant=-0.079259405371,
        name="methanol_to_cyclohexane",
        standard_state="mole_fraction_activity",
        source="repo-contained model-consistent reactive LLE fixture",
    )
    return mix, feed, reaction


def _assert_reactive_phase_route_pending(excinfo: pytest.ExceptionInfo[epcsaft.InputError]) -> None:
    message = str(excinfo.value)
    assert "native Ipopt reactive phase-equilibrium NLP route" in message


def test_neutral_reactive_lle_public_route_requires_native_ipopt() -> None:
    mix, feed, reaction = _neutral_reactive_lle_fixture()

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="reactive_lle",
            T=298.15,
            P=1.013e5,
            z=feed,
            balances={"total": {"Methanol": 1.0, "Cyclohexane": 1.0}},
            totals={"total": 1.0},
            reactions=[reaction],
            phase_options=epcsaft.EquilibriumOptions(max_iterations=80, tolerance=1.0e-8, min_composition=1.0e-12),
        )

    _assert_reactive_phase_route_pending(excinfo)
