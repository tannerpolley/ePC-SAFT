from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from epcsaft import _core, ePCSAFTMixture


def _reactive_electrolyte_lle_fixture() -> tuple[
    ePCSAFTMixture,
    np.ndarray,
    epcsaft.ReactionDefinition,
]:
    species = ["H2O", "Butanol", "Na+", "Cl-"]
    aq = np.asarray([0.798324680201737, 0.016320352824141723, 0.09267748348706063, 0.09267748348706063])
    org = np.asarray([0.37006036048879404, 0.6214918588210971, 0.004223890345054407, 0.004223890345054407])
    beta_org = 0.613766575013417
    feed = (1.0 - beta_org) * aq + beta_org * org
    mix = ePCSAFTMixture.from_dataset("2022_Ascani", species, feed, 298.15)
    reaction = epcsaft.ReactionDefinition.from_literature_constant(
        {"H2O": -1.0, "Butanol": 1.0},
        log_equilibrium_constant=-1.265500953237746,
        name="water_to_butanol",
        standard_state="mole_fraction_activity",
        source="repo-contained model-consistent reactive electrolyte LLE fixture",
    )
    return mix, feed, reaction


def _assert_reactive_phase_native_ipopt_gate(excinfo: pytest.ExceptionInfo[epcsaft.InputError]) -> None:
    assert "native Ipopt reactive phase-equilibrium NLP route" in str(excinfo.value)


def test_reactive_electrolyte_lle_public_route_requires_native_ipopt(monkeypatch) -> None:
    mix, feed, reaction = _reactive_electrolyte_lle_fixture()
    captured: dict[str, object] = {}

    def native_gate(*args):
        captured["args"] = args
        return {
            "backend": "ipopt",
            "compiled": False,
            "adapter_available": False,
            "problem_name": "reactive_two_phase_eos",
            "derivative_backend": "analytic_cppad",
            "accepted": False,
            "status": "ipopt_dependency_required",
            "postsolve": {},
        }

    monkeypatch.setattr(_core, "_native_reactive_electrolyte_lle_eos_route_result", native_gate)
    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="reactive_electrolyte_lle",
            T=298.15,
            P=1.013e5,
            z=feed,
            balances={
                "solvent_total": {"H2O": 1.0, "Butanol": 1.0},
                "sodium": {"Na+": 1.0},
                "chloride": {"Cl-": 1.0},
            },
            totals={
                "solvent_total": float(feed[0] + feed[1]),
                "sodium": float(feed[2]),
                "chloride": float(feed[3]),
            },
            reactions=[reaction],
            phase_options=epcsaft.EquilibriumOptions(
                max_iterations=80,
                tolerance=1.0e-8,
                min_composition=1.0e-12,
            ),
        )

    _assert_reactive_phase_native_ipopt_gate(excinfo)
    args = captured["args"]
    assert args[0] is mix._native
    assert args[3] == pytest.approx(feed.tolist())
    assert args[4] == 3
    assert args[5] == pytest.approx([1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0])
    assert args[6] == pytest.approx([float(feed[0] + feed[1]), float(feed[2]), float(feed[3])])
    assert args[7] == 1
    assert args[8] == pytest.approx([-1.0, 1.0, 0.0, 0.0])
