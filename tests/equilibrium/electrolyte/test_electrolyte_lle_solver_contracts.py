from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from epcsaft import ePCSAFTMixture
from tests.equilibrium.electrolyte.test_electrolyte_lle_smokes import _assert_electrolyte_lle_native_ipopt_gate


def _case2_feed() -> np.ndarray:
    return np.asarray(
        [
            0.940373242284748,
            0.04879624542603625,
            0.0019339313461782701,
            0.003481324798429627,
            0.005415256144607897,
        ],
        dtype=float,
    )


def _case2_mixture(feed=None) -> ePCSAFTMixture:
    if feed is None:
        feed = _case2_feed()
    return ePCSAFTMixture.from_dataset("2022_Ascani", ["H2O", "Butanol", "Na+", "K+", "Cl-"], feed, 298.15)


def _source_like_nacl_feed() -> np.ndarray:
    aq = np.asarray([0.798324680201737, 0.016320352824141723, 0.09267748348706063, 0.09267748348706063])
    org = np.asarray([0.37006036048879404, 0.6214918588210971, 0.004223890345054407, 0.004223890345054407])
    beta_org = 0.613766575013417
    return (1.0 - beta_org) * aq + beta_org * org


def test_source_like_nacl_feed_returns_sensible_native_ipopt_phase_split() -> None:
    feed = _source_like_nacl_feed()
    mix = ePCSAFTMixture.from_dataset("2022_Ascani", ["H2O", "Butanol", "Na+", "Cl-"], feed, 298.15)

    result = mix.equilibrium(
        kind="electrolyte_lle",
        T=298.15,
        P=1.013e5,
        z=feed,
        options=epcsaft.EquilibriumOptions(max_iterations=500, tolerance=1.0e-8, min_composition=1.0e-12),
    )

    assert result.phase_labels == ["aq", "org"]
    assert result.split_detected is True
    assert result.diagnostics["route_status"] == "accepted"
    assert result.diagnostics["solver_status"] == "success"
    assert result.diagnostics["application_status"] == "solve_succeeded"
    assert result.diagnostics["phase_distance"] > 0.1
    assert result.diagnostics["material_balance_norm"] <= 1.0e-8
    assert result.diagnostics["charge_balance_norm"] <= 1.0e-8
    aq, org = result.phases
    assert aq.composition[2] + aq.composition[3] > org.composition[2] + org.composition[3]
    assert org.composition[1] > aq.composition[1]


def test_ascani_case2_mixed_salt_executes_native_ipopt_route() -> None:
    mix = _case2_mixture()

    result = mix.equilibrium(
        kind="electrolyte_lle",
        T=298.15,
        P=1.0e5,
        z=_case2_feed(),
        options=epcsaft.EquilibriumOptions(max_iterations=180, tolerance=1.0e-8),
    )

    assert result.backend == "native_equilibrium_nlp"
    assert result.problem_kind == "electrolyte_lle"
    assert result.diagnostics["route_status"] == "accepted"
    assert result.diagnostics["solver_status"] == "success"
    assert result.diagnostics["application_status"] == "solve_succeeded"
    assert result.diagnostics["charge_balance_norm"] <= 1.0e-8

def test_auto_kind_routes_explicit_ionic_feed_to_native_ipopt_lle_gate() -> None:
    mix = _case2_mixture()

    result = mix.equilibrium(
        kind="auto",
        T=298.15,
        P=1.0e5,
        z=_case2_feed(),
        options=epcsaft.EquilibriumOptions(max_iterations=180, tolerance=1.0e-8),
    )

    assert result.problem_kind == "electrolyte_lle"
    assert result.diagnostics["route_status"] == "accepted"
    assert result.diagnostics["solver_backend"] == "ipopt"
