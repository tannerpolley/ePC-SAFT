from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from epcsaft import ePCSAFTMixture
from tests.equilibrium.electrolyte.test_electrolyte_lle_smokes import _assert_electrolyte_lle_route_pending


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


def test_ascani_case2_mixed_salt_requires_native_ipopt_route() -> None:
    mix = _case2_mixture()

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="electrolyte_lle",
            T=298.15,
            P=1.0e5,
            z=_case2_feed(),
            options=epcsaft.EquilibriumOptions(max_iterations=180, tolerance=1.0e-8),
        )

    _assert_electrolyte_lle_route_pending(excinfo)


def test_auto_kind_routes_explicit_ionic_feed_to_pending_electrolyte_lle() -> None:
    mix = _case2_mixture()

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="auto",
            T=298.15,
            P=1.0e5,
            z=_case2_feed(),
            options=epcsaft.EquilibriumOptions(max_iterations=180, tolerance=1.0e-8),
        )

    _assert_electrolyte_lle_route_pending(excinfo)


def test_electrolyte_lle_solver_budget_options_validate_before_route_gate() -> None:
    mix = _case2_mixture()

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="electrolyte_lle",
            T=298.15,
            P=1.0e5,
            z=_case2_feed(),
            options=epcsaft.EquilibriumOptions(
                max_iterations=80,
                tolerance=1.0e-12,
                max_seed_attempts=1,
                max_total_objective_evaluations=1,
            ),
        )

    _assert_electrolyte_lle_route_pending(excinfo)


def test_electrolyte_lle_ignored_option_dict_validates_before_route_gate() -> None:
    mix = _case2_mixture()

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="electrolyte_lle",
            T=298.15,
            P=1.0e5,
            z=_case2_feed(),
            options={
                "max_nfev": 1,
                "solver_tol": 1.0e-12,
                "tpdf_global_trials": 1200,
                "tpdf_local_trials": 600,
                "charge_weight": 1000.0,
                "seed_x": [0.55, 0.40, 0.025, 0.025],
                "force_seed_solve": True,
            },
        )

    _assert_electrolyte_lle_route_pending(excinfo)
