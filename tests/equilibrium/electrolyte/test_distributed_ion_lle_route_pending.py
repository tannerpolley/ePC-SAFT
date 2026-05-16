from __future__ import annotations

import pytest

import epcsaft
from tests.equilibrium.electrolyte.test_electrolyte_lle_smokes import _assert_electrolyte_lle_route_pending
from tests.equilibrium.electrolyte.test_electrolyte_lle_solver_contracts import _case2_feed, _case2_mixture


def test_distributed_ion_lle_public_route_requires_native_ipopt() -> None:
    feed = _case2_feed()
    mix = _case2_mixture(feed)

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="electrolyte_lle",
            T=298.15,
            P=1.0e5,
            z=feed,
            options=epcsaft.EquilibriumOptions(max_iterations=180, tolerance=1.0e-8),
        )

    _assert_electrolyte_lle_route_pending(excinfo)
