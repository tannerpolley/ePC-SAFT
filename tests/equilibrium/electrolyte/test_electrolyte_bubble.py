from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from epcsaft import InputError


def _salt_mixture() -> epcsaft.ePCSAFTMixture:
    x = np.asarray([0.98, 0.01, 0.01], dtype=float)
    return epcsaft.ePCSAFTMixture.from_dataset("2026_Khudaida", ["H2O", "Na+", "Cl-"], x, 298.15)


def test_electrolyte_bubble_pressure_requires_native_ipopt_route_builder() -> None:
    mix = _salt_mixture()

    with pytest.raises(InputError, match="native Ipopt equilibrium route builder"):
        mix.equilibrium(
            kind="electrolyte_bubble_pressure",
            T=298.15,
            x_liq=[0.98, 0.01, 0.01],
            vapor_species=["H2O"],
            volatile_species=["H2O"],
            nonvolatile_species=["Na+", "Cl-"],
            backend="native",
        )


def test_electrolyte_bubble_pressure_rejects_python_backend_alias() -> None:
    mix = _salt_mixture()

    with pytest.raises(InputError, match="backend must be None or 'native'"):
        mix.equilibrium(
            kind="electrolyte_bubble_pressure",
            T=298.15,
            x_liq=[0.98, 0.01, 0.01],
            vapor_species=["H2O"],
            backend="python",
        )
