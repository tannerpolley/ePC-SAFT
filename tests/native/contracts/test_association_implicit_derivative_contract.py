from __future__ import annotations

import pytest

from epcsaft import _core
from epcsaft.epcsaft import create_struct
from tests.helpers.native_cases import _ionic_state


def _state():
    mix, _species, _pressure, density, temperature, composition = _ionic_state()
    return mix.state(T=temperature, x=composition, rho=density)


def test_active_association_reports_implicit_backend_not_direct_cppad() -> None:
    state = _state()

    association_rows = [row for row in state.derivative_coverage_matrix() if row["quantity"] == "association"]

    assert association_rows
    assert association_rows[0]["backend"] in {"analytic_implicit", "backend_unavailable"}
    assert association_rows[0]["backend"] != "cppad"


def test_direct_cppad_eos_contribution_recording_rejects_active_association() -> None:
    mix, _species, _pressure, density, temperature, composition = _ionic_state()
    args = create_struct(mix.parameters)

    if not _core._native_cppad_smoke()["cppad_compiled"]:
        return

    with pytest.raises(_core.NativeValueError, match="backend_unavailable"):
        _core._native_cppad_eos_contributions(temperature, density, composition.tolist(), args)
