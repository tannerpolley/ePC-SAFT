from __future__ import annotations

from dataclasses import dataclass

import pytest

from epcsaft.equilibrium_core.native_requests import neutral_two_phase_eos_tolerances


@dataclass(frozen=True)
class _Options:
    tolerance: float
    min_composition: float


def test_neutral_two_phase_eos_tolerances_scale_pressure_and_phase_distance() -> None:
    tolerances = neutral_two_phase_eos_tolerances(2.0e5, _Options(tolerance=1.0e-6, min_composition=1.0e-12))

    assert tolerances == pytest.approx((1.0e-6, 0.2, 1.0e-6, 1.0e-4))


def test_neutral_two_phase_eos_tolerances_keep_solver_away_from_collapsed_lle_phase() -> None:
    tolerances = neutral_two_phase_eos_tolerances(5.0e6, _Options(tolerance=1.0e-8, min_composition=1.0e-12))

    assert tolerances == pytest.approx((1.0e-8, 5.0e-2, 1.0e-7, 1.0e-4))
