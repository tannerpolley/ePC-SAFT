from __future__ import annotations

import numpy as np
import pytest

from epcsaft import _core


def test_association_mass_action_block_reports_exact_residual_and_jacobians() -> None:
    density = 0.5
    site_fractions = np.asarray([0.8, 0.6], dtype=float)
    site_composition = np.asarray([0.25, 0.75], dtype=float)
    delta = np.asarray([[0.0, 2.0], [3.0, 0.0]], dtype=float)

    payload = _core._native_association_mass_action_block(
        density,
        site_fractions.tolist(),
        site_composition.tolist(),
        delta.ravel().tolist(),
    )

    association_sums = density * (delta @ (site_composition * site_fractions))
    expected_residuals = site_fractions * (1.0 + association_sums) - 1.0
    expected_site_jacobian = np.diag(1.0 + association_sums) + density * site_fractions[:, None] * delta * site_composition
    expected_composition_jacobian = density * site_fractions[:, None] * delta * site_fractions
    expected_density_derivative = site_fractions * (delta @ (site_composition * site_fractions))

    assert payload["block"] == "association_mass_action"
    assert payload["derivative_backend"] == "analytic"
    assert payload["constraint_names"] == ["association_site_0", "association_site_1"]
    assert payload["residuals"] == pytest.approx(expected_residuals.tolist(), rel=1.0e-14, abs=1.0e-14)
    assert np.asarray(payload["site_fraction_jacobian_row_major"], dtype=float).reshape((2, 2)) == pytest.approx(
        expected_site_jacobian,
        rel=1.0e-14,
        abs=1.0e-14,
    )
    assert np.asarray(payload["site_composition_jacobian_row_major"], dtype=float).reshape((2, 2)) == pytest.approx(
        expected_composition_jacobian,
        rel=1.0e-14,
        abs=1.0e-14,
    )
    assert payload["density_derivative"] == pytest.approx(expected_density_derivative.tolist())
