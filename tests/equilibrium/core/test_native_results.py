from __future__ import annotations

import json

import numpy as np
import pytest

import epcsaft
from epcsaft.equilibrium_core.native_results import neutral_two_phase_payload_to_result


def _accepted_native_payload() -> dict[str, object]:
    diagnostics = {
        "derivative_backend": "analytic_cppad",
        "rejection_reason": "accepted",
        "material_balance_norm": 1.0e-12,
        "pressure_consistency_norm": 2.0e-8,
        "chemical_potential_consistency_norm": 3.0e-9,
        "phase_distance": 0.4,
    }
    return {
        "accepted": True,
        "backend": "native_equilibrium_nlp",
        "problem_kind": "neutral_two_phase_eos",
        "phase_labels": ["phase_0", "phase_1"],
        "stable": False,
        "split_detected": True,
        "diagnostics": diagnostics,
        "phases": [
            {
                "label": "phase_0",
                "composition": [0.7, 0.3],
                "density": 120.0,
                "temperature": 300.0,
                "pressure": 2.0e5,
                "phase_fraction": 0.5,
                "ln_fugacity_coefficient": [-0.01, 0.02],
                "fugacity_coefficient": np.exp([-0.01, 0.02]).tolist(),
                "diagnostics": {"volume": 1.0e-2},
            },
            {
                "label": "phase_1",
                "composition": [0.1, 0.9],
                "density": 120.0,
                "temperature": 300.0,
                "pressure": 2.0e5,
                "phase_fraction": 0.5,
                "ln_fugacity_coefficient": [0.03, -0.04],
                "fugacity_coefficient": np.exp([0.03, -0.04]).tolist(),
                "diagnostics": {"volume": 1.0e-2},
            },
        ],
    }


def test_neutral_native_payload_converts_to_public_equilibrium_result() -> None:
    result = neutral_two_phase_payload_to_result(_accepted_native_payload())

    assert isinstance(result, epcsaft.EquilibriumResult)
    assert result.backend == "native_equilibrium_nlp"
    assert result.problem_kind == "neutral_two_phase_eos"
    assert result.phase_labels == ["phase_0", "phase_1"]
    assert result.stable is False
    assert result.split_detected is True
    assert result.diagnostics["derivative_backend"] == "analytic_cppad"
    assert len(result.phases) == 2
    assert result.phases[0].label == "phase_0"
    assert np.allclose(result.phases[0].composition, [0.7, 0.3])
    assert np.allclose(result.phases[0].ln_fugacity_coefficient, [-0.01, 0.02])
    assert np.allclose(result.phases[0].fugacity_coefficient, np.exp(result.phases[0].ln_fugacity_coefficient))
    json.dumps(result.to_dict(), allow_nan=False)


def test_neutral_native_payload_rejection_raises_solution_error() -> None:
    payload = _accepted_native_payload()
    payload["accepted"] = False
    payload["rejection_reason"] = "phase_distance"
    payload["diagnostics"] = {"rejection_reason": "phase_distance", "phase_distance": 0.0}
    payload["phases"] = []

    with pytest.raises(epcsaft.SolutionError, match="rejected") as exc_info:
        neutral_two_phase_payload_to_result(payload)

    assert exc_info.value.diagnostics["rejection_reason"] == "phase_distance"
