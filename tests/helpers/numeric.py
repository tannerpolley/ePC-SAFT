from __future__ import annotations

import numpy as np


def assert_allclose(actual, expected, rtol=1.0e-7, atol=0.0) -> None:
    actual_array = np.asarray(actual, dtype=float)
    expected_array = np.asarray(expected, dtype=float)
    if actual_array.shape != expected_array.shape:
        raise AssertionError(f"shape mismatch: {actual_array.shape} != {expected_array.shape}")
    if not np.allclose(actual_array, expected_array, rtol=rtol, atol=atol):
        diff = np.max(np.abs(actual_array - expected_array)) if actual_array.size else 0.0
        raise AssertionError(f"array mismatch: max abs diff {diff}")
