"""Validation-local phase-state payload helpers."""

from __future__ import annotations

from typing import Any

import numpy as np

from epcsaft._types import InputError, SolutionError
from epcsaft.equilibrium import EquilibriumOptions


def phase_state(
    mixture: Any,
    T: float,
    P: float,
    composition: np.ndarray,
    label: str,
    options: EquilibriumOptions,
    context: str,
) -> dict[str, Any]:
    try:
        state = mixture.state(T=T, P=P, x=composition, phase=label)
    except SolutionError:
        raise
    except (InputError, ValueError, RuntimeError, ArithmeticError) as exc:
        raise SolutionError(f"Failed to construct {label} phase during {context}: {exc}") from exc
    diagnostics = state.state_diagnostics(species=mixture.species) if options.include_phase_diagnostics else None
    return {
        "state": state,
        "ln_phi": np.asarray(state.fugacity_coefficient(), dtype=float),
        "density": float(state.density()),
        "diagnostics": diagnostics,
    }
