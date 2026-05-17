"""Internal builders for native equilibrium request payloads."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np


def build_reactive_speciation_native_request(
    *,
    T: float,
    P: float,
    initial_x: Any,
    balance_matrix: Any,
    total_vector: Any,
    species: Sequence[str],
    reactions: Sequence[Any],
    options: Any,
) -> dict[str, Any]:
    labels = [str(label) for label in species]
    balance_values = np.asarray(balance_matrix, dtype=float)
    reaction_matrix = np.asarray(
        [[float(reaction.stoichiometry.get(label, 0.0)) for label in labels] for reaction in reactions],
        dtype=float,
    )
    return {
        "T": float(T),
        "P": float(P),
        "initial_x": np.asarray(initial_x, dtype=float).tolist(),
        "balance_matrix": balance_values.reshape(-1).tolist(),
        "balance_rows": int(balance_values.shape[0]),
        "total_vector": np.asarray(total_vector, dtype=float).tolist(),
        "reaction_stoichiometry": reaction_matrix.reshape(-1).tolist(),
        "reaction_rows": int(reaction_matrix.shape[0]),
        "log_equilibrium_constants": [float(reaction.log_equilibrium_constant) for reaction in reactions],
        "reaction_standard_states": [reaction.convention.native_standard_state_code for reaction in reactions],
        "options": {
            "max_iterations": int(options.max_iterations),
            "tolerance": float(options.tolerance),
            "min_mole_fraction": float(options.min_mole_fraction),
            "jacobian_backend": str(options.jacobian_backend),
            "solver_backend": str(options.solver_backend),
            "phase": str(options.phase),
            "activity_output": str(options.activity_output),
        },
    }
