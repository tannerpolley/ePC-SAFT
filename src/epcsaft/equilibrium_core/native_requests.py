"""Internal builders for native equilibrium request payloads."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np


_CHEMICAL_POTENTIAL_TOLERANCE_FLOOR = 1.0e-7
_PHASE_DISTANCE_TOLERANCE_FLOOR = 1.0e-4


def neutral_two_phase_eos_tolerances(P: float, options: Any) -> tuple[float, float, float, float]:
    material_tolerance = float(options.tolerance)
    pressure_tolerance = max(abs(float(P)) * material_tolerance, material_tolerance)
    chemical_potential_tolerance = max(material_tolerance, _CHEMICAL_POTENTIAL_TOLERANCE_FLOOR)
    phase_distance_tolerance = max(10.0 * float(options.min_composition), _PHASE_DISTANCE_TOLERANCE_FLOOR)
    return material_tolerance, pressure_tolerance, chemical_potential_tolerance, phase_distance_tolerance


def build_reactive_two_phase_eos_native_request(
    *,
    T: float,
    P: float,
    feed: Any,
    balance_matrix: Any,
    total_vector: Any,
    species: Sequence[str],
    reactions: Sequence[Any],
    reaction_phase_stoichiometry: Any = None,
) -> dict[str, Any]:
    labels = [str(label) for label in species]
    balance_values = np.asarray(balance_matrix, dtype=float)
    reaction_matrix = np.asarray(
        [[float(reaction.stoichiometry.get(label, 0.0)) for label in labels] for reaction in reactions],
        dtype=float,
    )
    phase_stoichiometry_values: list[float] = []
    if reaction_phase_stoichiometry is not None:
        phase_stoichiometry_values = np.asarray(reaction_phase_stoichiometry, dtype=float).reshape(-1).tolist()
    return {
        "T": float(T),
        "P": float(P),
        "feed_amounts": np.asarray(feed, dtype=float).tolist(),
        "balance_matrix": balance_values.reshape(-1).tolist(),
        "balance_rows": int(balance_values.shape[0]),
        "total_vector": np.asarray(total_vector, dtype=float).tolist(),
        "reaction_stoichiometry": reaction_matrix.reshape(-1).tolist(),
        "reaction_rows": int(reaction_matrix.shape[0]),
        "log_equilibrium_constants": [float(reaction.log_equilibrium_constant) for reaction in reactions],
        "reaction_standard_states": [reaction.convention.native_standard_state_code for reaction in reactions],
        "reaction_phase_stoichiometry": phase_stoichiometry_values,
    }


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
    continuation_state: Any = None,
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
            "hessian_mode": str(options.hessian_mode),
            "iteration_history_limit": int(options.ipopt_iteration_history_limit),
            "linear_solver": str(options.ipopt_linear_solver),
            "acceptable_tolerance": float(
                options.ipopt_acceptable_tolerance
                if options.ipopt_acceptable_tolerance is not None
                else max(100.0 * float(options.tolerance), 1.0e-10)
            ),
            "constraint_violation_tolerance": float(
                options.ipopt_constraint_violation_tolerance
                if options.ipopt_constraint_violation_tolerance is not None
                else float(options.tolerance)
            ),
            "dual_infeasibility_tolerance": float(
                options.ipopt_dual_infeasibility_tolerance
                if options.ipopt_dual_infeasibility_tolerance is not None
                else float(options.tolerance)
            ),
            "complementarity_tolerance": float(
                options.ipopt_complementarity_tolerance
                if options.ipopt_complementarity_tolerance is not None
                else float(options.tolerance)
            ),
            "continuation_state": None if continuation_state is None else dict(continuation_state),
        },
    }
