"""External dataset bundle validation helpers."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np

from .parameters import get_prop_dict


def validate_dataset_bundle(
    dataset: Mapping[str, Any] | str | Path,
    *,
    species: Sequence[str],
    x: Any | None = None,
    T: float = 298.15,
    reactions: Sequence[Any] | None = None,
    user_options: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Return a structured validation report for an external parameter bundle."""

    labels = [str(item) for item in species]
    params = _load_params(dataset, labels, x=x, T=T, user_options=user_options)
    errors: list[str] = []
    warnings: list[str] = []
    _check_vector_lengths(params, labels, errors)
    charges = _check_charges(params, labels, errors, warnings)
    _check_born(params, labels, charges, errors, warnings)
    _check_binary_matrices(params, labels, errors, warnings)
    _check_reactions(labels, reactions or (), errors)
    report = {
        "valid": not errors,
        "species": labels,
        "component_count": len(labels),
        "errors": errors,
        "warnings": warnings,
        "provenance": {
            "dataset": str(dataset) if not isinstance(dataset, Mapping) else "mapping",
            "temperature": float(T),
        },
    }
    return report


def _load_params(
    dataset: Mapping[str, Any] | str | Path,
    species: list[str],
    *,
    x: Any | None,
    T: float,
    user_options: Mapping[str, Any] | None,
) -> Mapping[str, Any]:
    if isinstance(dataset, Mapping):
        return dataset
    composition = np.full(len(species), 1.0 / max(len(species), 1), dtype=float) if x is None else x
    return get_prop_dict(dataset, species, composition, T, user_options=user_options)


def _check_vector_lengths(params: Mapping[str, Any], species: list[str], errors: list[str]) -> None:
    n = len(species)
    for key in ("m", "s", "e"):
        if key not in params:
            errors.append(f"missing required parameter vector '{key}'")
            continue
        values = np.asarray(params[key], dtype=float).flatten()
        if values.size != n:
            errors.append(f"parameter vector '{key}' length {values.size} does not match species length {n}")
        elif np.any(~np.isfinite(values)):
            errors.append(f"parameter vector '{key}' contains non-finite values")


def _check_charges(
    params: Mapping[str, Any],
    species: list[str],
    errors: list[str],
    warnings: list[str],
) -> np.ndarray:
    if "z" not in params or np.asarray(params.get("z"), dtype=float).size == 0:
        warnings.append("charge vector 'z' is missing or empty; treating all species as neutral")
        return np.zeros(len(species), dtype=float)
    charges = np.asarray(params["z"], dtype=float).flatten()
    if charges.size != len(species):
        errors.append(f"charge vector 'z' length {charges.size} does not match species length {len(species)}")
        return np.zeros(len(species), dtype=float)
    for label, charge in zip(species, charges):
        expected = _charge_hint(label)
        if expected is not None and np.sign(charge) != np.sign(expected):
            errors.append(f"charge sign mismatch for {label}: expected {expected:+.0f}, got {charge:+g}")
    return charges


def _check_born(
    params: Mapping[str, Any],
    species: list[str],
    charges: np.ndarray,
    errors: list[str],
    warnings: list[str],
) -> None:
    if not np.any(np.abs(charges) > 1.0e-12):
        return
    if "d_born" not in params:
        warnings.append("charged species present but d_born vector is missing")
        return
    values = np.asarray(params["d_born"], dtype=float).flatten()
    if values.size != len(species):
        errors.append(f"d_born length {values.size} does not match species length {len(species)}")
        return
    for label, charge, value in zip(species, charges, values):
        if abs(charge) > 1.0e-12 and not np.isfinite(value):
            errors.append(f"non-finite d_born for charged species {label}")


def _check_binary_matrices(
    params: Mapping[str, Any],
    species: list[str],
    errors: list[str],
    warnings: list[str],
) -> None:
    n = len(species)
    for key in ("k_ij", "l_ij", "k_hb", "k_hb_ij"):
        if key not in params:
            if key == "k_ij":
                warnings.append("binary interaction matrix 'k_ij' is missing")
            continue
        matrix = np.asarray(params[key], dtype=float)
        if matrix.shape != (n, n):
            errors.append(f"binary interaction matrix '{key}' shape {matrix.shape} does not match ({n}, {n})")
        elif np.any(~np.isfinite(matrix)):
            errors.append(f"binary interaction matrix '{key}' contains non-finite values")


def _check_reactions(species: list[str], reactions: Sequence[Any], errors: list[str]) -> None:
    species_set = set(species)
    for reaction in reactions:
        stoich = getattr(reaction, "stoichiometry", {})
        for label in stoich:
            if str(label) not in species_set:
                errors.append(f"Unknown species '{label}' in reaction stoichiometry")


def _charge_hint(label: str) -> float | None:
    token = label.strip()
    if token.endswith("+"):
        return 1.0
    if token.endswith("-"):
        return -1.0
    if "^2-" in token or "2-" in token:
        return -1.0
    if "^2+" in token or "2+" in token:
        return 1.0
    return None
