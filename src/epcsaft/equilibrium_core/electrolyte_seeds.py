"""Helpers for constructing charge-neutral electrolyte LLE initial phases."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from epcsaft._types import InputError


@dataclass(frozen=True, slots=True)
class ElectrolyteInitialPhases:
    aq: np.ndarray
    org: np.ndarray
    phase_fraction: float
    diagnostics: dict[str, float | str]

    def to_initial_phases(self) -> dict[str, object]:
        return {
            "aq": self.aq.copy(),
            "org": self.org.copy(),
            "phase_fraction": float(self.phase_fraction),
        }


def charge_neutral_lle_seed_from_org_phase(
    feed: Sequence[float],
    org_phase: Sequence[float],
    phase_fraction: float,
    charges: Sequence[float],
    *,
    min_composition: float = 1.0e-12,
    seed_name: str = "org_phase_complement",
) -> ElectrolyteInitialPhases:
    feed_arr = _normalize_vector(feed, "feed", min_composition)
    org_arr = _normalize_vector(org_phase, "org_phase", min_composition)
    charges_arr = np.asarray(charges, dtype=float).flatten()
    beta = float(phase_fraction)
    if charges_arr.size != feed_arr.size:
        raise InputError("charges length must match feed length for electrolyte_lle initial phase construction.")
    if not np.isfinite(beta) or beta <= 0.0 or beta >= 1.0:
        raise InputError("phase_fraction must be > 0 and < 1 for electrolyte_lle initial phase construction.")
    org_charge = abs(float(np.dot(org_arr, charges_arr)))
    if org_charge > 1.0e-10:
        raise InputError("org_phase must be charge neutral for electrolyte_lle initial phase construction.")
    aq = (feed_arr - beta * org_arr) / (1.0 - beta)
    if np.any(aq < min_composition):
        raise InputError("org_phase and phase_fraction produce a non-positive aqueous complement.")
    aq = _normalize_vector(aq, "aq_complement", min_composition)
    material_error = float(np.max(np.abs((1.0 - beta) * aq + beta * org_arr - feed_arr)))
    aq_charge = abs(float(np.dot(aq, charges_arr)))
    if aq_charge > 1.0e-8:
        raise InputError("constructed aq phase is not charge neutral for electrolyte_lle.")
    return ElectrolyteInitialPhases(
        aq=aq,
        org=org_arr,
        phase_fraction=beta,
        diagnostics={
            "seed_name": seed_name,
            "initial_phase_material_balance_error": material_error,
            "initial_phase_charge_balance_error": max(aq_charge, org_charge),
        },
    )


def solvent_endpoint_seed(
    feed: Sequence[float],
    charges: Sequence[float],
    species: Sequence[str],
    organic_components: Mapping[str, float],
    salt_in_org: float,
    phase_fraction: float,
    *,
    min_composition: float = 1.0e-12,
) -> ElectrolyteInitialPhases:
    species_list = [str(name) for name in species]
    charges_arr = np.asarray(charges, dtype=float).flatten()
    org = np.full(len(species_list), min_composition, dtype=float)
    for name, amount in organic_components.items():
        org[species_list.index(name)] = float(amount)
    cation_indices = [i for i, charge in enumerate(charges_arr) if charge > 0.0]
    anion_indices = [i for i, charge in enumerate(charges_arr) if charge < 0.0]
    if len(cation_indices) != 1 or len(anion_indices) != 1:
        raise InputError("solvent_endpoint_seed currently supports one cation and one anion.")
    org[cation_indices[0]] = float(salt_in_org)
    org[anion_indices[0]] = float(salt_in_org)
    total = float(np.sum(org))
    if total <= 0.0:
        raise InputError("organic endpoint seed has non-positive total composition.")
    return charge_neutral_lle_seed_from_org_phase(
        feed,
        org / total,
        phase_fraction,
        charges_arr,
        min_composition=min_composition,
        seed_name="solvent_endpoint_seed",
    )


def _normalize_vector(values: Sequence[float], label: str, min_composition: float) -> np.ndarray:
    arr = np.asarray(values, dtype=float).flatten()
    if arr.size == 0:
        raise InputError(f"{label} must not be empty.")
    if not np.all(np.isfinite(arr)):
        raise InputError(f"{label} must contain only finite values.")
    if np.any(arr < 0.0):
        raise InputError(f"{label} must be non-negative.")
    total = float(np.sum(arr))
    if total <= 0.0:
        raise InputError(f"{label} must have a positive sum.")
    arr = arr / total
    arr = np.maximum(arr, min_composition)
    return arr / float(np.sum(arr))
