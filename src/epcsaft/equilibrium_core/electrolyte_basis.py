"""Ascani-style transformed variables for electrolyte liquid splits."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from epcsaft._types import InputError


@dataclass(frozen=True, slots=True)
class ElectrolyteBasis:
    """Charge-constrained basis used by electrolyte phase-equilibrium solvers."""

    species: tuple[str, ...]
    charges: np.ndarray
    neutral_indices: tuple[int, ...]
    charged_indices: tuple[int, ...]
    cation_indices: tuple[int, ...]
    anion_indices: tuple[int, ...]
    e_matrix: np.ndarray
    salt_pairs: tuple[dict[str, Any], ...]
    formula_feed: np.ndarray
    rank: int
    variable_model: str = "ascani_transformed_salt_pairs"

    def to_dict(self) -> dict[str, Any]:
        return {
            "variable_model": self.variable_model,
            "neutral_indices": list(self.neutral_indices),
            "charged_indices": list(self.charged_indices),
            "cation_indices": list(self.cation_indices),
            "anion_indices": list(self.anion_indices),
            "e_matrix": self.e_matrix.tolist(),
            "salt_pairs": [dict(pair) for pair in self.salt_pairs],
            "formula_feed": self.formula_feed.tolist(),
            "basis_rank": int(self.rank),
        }


def build_electrolyte_basis(
    species: Any,
    charges: Any,
    feed: Any,
    *,
    salt_labels: tuple[str, ...] = (),
) -> ElectrolyteBasis:
    """Build the independent counterion-pair matrix for an explicit-ion feed."""
    labels = tuple(str(item) for item in species)
    z = np.asarray(charges, dtype=float).flatten()
    x = np.asarray(feed, dtype=float).flatten()
    if z.size != len(labels) or x.size != len(labels):
        raise InputError("species, charges, and feed must have matching lengths for electrolyte basis construction.")

    neutral_indices = tuple(int(i) for i, charge in enumerate(z) if abs(float(charge)) <= 1.0e-12)
    cation_indices = tuple(int(i) for i, charge in enumerate(z) if float(charge) > 1.0e-12)
    anion_indices = tuple(int(i) for i, charge in enumerate(z) if float(charge) < -1.0e-12)
    if len(neutral_indices) < 1 or not cation_indices or not anion_indices:
        raise InputError("electrolyte basis requires neutral species, cations, and anions.")

    charged_indices = cation_indices + anion_indices
    pairs = _independent_counterion_pairs(labels, z, x, cation_indices, anion_indices, salt_labels)
    e_matrix = np.zeros((len(charged_indices) - 1, len(charged_indices)), dtype=float)
    charged_pos = {index: pos for pos, index in enumerate(charged_indices)}
    for row, pair in enumerate(pairs):
        e_matrix[row, charged_pos[int(pair["cation"])]] = float(pair["cation_stoich"])
        e_matrix[row, charged_pos[int(pair["anion"])]] = float(pair["anion_stoich"])

    rank = int(np.linalg.matrix_rank(e_matrix, tol=1.0e-12))
    if rank != len(charged_indices) - 1:
        raise InputError("electrolyte counterion-pair matrix is rank deficient.")

    formula_moles = np.asarray(
        [*(x[i] for i in neutral_indices), *(x[int(pair["cation"])] / float(pair["cation_stoich"]) for pair in pairs)],
        dtype=float,
    )
    formula_total = float(np.sum(formula_moles))
    if formula_total <= 0.0:
        raise InputError("electrolyte formula-basis feed has non-positive total.")

    return ElectrolyteBasis(
        species=labels,
        charges=z,
        neutral_indices=neutral_indices,
        charged_indices=tuple(int(i) for i in charged_indices),
        cation_indices=cation_indices,
        anion_indices=anion_indices,
        e_matrix=e_matrix,
        salt_pairs=tuple(pairs),
        formula_feed=formula_moles / formula_total,
        rank=rank,
    )


def _independent_counterion_pairs(
    species: tuple[str, ...],
    charges: np.ndarray,
    feed: np.ndarray,
    cation_indices: tuple[int, ...],
    anion_indices: tuple[int, ...],
    salt_labels: tuple[str, ...],
) -> list[dict[str, Any]]:
    if salt_labels:
        return [_pair_for_salt_label(species, charges, cation_indices, anion_indices, label) for label in salt_labels]

    pairs: list[dict[str, Any]] = []
    if len(anion_indices) == 1:
        anion_i = anion_indices[0]
        for cation_i in cation_indices:
            pairs.append(_pair_payload(species, charges, cation_i, anion_i))
        return pairs
    if len(cation_indices) == 1:
        cation_i = cation_indices[0]
        for anion_i in anion_indices:
            pairs.append(_pair_payload(species, charges, cation_i, anion_i))
        return pairs

    cations = sorted(cation_indices, key=lambda idx: (-float(feed[idx]), int(idx)))
    anions = sorted(anion_indices, key=lambda idx: (-float(feed[idx]), int(idx)))
    anchor_cation = cations[0]
    for anion_i in anions:
        pairs.append(_pair_payload(species, charges, anchor_cation, anion_i))
    anchor_anion = anions[0]
    for cation_i in cations[1:]:
        pairs.append(_pair_payload(species, charges, cation_i, anchor_anion))
    return pairs


def _pair_for_salt_label(
    species: tuple[str, ...],
    charges: np.ndarray,
    cation_indices: tuple[int, ...],
    anion_indices: tuple[int, ...],
    salt_label: str,
) -> dict[str, Any]:
    token = _salt_token(salt_label)
    matches = []
    for cation_i in cation_indices:
        for anion_i in anion_indices:
            candidate = _pair_payload(species, charges, cation_i, anion_i)
            if _salt_token(candidate["label"]) == token:
                matches.append((cation_i, anion_i))
    if len(matches) != 1:
        raise InputError("Could not uniquely map salt '{}' onto independent electrolyte ions.".format(salt_label))
    cation_i, anion_i = matches[0]
    return _pair_payload(species, charges, cation_i, anion_i)


def _pair_payload(species: tuple[str, ...], charges: np.ndarray, cation_i: int, anion_i: int) -> dict[str, Any]:
    cation_charge = abs(float(charges[int(cation_i)]))
    anion_charge = abs(float(charges[int(anion_i)]))
    cation_stoich, anion_stoich = _neutral_stoichiometry(cation_charge, anion_charge)
    cation_label = _ion_stem(species[int(cation_i)], cation_charge)
    anion_label = _ion_stem(species[int(anion_i)], anion_charge)
    cation_suffix = "" if cation_stoich == 1 else str(cation_stoich)
    anion_suffix = "" if anion_stoich == 1 else str(anion_stoich)
    return {
        "label": cation_label + cation_suffix + anion_label + anion_suffix,
        "cation": int(cation_i),
        "anion": int(anion_i),
        "cation_stoich": int(cation_stoich),
        "anion_stoich": int(anion_stoich),
        "cation_charge": float(charges[int(cation_i)]),
        "anion_charge": float(charges[int(anion_i)]),
    }


def _salt_token(label: Any) -> str:
    return "".join(ch for ch in str(label) if ch.isalnum()).lower()


def _neutral_stoichiometry(cation_charge: float, anion_charge: float) -> tuple[int, int]:
    cation_int = int(round(cation_charge))
    anion_int = int(round(anion_charge))
    if cation_int <= 0 or anion_int <= 0:
        raise InputError("electrolyte salt stoichiometry requires non-zero ion charges.")
    if abs(float(cation_int) - cation_charge) > 1.0e-12 or abs(float(anion_int) - anion_charge) > 1.0e-12:
        raise InputError("electrolyte salt stoichiometry currently requires integer ion charges.")
    gcd = int(np.gcd(cation_int, anion_int))
    return anion_int // gcd, cation_int // gcd


def _ion_stem(label: str, charge_magnitude: float | None = None) -> str:
    text = str(label)
    stripped = text.replace("+", "").replace("-", "")
    if charge_magnitude is not None and ("+" in text or "-" in text):
        charge_int = int(round(abs(float(charge_magnitude))))
        suffix = str(charge_int)
        if charge_int > 1 and stripped.endswith(suffix):
            stripped = stripped[: -len(suffix)]
    return stripped
