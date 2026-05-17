"""Ascani-style transformed variables for electrolyte liquid splits."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from epcsaft._types import InputError, SolutionError


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


def electrolyte_feed_from_molality_inputs(
    species: Any,
    charges: Any,
    molar_masses: Any,
    *,
    solvent_feed: Any,
    salt_molality: Any,
    basis_mass_kg: float = 1.0,
) -> np.ndarray:
    """Convert mixed-solvent salt molality input into explicit species mole fractions."""
    labels = [str(item) for item in species]
    z = np.asarray(charges, dtype=float).flatten()
    mw = np.asarray(molar_masses, dtype=float).flatten()
    if z.size != len(labels):
        raise InputError("species and charges must have matching lengths for molality feed conversion.")
    if mw.size != len(labels):
        raise InputError("molar_masses must include one value per species.")
    basis_mass = _positive_scalar(basis_mass_kg, "basis_mass_kg")
    solvent_x = normalize_solvent_feed(labels, z, solvent_feed)
    salt_items = normalize_salt_molality(salt_molality)
    neutral_mw = float(np.sum(solvent_x * mw))
    if not np.isfinite(neutral_mw) or neutral_mw <= 0.0:
        raise InputError("solvent_feed produced an invalid salt-free solvent molecular weight.")
    moles = np.zeros(len(labels), dtype=float)
    moles += solvent_x * (basis_mass / neutral_mw)
    cation_indices = tuple(int(i) for i, charge in enumerate(z) if float(charge) > 1.0e-12)
    anion_indices = tuple(int(i) for i, charge in enumerate(z) if float(charge) < -1.0e-12)
    for salt_label, molality in salt_items.items():
        pair = _pair_for_salt_label(tuple(labels), z, cation_indices, anion_indices, salt_label)
        salt_moles = float(molality) * basis_mass
        moles[int(pair["cation"])] += salt_moles * float(pair["cation_stoich"])
        moles[int(pair["anion"])] += salt_moles * float(pair["anion_stoich"])
    total = float(np.sum(moles))
    if total <= 0.0:
        raise InputError("Computed electrolyte feed has non-positive total moles.")
    feed = moles / total
    charge = float(np.dot(feed, z))
    if abs(charge) > 1.0e-10:
        raise InputError(f"molality-derived electrolyte feed must be charge neutral; charge balance is {charge}.")
    return feed


def normalize_solvent_feed(species: list[str], charges: np.ndarray, solvent_feed: Any) -> np.ndarray:
    neutral_indices = [i for i, charge in enumerate(charges) if abs(float(charge)) <= 1.0e-12]
    if not neutral_indices:
        raise InputError("solvent_feed requires at least one neutral solvent species.")
    solvent_x = np.zeros(len(species), dtype=float)
    if isinstance(solvent_feed, dict):
        for label, value in solvent_feed.items():
            try:
                index = species.index(str(label))
            except ValueError as exc:
                raise InputError(f"solvent_feed species '{label}' is not present in the mixture.") from exc
            if index not in neutral_indices:
                raise InputError(f"solvent_feed species '{label}' is ionic; expected neutral solvents only.")
            solvent_x[index] = float(value)
    else:
        values = np.asarray(solvent_feed, dtype=float).flatten()
        if values.size != len(neutral_indices):
            raise InputError("solvent_feed must be a dict or a vector with one entry per neutral solvent species.")
        for index, value in zip(neutral_indices, values):
            solvent_x[index] = float(value)
    if not np.all(np.isfinite(solvent_x)) or np.any(solvent_x < 0.0):
        raise InputError("solvent_feed must contain finite non-negative values.")
    total = float(np.sum(solvent_x))
    if total <= 0.0:
        raise InputError("solvent_feed must have a positive sum.")
    return solvent_x / total


def normalize_salt_molality(salt_molality: Any) -> dict[str, float]:
    if not isinstance(salt_molality, dict) or not salt_molality:
        raise InputError("salt_molality must be a non-empty dict like {'NaCl': 1.0}.")
    out: dict[str, float] = {}
    for label, value in salt_molality.items():
        molality = float(value)
        if not np.isfinite(molality) or molality < 0.0:
            raise InputError("salt_molality values must be finite and non-negative.")
        out[str(label)] = molality
    return out


def electrolyte_formula_basis(
    species: Any,
    charges: Any,
    feed: Any,
    *,
    salt_labels: tuple[str, ...] = (),
) -> dict[str, Any]:
    z = np.asarray(charges, dtype=float).flatten()
    labels = tuple(str(item) for item in species)
    neutral_indices = tuple(int(i) for i, charge in enumerate(z) if abs(float(charge)) <= 1.0e-12)
    cation_indices = tuple(int(i) for i, charge in enumerate(z) if float(charge) > 1.0e-12)
    anion_indices = tuple(int(i) for i, charge in enumerate(z) if float(charge) < -1.0e-12)
    if len(neutral_indices) < 2:
        raise InputError("electrolyte_lle requires at least two neutral solvent species.")
    if not cation_indices or not anion_indices:
        raise InputError("electrolyte_lle requires at least one cation and one anion.")
    payload = build_electrolyte_basis(labels, z, feed, salt_labels=salt_labels).to_dict()
    return {
        "neutral_indices": tuple(payload["neutral_indices"]),
        "charged_indices": tuple(payload["charged_indices"]),
        "cation_indices": tuple(payload["cation_indices"]),
        "anion_indices": tuple(payload["anion_indices"]),
        "salt_pairs": tuple(payload["salt_pairs"]),
        "formula_feed": np.asarray(payload["formula_feed"], dtype=float),
        "e_matrix": np.asarray(payload["e_matrix"], dtype=float),
        "basis_rank": int(payload["basis_rank"]),
        "variable_model": str(payload["variable_model"]),
    }


def formula_to_explicit_composition(
    formula_composition: Any,
    basis: dict[str, Any],
    ncomp: int,
) -> tuple[np.ndarray, float]:
    formula = np.asarray(formula_composition, dtype=float)
    explicit = np.zeros(int(ncomp), dtype=float)
    neutral_indices = basis["neutral_indices"]
    salt_pairs = basis["salt_pairs"]
    for pos, index in enumerate(neutral_indices):
        explicit[int(index)] += float(formula[pos])
    offset = len(neutral_indices)
    for salt_pos, pair in enumerate(salt_pairs):
        amount = float(formula[offset + salt_pos])
        explicit[int(pair["cation"])] += amount * float(pair.get("cation_stoich", 1.0))
        explicit[int(pair["anion"])] += amount * float(pair.get("anion_stoich", 1.0))
    total = float(np.sum(explicit))
    if total <= 0.0:
        raise SolutionError("Formula-basis electrolyte phase expanded to a non-positive explicit composition.")
    return explicit / total, total


def explicit_to_formula_composition(composition: Any, basis: dict[str, Any]) -> np.ndarray:
    comp = np.asarray(composition, dtype=float)
    values = [float(comp[index]) for index in basis["neutral_indices"]]
    values.extend(
        float(comp[int(pair["cation"])]) / float(pair.get("cation_stoich", 1.0)) for pair in basis["salt_pairs"]
    )
    out = np.asarray(values, dtype=float)
    total = float(np.sum(out))
    if total <= 0.0:
        raise InputError("Explicit electrolyte composition cannot be represented on the formula basis.")
    return out / total


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
        raise InputError(f"Could not uniquely map salt '{salt_label}' onto independent electrolyte ions.")
    cation_i, anion_i = matches[0]
    return _pair_payload(species, charges, cation_i, anion_i)


def _positive_scalar(value: Any, label: str) -> float:
    if value is None:
        raise InputError(f"{label} is required.")
    out = float(value)
    if not np.isfinite(out) or out <= 0.0:
        raise InputError(f"{label} must be a positive finite scalar.")
    return out


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
    cation_int = round(cation_charge)
    anion_int = round(anion_charge)
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
        charge_int = round(abs(float(charge_magnitude)))
        suffix = str(charge_int)
        if charge_int > 1 and stripped.endswith(suffix):
            stripped = stripped[: -len(suffix)]
    return stripped
