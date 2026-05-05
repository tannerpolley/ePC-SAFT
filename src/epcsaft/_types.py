"""Small public data types and helpers for the native ePC-SAFT runtime."""

from dataclasses import dataclass, field

import numpy as np


class InputError(Exception):
    """Exception raised for invalid user input."""

    def __init__(self, message):
        self.message = message
        super().__init__(message)


class SolutionError(Exception):
    """Exception raised when a solver fails to converge or returns invalid data."""

    def __init__(self, message, diagnostics=None):
        self.message = message
        self.diagnostics = diagnostics
        if diagnostics is None:
            super().__init__(message)
        else:
            super().__init__(message, diagnostics)


def phase_to_int(phase):
    """Normalize a phase token to the native liquid/vapor integer flag."""
    if phase in (0, "liq", "liquid"):
        return 0
    if phase in (1, "vap", "vapor", "gas"):
        return 1
    raise InputError("phase must be 0/'liq' or 1/'vap'.")


def vector_to_array(values):
    """Convert a vector-like payload into a NumPy float array."""
    return np.asarray(values, dtype=float)


def pair_labels_from_species(params, species):
    """Build ordered mean-ionic pair labels from species metadata."""
    z = np.asarray(params.get("z", []), dtype=float).flatten()
    if z.size == 0 or np.allclose(z, 0.0):
        raise InputError("MIAC calculations require ionic species (non-zero z).")
    if species is None or len(species) != len(z):
        raise InputError("species list (matching x order) is required to label salts.")
    idx_cat = np.where(z > 0)[0]
    idx_an = np.where(z < 0)[0]
    if len(idx_cat) == 0 or len(idx_an) == 0:
        raise InputError("MIAC calculations need at least one cation and one anion.")
    return [species[ic] + species[ia] for ic in idx_cat for ia in idx_an]


def ion_labels_from_species(params, species):
    """Build ordered ion labels from species metadata."""
    z = np.asarray(params.get("z", []), dtype=float).flatten()
    if z.size == 0 or np.allclose(z, 0.0):
        raise InputError("gsolv calculations require ionic species in params['z'].")
    if species is None or len(species) != len(z):
        raise InputError("species list (matching x order) is required to label ions.")
    return [species[i] for i in np.where(np.abs(z) > 1e-12)[0]]


@dataclass(frozen=True, slots=True)
class ActivityCoefficientResult:
    """Bundled activity-coefficient outputs for a single state."""

    species: tuple[str, ...]
    component_activity_coefficients: np.ndarray = field(repr=False)
    solvation_free_energy_values: np.ndarray = field(repr=False)
    mean_ionic_activity_coefficients_mole_fraction_values: np.ndarray = field(repr=False)
    mean_ionic_activity_coefficients_molality_values: np.ndarray = field(repr=False)
    pair_labels: tuple[str, ...] = field(default_factory=tuple)
    ion_labels: tuple[str, ...] = field(default_factory=tuple)
    ion_indices: np.ndarray = field(repr=False, default_factory=lambda: np.asarray([], dtype=int))
    cation_indices: np.ndarray = field(repr=False, default_factory=lambda: np.asarray([], dtype=int))
    anion_indices: np.ndarray = field(repr=False, default_factory=lambda: np.asarray([], dtype=int))
    solvent_indices: np.ndarray = field(repr=False, default_factory=lambda: np.asarray([], dtype=int))
    pair_cation_indices: np.ndarray = field(repr=False, default_factory=lambda: np.asarray([], dtype=int))
    pair_anion_indices: np.ndarray = field(repr=False, default_factory=lambda: np.asarray([], dtype=int))
    pair_nu_cation: np.ndarray = field(repr=False, default_factory=lambda: np.asarray([], dtype=int))
    pair_nu_anion: np.ndarray = field(repr=False, default_factory=lambda: np.asarray([], dtype=int))
    pair_molality: np.ndarray = field(repr=False, default_factory=lambda: np.asarray([], dtype=float))
    pair_conversion_factor: np.ndarray = field(repr=False, default_factory=lambda: np.asarray([], dtype=float))
    solvent_index: int = -1
    osmotic_coefficient: float = float("nan")

    def __post_init__(self):
        object.__setattr__(self, "species", tuple(self.species))
        object.__setattr__(
            self, "component_activity_coefficients", np.asarray(self.component_activity_coefficients, dtype=float)
        )
        object.__setattr__(
            self, "solvation_free_energy_values", np.asarray(self.solvation_free_energy_values, dtype=float)
        )
        object.__setattr__(
            self,
            "mean_ionic_activity_coefficients_mole_fraction_values",
            np.asarray(self.mean_ionic_activity_coefficients_mole_fraction_values, dtype=float),
        )
        object.__setattr__(
            self,
            "mean_ionic_activity_coefficients_molality_values",
            np.asarray(self.mean_ionic_activity_coefficients_molality_values, dtype=float),
        )
        object.__setattr__(self, "pair_labels", tuple(self.pair_labels))
        object.__setattr__(self, "ion_labels", tuple(self.ion_labels))
        object.__setattr__(self, "ion_indices", np.asarray(self.ion_indices, dtype=int))
        object.__setattr__(self, "cation_indices", np.asarray(self.cation_indices, dtype=int))
        object.__setattr__(self, "anion_indices", np.asarray(self.anion_indices, dtype=int))
        object.__setattr__(self, "solvent_indices", np.asarray(self.solvent_indices, dtype=int))
        object.__setattr__(self, "pair_cation_indices", np.asarray(self.pair_cation_indices, dtype=int))
        object.__setattr__(self, "pair_anion_indices", np.asarray(self.pair_anion_indices, dtype=int))
        object.__setattr__(self, "pair_nu_cation", np.asarray(self.pair_nu_cation, dtype=int))
        object.__setattr__(self, "pair_nu_anion", np.asarray(self.pair_nu_anion, dtype=int))
        object.__setattr__(self, "pair_molality", np.asarray(self.pair_molality, dtype=float))
        object.__setattr__(self, "pair_conversion_factor", np.asarray(self.pair_conversion_factor, dtype=float))
        object.__setattr__(self, "solvent_index", int(self.solvent_index))
        object.__setattr__(self, "osmotic_coefficient", float(self.osmotic_coefficient))

    def component_activity_coefficients_map(self):
        """Return component activity coefficients keyed by species label."""
        return {label: float(value) for label, value in zip(self.species, self.component_activity_coefficients)}

    def solvation_free_energy(self):
        """Return ion solvation free-energy values keyed by ion label."""
        out = {}
        for idx, label in zip(self.ion_indices, self.ion_labels):
            out[str(label)] = float(self.solvation_free_energy_values[int(idx)])
        return out

    def mean_ionic_activity_coefficients_map(self, basis="mole"):
        """Return mean-ionic activity coefficients in the requested basis."""
        token = str(basis).strip().lower()
        if token in {"mole", "mole_fraction", "molefraction", "x"}:
            values = self.mean_ionic_activity_coefficients_mole_fraction_values
        elif token in {"molality", "m"}:
            values = self.mean_ionic_activity_coefficients_molality_values
        else:
            raise InputError("basis must be one of: 'mole', 'mole_fraction', 'x', 'molality', 'm'.")
        return {label: float(value) for label, value in zip(self.pair_labels, values)}
