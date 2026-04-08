"""Small public data types and helpers for the native PC-SAFT runtime."""

from dataclasses import dataclass, field

import numpy as np


class InputError(Exception):
    """Exception raised for invalid user input."""

    def __init__(self, message):
        self.message = message
        super().__init__(message)


class SolutionError(Exception):
    """Exception raised when a solver fails to converge or returns invalid data."""

    def __init__(self, message):
        self.message = message
        super().__init__(message)


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
class PhaseResult:
    """Single phase returned by flash and LLE solvers."""

    beta: float
    x: np.ndarray = field(repr=False)
    rho: float
    lnfugcoef: np.ndarray = field(repr=False)
    lnfug: np.ndarray = field(repr=False)

    def __post_init__(self):
        object.__setattr__(self, "beta", float(self.beta))
        object.__setattr__(self, "x", np.asarray(self.x, dtype=float))
        object.__setattr__(self, "rho", float(self.rho))
        object.__setattr__(self, "lnfugcoef", np.asarray(self.lnfugcoef, dtype=float))
        object.__setattr__(self, "lnfug", np.asarray(self.lnfug, dtype=float))


@dataclass(frozen=True, slots=True)
class FlashResult:
    """Two-phase flash result with the scalar target value and phase objects."""

    value: float
    phases: tuple[PhaseResult, ...]
    kind: str

    def __post_init__(self):
        object.__setattr__(self, "value", float(self.value))
        object.__setattr__(self, "phases", tuple(self.phases))
        object.__setattr__(self, "kind", str(self.kind))

    @property
    def pressure(self):
        """Return the phase pressure for `TQ` flash results."""
        return self.value if self.kind == "TQ" else None

    @property
    def temperature(self):
        """Return the phase temperature for `PQ` flash results."""
        return self.value if self.kind == "PQ" else None


@dataclass(frozen=True, slots=True)
class VaporizationResult:
    """Vaporization-pressure result."""

    value: float
    pressure: float

    def __post_init__(self):
        object.__setattr__(self, "value", float(self.value))
        object.__setattr__(self, "pressure", float(self.pressure))
@dataclass(frozen=True, slots=True)
class ActivityCoeffResult:
    """Bundled activity-coefficient outputs for a single state."""

    species: tuple[str, ...]
    component_gamma: np.ndarray = field(repr=False)
    gsolv_values: np.ndarray = field(repr=False)
    gamma_mean_ionic_x_values: np.ndarray = field(repr=False)
    gamma_mean_ionic_m_values: np.ndarray = field(repr=False)
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
    osmotic_c: float = float("nan")

    def __post_init__(self):
        object.__setattr__(self, "species", tuple(self.species))
        object.__setattr__(self, "component_gamma", np.asarray(self.component_gamma, dtype=float))
        object.__setattr__(self, "gsolv_values", np.asarray(self.gsolv_values, dtype=float))
        object.__setattr__(self, "gamma_mean_ionic_x_values", np.asarray(self.gamma_mean_ionic_x_values, dtype=float))
        object.__setattr__(self, "gamma_mean_ionic_m_values", np.asarray(self.gamma_mean_ionic_m_values, dtype=float))
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
        object.__setattr__(self, "osmotic_c", float(self.osmotic_c))

    def component(self):
        """Return component activity values keyed by species label."""
        return {label: float(value) for label, value in zip(self.species, self.component_gamma)}

    def ion(self):
        """Return ion solvation/free-energy values keyed by ion label."""
        out = {}
        for idx, label in zip(self.ion_indices, self.ion_labels):
            out[str(label)] = float(self.gsolv_values[int(idx)])
        return out

    def mean_ionic_x(self):
        """Return mean-ionic activity values on the mole-fraction basis."""
        return {label: float(value) for label, value in zip(self.pair_labels, self.gamma_mean_ionic_x_values)}

    def mean_ionic_m(self):
        """Return mean-ionic activity values on the molality basis."""
        return {label: float(value) for label, value in zip(self.pair_labels, self.gamma_mean_ionic_m_values)}

    def as_basis(self, basis):
        """Return the mean-ionic values in the requested basis."""
        token = str(basis).strip().lower()
        if token in {"x", "mole_fraction", "molefraction"}:
            return self.mean_ionic_x()
        if token in {"m", "molality"}:
            return self.mean_ionic_m()
        raise InputError("basis must be one of: 'mole_fraction', 'x', 'molality', 'm'.")

