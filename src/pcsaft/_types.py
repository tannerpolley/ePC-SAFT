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
    if phase in (0, "liq", "liquid"):
        return 0
    if phase in (1, "vap", "vapor", "gas"):
        return 1
    raise InputError("phase must be 0/'liq' or 1/'vap'.")


def vector_to_array(values):
    return np.asarray(values, dtype=float)


def pair_labels_from_species(params, species):
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
    z = np.asarray(params.get("z", []), dtype=float).flatten()
    if z.size == 0 or np.allclose(z, 0.0):
        raise InputError("gsolv calculations require ionic species in params['z'].")
    if species is None or len(species) != len(z):
        raise InputError("species list (matching x order) is required to label ions.")
    return [species[i] for i in np.where(np.abs(z) > 1e-12)[0]]


@dataclass(frozen=True, slots=True)
class PhaseResult:
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
    value: float
    phases: tuple[PhaseResult, ...]
    kind: str

    def __post_init__(self):
        object.__setattr__(self, "value", float(self.value))
        object.__setattr__(self, "phases", tuple(self.phases))
        object.__setattr__(self, "kind", str(self.kind))

    @property
    def pressure(self):
        return self.value if self.kind == "TQ" else None

    @property
    def temperature(self):
        return self.value if self.kind == "PQ" else None


@dataclass(frozen=True, slots=True)
class VaporizationResult:
    value: float
    pressure: float

    def __post_init__(self):
        object.__setattr__(self, "value", float(self.value))
        object.__setattr__(self, "pressure", float(self.pressure))


@dataclass(frozen=True, slots=True)
class MultiphaseLLEResult:
    n_phases: int
    phases: tuple[PhaseResult, ...]
    tpdf_min: float
    tpdf_seed_x: np.ndarray = field(repr=False)
    converged: bool
    status: int
    message: str
    residual_norm: float
    e_matrix: np.ndarray = field(repr=False)
    ion_pair_rows: tuple = field(default_factory=tuple)
    charged_species: tuple[str, ...] = field(default_factory=tuple)
    charged_species_indices: np.ndarray = field(repr=False, default_factory=lambda: np.asarray([], dtype=int))
    solver_result: object | None = None

    def __post_init__(self):
        object.__setattr__(self, "n_phases", int(self.n_phases))
        object.__setattr__(self, "phases", tuple(self.phases))
        object.__setattr__(self, "tpdf_min", float(self.tpdf_min))
        object.__setattr__(self, "tpdf_seed_x", np.asarray(self.tpdf_seed_x, dtype=float))
        object.__setattr__(self, "converged", bool(self.converged))
        object.__setattr__(self, "status", int(self.status))
        object.__setattr__(self, "message", str(self.message))
        object.__setattr__(self, "residual_norm", float(self.residual_norm))
        object.__setattr__(self, "e_matrix", np.asarray(self.e_matrix, dtype=float))
        object.__setattr__(self, "ion_pair_rows", tuple(self.ion_pair_rows))
        object.__setattr__(self, "charged_species", tuple(self.charged_species))
        object.__setattr__(self, "charged_species_indices", np.asarray(self.charged_species_indices, dtype=int))

    @classmethod
    def from_payload(cls, payload):
        phases = tuple(
            PhaseResult(
                phase.get("beta", 0.0),
                phase.get("x", []),
                phase.get("rho", 0.0),
                phase.get("lnfugcoef", []),
                phase.get("lnfug", []),
            )
            for phase in payload.get("phases", [])
        )
        return cls(
            n_phases=payload.get("n_phases", 0),
            phases=phases,
            tpdf_min=payload.get("tpdf_min", float("nan")),
            tpdf_seed_x=payload.get("tpdf_seed_x", []),
            converged=payload.get("converged", False),
            status=payload.get("status", 0),
            message=payload.get("message", ""),
            residual_norm=payload.get("residual_norm", float("nan")),
            e_matrix=payload.get("e_matrix", []),
            ion_pair_rows=tuple(payload.get("ion_pair_rows", [])),
            charged_species=tuple(payload.get("charged_species", [])),
            charged_species_indices=payload.get("charged_species_indices", []),
            solver_result=payload.get("solver_result"),
        )
