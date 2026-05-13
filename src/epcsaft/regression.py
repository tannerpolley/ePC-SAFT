"""Parameter regression helpers for user-owned ePC-SAFT datasets."""

from __future__ import annotations

import csv
import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from ._types import InputError, phase_to_int
from .epcsaft import (
    _core,
    _evaluate_generic_native_debug,
    _fit_generic_native_ceres,
    _fit_generic_native_least_squares,
    _fit_pure_neutral_native_debug,
    _fit_pure_neutral_native_least_squares,
    check_association,
    create_struct,
    vector_to_array,
)
from .parameter_templates import _infer_pure_template_name
from .parameters import (
    _MISSING,
    _deterministic_default,
    _invalidate_dataset_cache,
    _load_dataset,
    _matrix_value,
    _normalize_component,
    _resolve_component_field_with_source,
    _solvent_token_for_component,
    get_prop_dict,
    molality_to_molefraction,
)

PURE_NEUTRAL_MODE = "pure_neutral"
PURE_ION_MODE = "pure_ion"
BINARY_PAIR_MODE = "binary_pair"
LIQUID_ELECTROLYTE_MODE = "liquid_electrolyte"

TERM_DENSITY = "density"
TERM_PURE_VLE = "pure_vle_fugacity_balance"
TERM_OSMOTIC = "osmotic_coefficient"
TERM_MIAC = "mean_ionic_activity"
TERM_BINARY_VLE = "binary_vle_fugacity_balance"
TERM_BINARY_LLE = "binary_lle_fugacity_balance"
TERM_RELATIVE_PERMITTIVITY = "relative_permittivity"
TERM_MEA_CO2_H2O_DENSITY = "mea_co2_h2o_density"
TERM_MEA_CO2_H2O_CO2_FUGACITY = "mea_co2_h2o_co2_fugacity"
TERM_MEA_CO2_H2O_OSMOTIC = "mea_co2_h2o_osmotic_coefficient"

TARGET_ROW_FAMILY_ALIASES = {
    "density": "pure_density",
    "pure_density": "pure_density",
    "pure_vapor_pressure": "pure_vapor_pressure",
    "pure_vle": "pure_vapor_pressure",
    "p_rho_t": "p_rho_t",
    "prt": "p_rho_t",
    "binary_vle": "binary_vle",
    "binary_lle": "binary_lle",
    "osmotic": TERM_OSMOTIC,
    "osmotic_coefficient": TERM_OSMOTIC,
    "miac": TERM_MIAC,
    "mean_ionic_activity": TERM_MIAC,
    "relative_permittivity": TERM_RELATIVE_PERMITTIVITY,
    "dielectric": TERM_RELATIVE_PERMITTIVITY,
    "activity": "activity",
    "fugacity": "fugacity",
    "speciation": "speciation",
    "vle_partial_pressure": "vle_partial_pressure",
    "partial_pressure": "vle_partial_pressure",
    "lle_phase_composition": "lle_phase_composition",
    "regularization": "regularization",
}
TARGET_ROW_FAMILIES = tuple(dict.fromkeys(TARGET_ROW_FAMILY_ALIASES.values()))

PURE_DENSITY_KEYS_MOLAR = ("rho",)
PURE_DENSITY_KEYS_MASS = ("rho_kg_m3", "rho_mass_kg_m3", "rho_liq_kg_m3", "rho_sat_liq_kg_m3")

PURE_REQUIRED_FIELDS = (
    "m",
    "s",
    "e",
    "e_assoc",
    "vol_a",
    "z",
    "dielc",
    "d_born",
    "f_solv",
    "MW",
)

MATRIX_FILE_NAMES = {
    "k_ij": "k_ij.csv",
    "l_ij": "l_ij.csv",
    "k_hb_ij": "k_hb_ij.csv",
}
DEFAULT_TARGETS = {
    PURE_NEUTRAL_MODE: {
        "nonassociating": ("m", "s", "e"),
        "associating": ("m", "s", "e", "e_assoc", "vol_a"),
    },
    PURE_ION_MODE: ("s", "e"),
    BINARY_PAIR_MODE: ("k_ij",),
}

DEFAULT_BOUNDS = {
    "m": (0.1, 25.0),
    "s": (1.0, 10.0),
    "e": (1.0, 12000.0),
    "e_assoc": (0.0, 20000.0),
    "vol_a": (0.0, 1.0),
    "d_born": (0.1, 10.0),
    "k_ij": (-2.0, 2.0),
    "l_ij": (-2.0, 2.0),
    "k_hb_ij": (-2.0, 2.0),
    "k_ij_slope": (-1.0, 1.0),
    "k_ij_intercept": (-2.0, 2.0),
    "l_ij_slope": (-1.0, 1.0),
    "l_ij_intercept": (-2.0, 2.0),
    "k_hb_ij_slope": (-1.0, 1.0),
    "k_hb_ij_intercept": (-2.0, 2.0),
}

NATIVE_TARGET_KINDS = {
    "m": 0,
    "s": 1,
    "e": 2,
    "e_assoc": 3,
    "vol_a": 4,
    "d_born": 5,
    "k_ij": 6,
    "l_ij": 7,
    "k_hb_ij": 8,
}

NATIVE_TERM_KINDS = {
    TERM_DENSITY: 1,
    TERM_PURE_VLE: 2,
    TERM_OSMOTIC: 3,
    TERM_MIAC: 4,
    TERM_BINARY_VLE: 5,
    TERM_MEA_CO2_H2O_DENSITY: 1,
    TERM_MEA_CO2_H2O_OSMOTIC: 3,
    TERM_MEA_CO2_H2O_CO2_FUGACITY: 6,
}

HESSIAN_NOT_IMPLEMENTED_REASON = "Hessian support is a skeleton for future IPOPT-compatible optimizer integration."
LIQUID_ELECTROLYTE_BACKEND_UNAVAILABLE_REASON = (
    "backend_unavailable: liquid-electrolyte parameter fitting needs native optimizer internals that are not part "
    "of this API-contract slice."
)
BINARY_CERES_BACKEND_UNAVAILABLE_REASON = (
    "backend_unavailable: binary Ceres regression requires real k_ij residual derivatives through fugacity "
    "and density-root sensitivities; the legacy generic binary optimizer is finite-difference only."
)


def _copy_mapping(mapping: Mapping[str, Any] | None) -> dict[str, Any]:
    return {} if mapping is None else {str(k): v for k, v in mapping.items()}


def _copy_bounds_contract(
    bounds: Mapping[str, tuple[float | None, float | None]] | FitBounds | None,
) -> dict[str, tuple[float | None, float | None]]:
    if bounds is None:
        return {}
    if isinstance(bounds, FitBounds):
        names = sorted(set(bounds.lower) | set(bounds.upper))
        return {name: (bounds.lower.get(name), bounds.upper.get(name)) for name in names}
    return {
        str(name): (None if pair[0] is None else float(pair[0]), None if pair[1] is None else float(pair[1]))
        for name, pair in bounds.items()
    }


_DBORN_DIRECT_SOURCES = {
    "dielectric",
    "dielectric_measurement",
    "dielectric_or_ion_activity",
    "explicit_override",
    "ion_activity",
    "mean_ionic_activity",
    "osmotic",
    "relative_permittivity",
}
_DIRECT_BINARY_SOURCES = {
    "direct_binary_activity",
    "direct_binary_lle",
    "direct_binary_vle",
    "explicit_override",
}
_DIRECT_ELECTROLYTE_BINARY_SOURCES = {
    "direct_binary_activity",
    "direct_electrolyte_activity",
    "direct_electrolyte_osmotic",
    "direct_salt_pair_data",
    "explicit_override",
}
_DIRECT_NEUTRAL_ION_SOURCES = {
    "direct_neutral_ion_activity",
    "direct_neutral_ion_data",
    "direct_electrolyte_activity",
    "explicit_override",
}


def _source_token(value: str | None) -> str:
    return str(value or "").strip().lower()


@dataclass(frozen=True, slots=True)
class FitParameter:
    """One pure/species parameter requested for provenance-aware fitting."""

    component: str
    parameter: str
    source: str = ""
    allow_without_direct_data: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(self, "component", _normalize_component(self.component))
        object.__setattr__(self, "parameter", str(self.parameter))
        object.__setattr__(self, "source", _source_token(self.source))
        object.__setattr__(self, "allow_without_direct_data", bool(self.allow_without_direct_data))

    @property
    def key(self) -> str:
        """Return the stable report key for this parameter."""
        return f"{self.component}.{self.parameter}"


@dataclass(frozen=True, slots=True)
class BinaryInteraction:
    """One binary-interaction parameter requested for provenance-aware fitting."""

    pair: tuple[str, str]
    parameter: str = "k_ij"
    source: str = ""
    allow_without_direct_data: bool = False

    def __post_init__(self) -> None:
        pair = tuple(_normalize_component(str(name)) for name in self.pair)
        if len(pair) != 2 or pair[0] == pair[1]:
            raise InputError("BinaryInteraction requires two distinct pair components.")
        object.__setattr__(self, "pair", pair)
        object.__setattr__(self, "parameter", str(self.parameter))
        object.__setattr__(self, "source", _source_token(self.source))
        object.__setattr__(self, "allow_without_direct_data", bool(self.allow_without_direct_data))

    @property
    def key(self) -> str:
        """Return the stable report key for this interaction."""
        return f"{self.pair[0]}:{self.pair[1]}.{self.parameter}"


@dataclass(frozen=True, slots=True)
class RelativePermittivityResidual:
    """One relative-permittivity data residual for regression/provenance workflows."""

    T: float
    P: float
    composition: Mapping[str, float]
    epsilon_r_exp: float
    weight: float = 1.0
    source: str = "dielectric_measurement"

    def __post_init__(self) -> None:
        composition = {str(k): float(v) for k, v in self.composition.items()}
        if not composition:
            raise InputError("RelativePermittivityResidual requires a non-empty composition mapping.")
        total = float(sum(composition.values()))
        if total <= 0.0 or not math.isfinite(total):
            raise InputError("RelativePermittivityResidual composition must have a positive finite sum.")
        normalized = {name: value / total for name, value in composition.items()}
        if any(value < -1.0e-12 or not math.isfinite(value) for value in normalized.values()):
            raise InputError("RelativePermittivityResidual composition values must be finite and non-negative.")
        object.__setattr__(self, "T", float(self.T))
        object.__setattr__(self, "P", float(self.P))
        object.__setattr__(self, "composition", normalized)
        object.__setattr__(self, "epsilon_r_exp", float(self.epsilon_r_exp))
        object.__setattr__(self, "weight", float(self.weight))
        object.__setattr__(self, "source", _source_token(self.source))
        if not math.isfinite(self.T) or self.T <= 0.0:
            raise InputError("RelativePermittivityResidual T must be positive and finite.")
        if not math.isfinite(self.P) or self.P <= 0.0:
            raise InputError("RelativePermittivityResidual P must be positive and finite.")
        if not math.isfinite(self.epsilon_r_exp) or self.epsilon_r_exp <= 0.0:
            raise InputError("RelativePermittivityResidual epsilon_r_exp must be positive and finite.")
        if not math.isfinite(self.weight) or self.weight <= 0.0:
            raise InputError("RelativePermittivityResidual weight must be positive and finite.")

    def to_record(self, *, species: Sequence[str] | None = None) -> dict[str, Any]:
        """Return a flat regression record with ``x_*`` composition columns."""
        labels = [str(label) for label in (species or self.composition.keys())]
        record: dict[str, Any] = {
            "T": self.T,
            "P": self.P,
            "epsilon_r_exp": self.epsilon_r_exp,
            "source": self.source,
        }
        for label in labels:
            record[f"x_{label}"] = float(self.composition.get(label, 0.0))
        return record

    def to_fit_term(self, *, species: Sequence[str] | None = None) -> FitTerm:
        """Return this residual as a first-class regression term descriptor."""
        return FitTerm(
            term_type=TERM_RELATIVE_PERMITTIVITY,
            records=(self.to_record(species=species),),
            weight=self.weight,
            residual_count=1,
        )


def validate_regression_provenance(
    parameters: Sequence[FitParameter | BinaryInteraction],
    *,
    terms: Sequence[FitTerm] | None = None,
    species: Sequence[str] | None = None,
    charges: Sequence[float] | Mapping[str, float] | None = None,
    strict: bool = True,
) -> dict[str, Any]:
    """Validate provenance declarations for fitted parameters.

    The validator is chemistry-agnostic: callers declare what data supports a
    parameter, and the package rejects unsupported electrostatic/Born/binary
    targets unless an explicit override is supplied.
    """
    _ = terms
    labels = [] if species is None else [_normalize_component(label) for label in species]
    charge_map = _charge_map(labels, charges)
    fitted: list[str] = []
    parameter_sources: dict[str, str] = {}
    data_sources: dict[str, list[str]] = {}
    warnings: list[str] = []
    errors: list[str] = []

    for item in parameters:
        if isinstance(item, FitParameter):
            fitted.append(item.key)
            parameter_sources[item.key] = item.source
            data_sources[item.key] = [item.source] if item.source else []
            if item.parameter == "d_born" and item.source not in _DBORN_DIRECT_SOURCES:
                message = (
                    f"d_born for {item.component} requires dielectric, relative-permittivity, ion-activity, "
                    "osmotic, or explicit_override provenance before fitting."
                )
                if item.allow_without_direct_data:
                    warnings.append(message)
                else:
                    errors.append(message)
            continue

        if not isinstance(item, BinaryInteraction):
            raise InputError("parameters must contain FitParameter or BinaryInteraction declarations.")
        fitted.append(item.key)
        parameter_sources[item.key] = item.source
        data_sources[item.key] = [item.source] if item.source else []
        left_charge = charge_map.get(item.pair[0])
        right_charge = charge_map.get(item.pair[1])
        if item.parameter in {"k_ij", "l_ij", "k_hb_ij"}:
            if left_charge is not None and right_charge is not None:
                if left_charge * right_charge > 0.0 and item.parameter == "k_ij":
                    errors.append(
                        f"{item.parameter} for same-sign ionic pair {item.pair[0]}:{item.pair[1]} is not a "
                        "meaningful fit target because same-sign short-range dispersion is suppressed."
                    )
                elif left_charge * right_charge < 0.0 and item.source not in _DIRECT_ELECTROLYTE_BINARY_SOURCES:
                    message = (
                        f"{item.parameter} for opposite-sign ionic pair {item.pair[0]}:{item.pair[1]} requires "
                        "direct electrolyte activity/osmotic/salt-pair provenance or explicit_override."
                    )
                    if item.allow_without_direct_data:
                        warnings.append(message)
                    else:
                        errors.append(message)
                elif (left_charge == 0.0) != (right_charge == 0.0) and item.source not in _DIRECT_NEUTRAL_ION_SOURCES:
                    message = (
                        f"{item.parameter} for neutral-ion pair {item.pair[0]}:{item.pair[1]} requires direct "
                        "neutral-ion or electrolyte provenance or explicit_override."
                    )
                    if item.allow_without_direct_data:
                        warnings.append(message)
                    else:
                        errors.append(message)
                elif left_charge == 0.0 and right_charge == 0.0 and item.source not in _DIRECT_BINARY_SOURCES:
                    message = f"{item.parameter} for neutral pair {item.pair[0]}:{item.pair[1]} requires direct binary provenance."
                    if item.allow_without_direct_data:
                        warnings.append(message)
                    else:
                        errors.append(message)
            elif item.source not in _DIRECT_BINARY_SOURCES and not item.allow_without_direct_data:
                errors.append(f"{item.parameter} for {item.pair[0]}:{item.pair[1]} requires declared data provenance.")

    if errors and strict:
        raise InputError("; ".join(errors))
    if errors:
        warnings.extend(errors)
    return {
        "fitted_parameters": fitted,
        "fixed_parameters": [],
        "parameter_sources": parameter_sources,
        "data_sources_by_parameter": data_sources,
        "warnings": warnings,
        "backend_policy": "native_regression_only",
    }


def _charge_map(species: Sequence[str], charges: Sequence[float] | Mapping[str, float] | None) -> dict[str, float]:
    if charges is None:
        return {}
    if isinstance(charges, Mapping):
        return {_normalize_component(str(label)): float(value) for label, value in charges.items()}
    values = np.asarray(charges, dtype=float).flatten()
    if values.size != len(species):
        raise InputError("charges must match species length when provided as a sequence.")
    return {label: float(value) for label, value in zip(species, values)}


@dataclass(slots=True)
class FitBounds:
    """Bounds for named regression variables."""

    lower: dict[str, float] = field(default_factory=dict)
    upper: dict[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.lower = {str(k): float(v) for k, v in self.lower.items()}
        self.upper = {str(k): float(v) for k, v in self.upper.items()}

    def arrays_for(self, names: Sequence[str]) -> tuple[np.ndarray, np.ndarray]:
        lower: list[float] = []
        upper: list[float] = []
        for name in names:
            lo, hi = DEFAULT_BOUNDS[name]
            lower.append(self.lower.get(name, lo))
            upper.append(self.upper.get(name, hi))
        return np.asarray(lower, dtype=float), np.asarray(upper, dtype=float)


@dataclass(slots=True)
class FitTerm:
    """One weighted family of regression residuals."""

    term_type: str
    records: tuple[dict[str, Any], ...] = field(default_factory=tuple, repr=False)
    weight: float = 1.0
    residual_count: int = 0

    def __post_init__(self) -> None:
        self.term_type = str(self.term_type)
        self.records = tuple(dict(record) for record in self.records)
        self.weight = float(self.weight)
        self.residual_count = int(self.residual_count)


def _normalize_target_row_family(row_family: str) -> str:
    token = str(row_family).strip().lower().replace("-", "_").replace(" ", "_")
    if token not in TARGET_ROW_FAMILY_ALIASES:
        allowed = ", ".join(TARGET_ROW_FAMILIES)
        raise InputError(f"Unsupported target row family '{row_family}'. Supported families: {allowed}.")
    return TARGET_ROW_FAMILY_ALIASES[token]


def _record_has_any(record: Mapping[str, Any], keys: Sequence[str]) -> bool:
    return any(key in record and record[key] not in (None, "") for key in keys)


def _record_has_prefixed(record: Mapping[str, Any], prefixes: Sequence[str]) -> bool:
    return any(any(str(key).startswith(prefix) for key in record) for prefix in prefixes)


def _require_target_row_fields(row_family: str, values: Mapping[str, Any]) -> None:
    temperature_families = {
        "pure_density",
        "pure_vapor_pressure",
        "p_rho_t",
        "binary_vle",
        "binary_lle",
        TERM_OSMOTIC,
        TERM_MIAC,
        TERM_RELATIVE_PERMITTIVITY,
        "activity",
        "fugacity",
        "speciation",
        "vle_partial_pressure",
        "lle_phase_composition",
    }
    if row_family == "regularization":
        if not _record_has_any(values, ("parameter", "target_parameter", "name")):
            raise InputError("regularization target rows require a parameter or target_parameter field.")
        if not _record_has_any(values, ("target_value", "value", "prior", "center")):
            raise InputError("regularization target rows require a target_value, value, prior, or center field.")
        return

    if row_family in temperature_families and not _record_has_any(values, ("T", "temperature", "temperature_K")):
        raise InputError(f"{row_family} target rows require T or temperature_K.")

    if row_family == "pure_density" and not _record_has_any(
        values, (*PURE_DENSITY_KEYS_MOLAR, *PURE_DENSITY_KEYS_MASS)
    ):
        raise InputError("pure_density target rows require rho or a mass-density value.")
    if row_family == "pure_vapor_pressure" and not _record_has_any(values, ("P", "vapor_pressure", "P_sat")):
        raise InputError("pure_vapor_pressure target rows require P, vapor_pressure, or P_sat.")
    if row_family == "p_rho_t" and (
        not _record_has_any(values, ("P", "pressure", "pressure_Pa"))
        or not _record_has_any(values, (*PURE_DENSITY_KEYS_MOLAR, *PURE_DENSITY_KEYS_MASS))
    ):
        raise InputError("p_rho_t target rows require pressure and density values.")
    if row_family == "binary_vle" and not (
        _record_has_prefixed(values, ("x_",)) and _record_has_prefixed(values, ("y_",))
    ):
        raise InputError("binary_vle target rows require x_* and y_* composition columns.")
    if row_family == "binary_lle" and not (
        (_record_has_prefixed(values, ("x_alpha_",)) and _record_has_prefixed(values, ("x_beta_",)))
        or (_record_has_prefixed(values, ("x_phase1_",)) and _record_has_prefixed(values, ("x_phase2_",)))
    ):
        raise InputError("binary_lle target rows require alpha/beta or phase1/phase2 composition columns.")
    if row_family == TERM_OSMOTIC and not _record_has_any(values, ("osmotic_coefficient", "osmotic")):
        raise InputError("osmotic_coefficient target rows require osmotic_coefficient or osmotic.")
    if row_family == TERM_MIAC and not _record_has_any(
        values, ("mean_ionic_activity", "mean_ionic_activity_coefficient", "miac")
    ):
        raise InputError("mean_ionic_activity target rows require a MIAC value.")
    if row_family == TERM_RELATIVE_PERMITTIVITY and not _record_has_any(
        values, ("epsilon_r_exp", "relative_permittivity", "epsilon_r")
    ):
        raise InputError("relative_permittivity target rows require epsilon_r_exp or relative_permittivity.")
    if row_family == "activity" and not (
        _record_has_prefixed(values, ("activity_",)) or _record_has_any(values, ("target_activities", "activities"))
    ):
        raise InputError("activity target rows require activity_* columns or target_activities.")
    if row_family == "fugacity" and not (
        _record_has_prefixed(values, ("fugacity_",)) or _record_has_any(values, ("target_fugacities", "fugacities"))
    ):
        raise InputError("fugacity target rows require fugacity_* columns or target_fugacities.")
    if row_family == "speciation" and not (
        _record_has_prefixed(values, ("x_",))
        or _record_has_any(values, ("target_x", "target_composition", "speciation_targets"))
    ):
        raise InputError("speciation target rows require x_* columns, target_x, or speciation_targets.")
    if row_family == "vle_partial_pressure" and not _record_has_any(
        values, ("target_partial_pressures", "partial_pressures", "target_pressure", "target_partial_pressure")
    ):
        raise InputError("vle_partial_pressure target rows require partial-pressure targets.")
    if row_family == "lle_phase_composition" and not (
        _record_has_prefixed(values, ("x_alpha_", "x_beta_", "x_phase1_", "x_phase2_"))
        or _record_has_any(values, ("target_x_alpha", "target_x_beta", "phase_compositions"))
    ):
        raise InputError("lle_phase_composition target rows require phase-composition targets.")


@dataclass(slots=True)
class TargetRow:
    """Application-neutral regression or validation target row."""

    row_family: str
    values: dict[str, Any] = field(default_factory=dict, repr=False)
    row_id: str = ""
    weight: float = 1.0
    phase: str = ""
    species: tuple[str, ...] = field(default_factory=tuple)
    source: str = ""
    metadata: dict[str, Any] = field(default_factory=dict, repr=False)

    def __post_init__(self) -> None:
        self.row_family = _normalize_target_row_family(self.row_family)
        self.values = dict(self.values)
        self.row_id = str(self.row_id or self.values.get("row_id", self.values.get("id", "")))
        self.weight = float(self.weight)
        self.phase = str(self.phase or self.values.get("phase", ""))
        self.species = tuple(str(label) for label in self.species)
        self.source = str(self.source or self.values.get("source", ""))
        self.metadata = _copy_mapping(self.metadata)
        if not math.isfinite(self.weight) or self.weight <= 0.0:
            raise InputError("TargetRow weight must be positive and finite.")
        _require_target_row_fields(self.row_family, self.values)

    @classmethod
    def from_record(cls, record: Mapping[str, Any], *, default_family: str | None = None) -> TargetRow:
        payload = dict(record)
        family = payload.pop("row_family", payload.pop("target_family", payload.pop("family", default_family)))
        if family is None:
            raise InputError("Target rows require row_family, target_family, family, or a default_family.")
        row_id = str(payload.pop("row_id", payload.pop("id", "")))
        weight = float(payload.pop("weight", 1.0))
        phase = str(payload.pop("phase", ""))
        source = str(payload.pop("source", ""))
        metadata = payload.pop("metadata", {})
        species = payload.pop("species", ())
        return cls(
            str(family),
            payload,
            row_id=row_id,
            weight=weight,
            phase=phase,
            species=tuple(species),
            source=source,
            metadata=metadata,
        )

    def to_record(self) -> dict[str, Any]:
        record = {"row_family": self.row_family, **dict(self.values)}
        if self.row_id:
            record["row_id"] = self.row_id
        if self.weight != 1.0:
            record["weight"] = self.weight
        if self.phase:
            record["phase"] = self.phase
        if self.species:
            record["species"] = list(self.species)
        if self.source:
            record["source"] = self.source
        if self.metadata:
            record["metadata"] = _json_like_regression(self.metadata)
        return record


@dataclass(slots=True)
class TargetDataset:
    """Collection of generic target rows for regression and validation workflows."""

    rows: tuple[TargetRow, ...] = field(default_factory=tuple)
    name: str = ""
    species: tuple[str, ...] = field(default_factory=tuple)
    parameter_set: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict, repr=False)

    def __post_init__(self) -> None:
        self.rows = tuple(row if isinstance(row, TargetRow) else TargetRow.from_record(row) for row in self.rows)
        if not self.rows:
            raise InputError("TargetDataset requires at least one TargetRow.")
        self.name = str(self.name)
        self.species = tuple(str(label) for label in self.species)
        self.parameter_set = None if self.parameter_set is None else str(self.parameter_set)
        self.metadata = _copy_mapping(self.metadata)

    @classmethod
    def from_records(
        cls,
        records: Any,
        *,
        default_family: str | None = None,
        name: str = "",
        species: Sequence[str] = (),
        parameter_set: str | Path | None = None,
        metadata: Mapping[str, Any] | None = None,
    ) -> TargetDataset:
        rows = tuple(
            TargetRow.from_record(record, default_family=default_family) for record in load_regression_records(records)
        )
        return cls(
            rows=rows,
            name=name,
            species=tuple(str(label) for label in species),
            parameter_set=_source_dataset_label(parameter_set),
            metadata=_copy_mapping(metadata),
        )

    @property
    def families(self) -> tuple[str, ...]:
        return tuple(dict.fromkeys(row.row_family for row in self.rows))

    def rows_for(self, row_family: str) -> tuple[TargetRow, ...]:
        normalized = _normalize_target_row_family(row_family)
        return tuple(row for row in self.rows if row.row_family == normalized)

    def to_records(self) -> list[dict[str, Any]]:
        return [row.to_record() for row in self.rows]


@dataclass(slots=True)
class FitProblem:
    """Normalized description of a regression problem."""

    mode: str
    records: tuple[dict[str, Any], ...] = field(default_factory=tuple, repr=False)
    component: str | None = None
    pair: tuple[str, str] | None = None
    solvent: str | None = None
    dataset: str | None = None
    fit_targets: tuple[str, ...] = field(default_factory=tuple)
    optimization_parameters: tuple[str, ...] = field(default_factory=tuple)
    fixed_parameters: dict[str, Any] = field(default_factory=dict, repr=False)
    initial_guess: dict[str, float] = field(default_factory=dict, repr=False)
    bounds: dict[str, tuple[float | None, float | None]] = field(default_factory=dict)
    weights: dict[str, float] = field(default_factory=dict)
    loss: str = "linear"
    solver_options: dict[str, Any] = field(default_factory=dict, repr=False)
    output_report: bool = False
    assoc_scheme: str = ""
    temperature_model: str = "constant"
    terms: tuple[FitTerm, ...] = field(default_factory=tuple)
    target_dataset: TargetDataset | None = None
    pure_file_hint: str | None = None

    def __post_init__(self) -> None:
        self.mode = str(self.mode)
        self.records = tuple(dict(record) for record in self.records)
        self.component = None if self.component is None else str(self.component)
        self.pair = None if self.pair is None else (str(self.pair[0]), str(self.pair[1]))
        self.solvent = None if self.solvent is None else str(self.solvent)
        self.dataset = None if self.dataset is None else str(self.dataset)
        self.fit_targets = tuple(str(name) for name in self.fit_targets)
        self.optimization_parameters = tuple(str(name) for name in self.optimization_parameters)
        self.fixed_parameters = _copy_mapping(self.fixed_parameters)
        self.initial_guess = {str(k): float(v) for k, v in self.initial_guess.items()}
        self.bounds = _copy_bounds_contract(self.bounds)
        self.weights = {str(k): float(v) for k, v in self.weights.items()}
        self.loss = str(self.loss)
        self.solver_options = _copy_mapping(self.solver_options)
        self.output_report = bool(self.output_report)
        self.assoc_scheme = str(self.assoc_scheme or "")
        self.temperature_model = str(self.temperature_model)
        self.terms = tuple(self.terms)
        if self.target_dataset is not None and not isinstance(self.target_dataset, TargetDataset):
            self.target_dataset = TargetDataset.from_records(self.target_dataset)
        self.pure_file_hint = None if self.pure_file_hint is None else str(self.pure_file_hint)


@dataclass(slots=True)
class FitResult:
    """Result payload returned by the package regression helpers."""

    problem: FitProblem
    fitted_values: dict[str, float] = field(default_factory=dict)
    rendered_values: dict[str, str | float] = field(default_factory=dict)
    metrics_by_term: dict[str, float] = field(default_factory=dict)
    cost: float = float("nan")
    residual_norm: float = float("nan")
    success: bool = False
    status: int = 0
    message: str = ""
    nfev: int = 0
    backend: str = "least_squares_native"
    optimizer_backend: str = "least_squares_native"
    derivative_backend: str = "not_available"
    objective_initial: float = float("nan")
    objective_final: float = float("nan")
    n_residual_evaluations: int = 0
    n_jacobian_evaluations: int = 0
    gradient_norm: float = float("nan")
    step_norm: float = float("nan")
    python_objective_used: bool = False
    parameter_map: dict[str, float] = field(default_factory=dict)
    row_diagnostics: list[dict[str, Any]] = field(default_factory=list)
    active_bounds: list[str] = field(default_factory=list)
    jacobian_available: bool = False
    jacobian_backend: str = "not_available"
    jacobian_fallback_used: bool = False
    jacobian_fallback_reason: str = ""
    backend_unavailable_reason: str = ""
    hessian_available: bool = False
    hessian_backend: str = "not_implemented"
    hessian_fallback_used: bool = False
    hessian_fallback_reason: str = HESSIAN_NOT_IMPLEMENTED_REASON
    provenance_report: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.fitted_values = {str(k): float(v) for k, v in self.fitted_values.items()}
        rendered: dict[str, str | float] = {}
        for key, value in self.rendered_values.items():
            rendered[str(key)] = (
                float(value) if isinstance(value, (int, float, np.integer, np.floating)) else str(value)
            )
        self.rendered_values = rendered
        self.metrics_by_term = {str(k): float(v) for k, v in self.metrics_by_term.items()}
        self.cost = float(self.cost)
        self.residual_norm = float(self.residual_norm)
        self.success = bool(self.success)
        self.status = int(self.status)
        self.message = str(self.message)
        self.nfev = int(self.nfev)
        self.backend = str(self.backend)
        self.optimizer_backend = str(self.optimizer_backend)
        self.derivative_backend = str(self.derivative_backend)
        self.objective_initial = float(self.objective_initial)
        self.objective_final = float(self.objective_final)
        self.n_residual_evaluations = int(self.n_residual_evaluations)
        self.n_jacobian_evaluations = int(self.n_jacobian_evaluations)
        self.gradient_norm = float(self.gradient_norm)
        self.step_norm = float(self.step_norm)
        self.python_objective_used = bool(self.python_objective_used)
        self.parameter_map = {str(k): float(v) for k, v in self.parameter_map.items()}
        self.row_diagnostics = [dict(row) for row in self.row_diagnostics]
        self.active_bounds = [str(item) for item in self.active_bounds]
        self.jacobian_available = bool(self.jacobian_available)
        self.jacobian_backend = str(self.jacobian_backend)
        self.jacobian_fallback_used = bool(self.jacobian_fallback_used)
        self.jacobian_fallback_reason = str(self.jacobian_fallback_reason)
        self.backend_unavailable_reason = str(self.backend_unavailable_reason)
        self.hessian_available = bool(self.hessian_available)
        self.hessian_backend = str(self.hessian_backend)
        self.hessian_fallback_used = bool(self.hessian_fallback_used)
        self.hessian_fallback_reason = str(self.hessian_fallback_reason)
        self.provenance_report = _copy_mapping(self.provenance_report)


@dataclass(slots=True)
class ReactiveElectrolyteRegressionResult:
    """Fixed-shape residual payload for coupled reactive-electrolyte objectives."""

    residuals: np.ndarray
    residual_names: tuple[str, ...] = field(default_factory=tuple)
    record_results: tuple[dict[str, Any], ...] = field(default_factory=tuple, repr=False)
    diagnostics: dict[str, Any] = field(default_factory=dict)
    success_count: int = 0
    failure_count: int = 0
    backend: str = "python_orchestrated_native_reactive_electrolyte"

    def __post_init__(self) -> None:
        self.residuals = np.asarray(self.residuals, dtype=float).reshape(-1)
        self.residual_names = tuple(str(name) for name in self.residual_names)
        self.record_results = tuple(dict(row) for row in self.record_results)
        self.diagnostics = _copy_mapping(self.diagnostics)
        self.success_count = int(self.success_count)
        self.failure_count = int(self.failure_count)
        self.backend = str(self.backend)
        if self.residual_names and len(self.residual_names) != int(self.residuals.size):
            raise InputError("ReactiveElectrolyteRegressionResult residual_names must match residual length.")

    @property
    def residual_norm(self) -> float:
        """Return the Euclidean norm of the fixed residual vector."""
        return float(np.linalg.norm(self.residuals))

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-like result payload."""
        return {
            "residuals": [float(value) for value in self.residuals],
            "residual_names": list(self.residual_names),
            "residual_norm": self.residual_norm,
            "success_count": int(self.success_count),
            "failure_count": int(self.failure_count),
            "backend": self.backend,
            "record_results": _json_like_regression(self.record_results),
            "diagnostics": _json_like_regression(self.diagnostics),
        }


def _fit_derivative_metadata(result: Mapping[str, Any]) -> dict[str, Any]:
    """Extract compact Jacobian/Hessian metadata from native regression payloads."""

    return {
        "optimizer_backend": str(result.get("optimizer_backend", result.get("backend", "not_available"))),
        "derivative_backend": str(result.get("derivative_backend", result.get("jacobian_backend", "not_available"))),
        "objective_initial": float(result.get("objective_initial", result.get("initial_cost", float("nan")))),
        "objective_final": float(result.get("objective_final", result.get("cost", float("nan")))),
        "n_residual_evaluations": int(result.get("n_residual_evaluations", result.get("objective_evaluations", 0))),
        "n_jacobian_evaluations": int(result.get("n_jacobian_evaluations", result.get("gradient_evaluations", 0))),
        "gradient_norm": float(result.get("gradient_norm", float("nan"))),
        "step_norm": float(result.get("step_norm", float("nan"))),
        "python_objective_used": bool(result.get("python_objective_used", False)),
        "jacobian_available": bool(result.get("jacobian_available", False)),
        "jacobian_backend": str(result.get("jacobian_backend", "not_available")),
        "jacobian_fallback_used": bool(result.get("jacobian_fallback_used", False)),
        "jacobian_fallback_reason": str(result.get("jacobian_fallback_reason", "")),
        "backend_unavailable_reason": str(result.get("backend_unavailable_reason", "")),
        "hessian_available": bool(result.get("hessian_available", False)),
        "hessian_backend": str(result.get("hessian_backend", "not_implemented")),
        "hessian_fallback_used": bool(result.get("hessian_fallback_used", False)),
        "hessian_fallback_reason": str(result.get("hessian_fallback_reason", HESSIAN_NOT_IMPLEMENTED_REASON)),
    }


def load_regression_records(records: Any) -> list[dict[str, Any]]:
    """Load flat regression records from CSV, tabular objects, or mappings."""

    if isinstance(records, (str, Path)):
        path = Path(records).expanduser()
        if not path.exists():
            raise FileNotFoundError(f"Regression record file not found: {path}")
        if path.suffix.lower() != ".csv":
            raise InputError("Only CSV file inputs are supported for file-driven regression records.")
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            return [dict(row) for row in reader]

    if hasattr(records, "to_dict"):
        try:
            payload = records.to_dict("records")
        except TypeError:
            payload = records.to_dict()
        if isinstance(payload, list):
            return [dict(row) for row in payload]

    if isinstance(records, Mapping):
        return [dict(records)]

    try:
        items = list(records)
    except TypeError as exc:
        raise InputError("records must be a CSV path, a tabular object, or an iterable of mappings.") from exc

    if not items:
        raise InputError("At least one regression record is required.")
    normalized: list[dict[str, Any]] = []
    for item in items:
        if not isinstance(item, Mapping):
            raise InputError("Regression record iterables must contain mapping-like items.")
        normalized.append(dict(item))
    return normalized


def _normalize_records(records: Any) -> list[dict[str, Any]]:
    normalized = load_regression_records(records)
    if not normalized:
        raise InputError("At least one regression record is required.")
    return normalized


def _assoc_is_enabled(assoc_scheme: str | None) -> bool:
    token = str(assoc_scheme or "").strip().lower()
    return token not in {"", "none", "null", "0"}


def _normalize_fit_targets(mode: str, fit_targets: Iterable[str] | None, assoc_scheme: str = "") -> tuple[str, ...]:
    if fit_targets is None:
        if mode == PURE_NEUTRAL_MODE:
            return tuple(DEFAULT_TARGETS[PURE_NEUTRAL_MODE]["nonassociating"])
        return tuple(DEFAULT_TARGETS[mode])
    names = tuple(str(name) for name in fit_targets)
    if not names:
        raise InputError("fit_targets must include at least one target name.")
    return names


def _normalize_temperature_model(token: str | None) -> str:
    value = str(token or "constant").strip().lower()
    if value not in {"constant", "linear"}:
        raise InputError("temperature_model must be 'constant' or 'linear'.")
    return value


def _optimization_parameter_names(mode: str, fit_targets: Sequence[str], temperature_model: str) -> tuple[str, ...]:
    if mode != BINARY_PAIR_MODE or temperature_model == "constant":
        return tuple(str(name) for name in fit_targets)
    names: list[str] = []
    for name in fit_targets:
        names.append(f"{name}_slope")
        names.append(f"{name}_intercept")
    return tuple(names)


def _coerce_bounds(bounds: FitBounds | Mapping[str, tuple[float | None, float | None]] | None) -> FitBounds:
    if bounds is None:
        return FitBounds()
    if isinstance(bounds, FitBounds):
        return FitBounds(lower=bounds.lower, upper=bounds.upper)
    lower: dict[str, float] = {}
    upper: dict[str, float] = {}
    for name, pair in bounds.items():
        lo, hi = pair
        if lo is not None:
            lower[str(name)] = float(lo)
        if hi is not None:
            upper[str(name)] = float(hi)
    return FitBounds(lower=lower, upper=upper)


def _value_from_record(record: Mapping[str, Any], *keys: str, required: bool = False) -> Any:
    for key in keys:
        if key in record and record[key] not in (None, ""):
            return record[key]
    if required:
        raise InputError(f"Regression record is missing one of the required keys: {', '.join(keys)}.")
    return None


def _float_from_record(record: Mapping[str, Any], *keys: str, required: bool = False) -> float | None:
    raw = _value_from_record(record, *keys, required=required)
    if raw is None:
        return None
    try:
        value = float(raw)
    except (TypeError, ValueError) as exc:
        raise InputError(f"Expected a numeric value for one of {keys}, got {raw!r}.") from exc
    if not math.isfinite(value):
        raise InputError(f"Expected a finite numeric value for one of {keys}, got {raw!r}.")
    return value


def _prefixed_species_values(record: Mapping[str, Any], prefix: str) -> dict[str, float]:
    values: dict[str, float] = {}
    for key, raw in record.items():
        if not str(key).startswith(prefix):
            continue
        name = _normalize_component(str(key)[len(prefix) :])
        if not name:
            continue
        try:
            value = float(raw)
        except (TypeError, ValueError) as exc:
            raise InputError(f"Expected a numeric value for column '{key}', got {raw!r}.") from exc
        if not math.isfinite(value):
            raise InputError(f"Expected a finite numeric value for column '{key}', got {raw!r}.")
        values[name] = value
    return values


def _infer_species_union(records: Sequence[Mapping[str, Any]], prefixes: Sequence[str]) -> tuple[str, ...]:
    ordered: list[str] = []
    seen: set[str] = set()
    for record in records:
        for prefix in prefixes:
            for species in _prefixed_species_values(record, prefix):
                if species not in seen:
                    seen.add(species)
                    ordered.append(species)
    if not ordered:
        raise InputError("Could not infer any species columns from the regression records.")
    return tuple(ordered)


def _composition_from_record(record: Mapping[str, Any], prefix: str, species: Sequence[str]) -> np.ndarray:
    raw = _prefixed_species_values(record, prefix)
    values = {str(name): raw.get(str(name)) for name in species}
    present = {name: value for name, value in values.items() if value is not None}
    if not present:
        raise InputError(f"Regression record is missing composition columns with prefix '{prefix}'.")
    if len(present) == len(species) - 1:
        missing = [name for name, value in values.items() if value is None]
        values[missing[0]] = 1.0 - sum(float(value) for value in present.values())
    elif len(present) != len(species):
        missing = [name for name, value in values.items() if value is None]
        raise InputError(
            f"Regression record is missing composition columns for: {', '.join(missing)} with prefix '{prefix}'."
        )

    array = np.asarray([float(values[name]) for name in species], dtype=float)
    if np.any(~np.isfinite(array)):
        raise InputError(f"Composition values for prefix '{prefix}' must be finite.")
    if np.any(array < -1.0e-12):
        raise InputError(f"Composition values for prefix '{prefix}' must be non-negative.")
    total = float(np.sum(array))
    if total <= 0.0:
        raise InputError(f"Composition values for prefix '{prefix}' must sum to a positive number.")
    array = np.clip(array / total, 0.0, None)
    return array


def _family_scale(term: FitTerm) -> float:
    if term.residual_count <= 0:
        return 1.0
    return math.sqrt(float(term.weight) / float(term.residual_count))


def _safe_log_fraction(value: float) -> float:
    if value <= 0.0:
        raise InputError("Fugacity-balance records require strictly positive composition values.")
    return math.log(value)


def _best_pair_label(mapping: Mapping[str, float], record: Mapping[str, Any]) -> float:
    explicit = _value_from_record(record, "pair_label", "mean_ionic_label", "salt_label", required=False)
    if explicit is not None:
        key = str(explicit)
        if key not in mapping:
            raise InputError(f"Requested mean-ionic label '{key}' is not present in the calculated state.")
        return float(mapping[key])
    if len(mapping) != 1:
        labels = ", ".join(sorted(mapping))
        raise InputError(
            f"Regression record must specify pair_label when multiple mean-ionic labels are available: {labels}."
        )
    return float(next(iter(mapping.values())))


def _seed_value(name: str, initial_guess: Mapping[str, float], current: Mapping[str, Any]) -> float:
    if name in initial_guess:
        value = float(initial_guess[name])
    elif name in current and current[name] not in (None, ""):
        value = float(current[name])
    else:
        raise InputError(f"An initial guess is required for regression target '{name}'.")
    if not math.isfinite(value):
        raise InputError(f"Initial guess for '{name}' must be finite.")
    return value


def _binary_seed_value(
    target: str,
    temperature_model: str,
    initial_guess: Mapping[str, float],
    current: Mapping[str, float],
) -> dict[str, float]:
    if temperature_model == "constant":
        return {target: float(initial_guess.get(target, current.get(target, 0.0)))}
    slope_name = f"{target}_slope"
    intercept_name = f"{target}_intercept"
    slope = float(initial_guess.get(slope_name, 0.0))
    intercept = float(initial_guess.get(intercept_name, current.get(target, 0.0)))
    return {slope_name: slope, intercept_name: intercept}


def _pure_seed_payload(
    component: str,
    T_ref: float,
    assoc_scheme: str,
    dataset: str | Path | None,
    pure_set: str | None,
) -> tuple[dict[str, Any], str | None]:
    source_key: str | None = None
    payload: dict[str, Any] = {}
    if dataset is not None:
        dataset_obj = _load_dataset(dataset)
        for field in PURE_REQUIRED_FIELDS:
            value, source = _resolve_component_field_with_source(
                dataset_obj,
                component,
                field,
                T_ref,
                pure_set_key=pure_set,
            )
            payload[field] = value
            if source_key is None and source is not None:
                source_key = source
    for field in PURE_REQUIRED_FIELDS:
        if field in payload and payload[field] not in (None, ""):
            continue
        default = _deterministic_default(component, field, T_ref)
        if default is not None and default is not _MISSING:
            payload[field] = default
    payload.setdefault("assoc_scheme", assoc_scheme)
    payload.setdefault("z", 0.0)
    payload.setdefault("dielc", 8.0)
    payload.setdefault("d_born", 0.0)
    payload.setdefault("f_solv", 1.0)
    return payload, source_key


def _build_single_component_params(component: str, values: Mapping[str, Any], assoc_scheme: str) -> dict[str, Any]:
    params: dict[str, Any] = {}
    for field in PURE_REQUIRED_FIELDS:
        value = values.get(field)
        if value is None:
            raise InputError(f"Missing required pure-component parameter '{field}' for {component}.")
        params[field] = np.asarray([float(value)], dtype=float)
    assoc_value = assoc_scheme or values.get("assoc_scheme", "")
    params["assoc_scheme"] = [None if not _assoc_is_enabled(str(assoc_value or "")) else str(assoc_value)]
    params["k_ij"] = np.zeros((1, 1), dtype=float)
    params["l_ij"] = np.zeros((1, 1), dtype=float)
    params["k_hb"] = np.zeros((1, 1), dtype=float)
    if abs(float(params["z"][0])) <= 1.0e-12:
        params["z"] = np.asarray([], dtype=float)
    return params


def _normalize_vector_map(names: Sequence[str], values: Sequence[float]) -> dict[str, float]:
    return {str(name): float(value) for name, value in zip(names, values)}


def _render_binary_values(
    vector_map: Mapping[str, float], fit_targets: Sequence[str], temperature_model: str
) -> dict[str, str | float]:
    rendered: dict[str, str | float] = {}
    if temperature_model == "constant":
        for target in fit_targets:
            rendered[str(target)] = float(vector_map[str(target)])
        return rendered
    for target in fit_targets:
        slope = float(vector_map[f"{target}_slope"])
        intercept = float(vector_map[f"{target}_intercept"])
        sign = "+" if intercept >= 0.0 else "-"
        rendered[str(target)] = f"{slope:.12g}*T {sign} {abs(intercept):.12g}"
    return rendered


def _term_summary(records: Sequence[dict[str, Any]], family: str, weight: float, residual_count: int) -> FitTerm:
    return FitTerm(term_type=family, records=tuple(records), weight=weight, residual_count=residual_count)


def _require_record_value(records: Sequence[dict[str, Any]], family: str, key: str) -> None:
    for record in records:
        if _value_from_record(record, key, required=False) in (None, ""):
            raise InputError(f"{family} regression records require a '{key}' value for every record in that family.")


def _build_pure_neutral_terms(records: Sequence[dict[str, Any]]) -> tuple[FitTerm, ...]:
    density_records = [
        record
        for record in records
        if _value_from_record(record, *PURE_DENSITY_KEYS_MOLAR, *PURE_DENSITY_KEYS_MASS, required=False) is not None
    ]
    saturation_records = [record for record in records if _value_from_record(record, "P", required=False) is not None]
    if not density_records:
        raise InputError(
            "pure_neutral regression requires at least one density record with a molar-density 'rho' value or a "
            "mass-density value such as 'rho_kg_m3' or 'rho_sat_liq_kg_m3'."
        )
    if not saturation_records:
        raise InputError(
            "pure_neutral regression requires at least one saturation record with an experimental 'P' value."
        )
    _require_record_value(density_records, PURE_NEUTRAL_MODE, "P")
    return (
        _term_summary(density_records, TERM_DENSITY, 1.0, len(density_records)),
        _term_summary(saturation_records, TERM_PURE_VLE, 1.0, len(saturation_records)),
    )


def _pure_neutral_density_molar(record: Mapping[str, Any], molecular_weight: float) -> float:
    rho_molar = _float_from_record(record, *PURE_DENSITY_KEYS_MOLAR, required=False)
    if rho_molar is not None:
        return rho_molar
    rho_mass = _float_from_record(record, *PURE_DENSITY_KEYS_MASS, required=True)
    if molecular_weight <= 0.0:
        raise InputError("pure_neutral regression requires a positive MW value when converting mass density data.")
    return rho_mass / molecular_weight


def _build_pure_ion_terms(records: Sequence[dict[str, Any]]) -> tuple[FitTerm, ...]:
    density_records = [
        record
        for record in records
        if _value_from_record(record, *PURE_DENSITY_KEYS_MOLAR, *PURE_DENSITY_KEYS_MASS, required=False) is not None
    ]
    osmotic_records = [
        record
        for record in records
        if _value_from_record(record, "osmotic_coefficient", "osmotic", required=False) is not None
    ]
    miac_records = [
        record
        for record in records
        if _value_from_record(
            record,
            "mean_ionic_activity",
            "mean_ionic_activity_coefficient",
            "miac",
            required=False,
        )
        is not None
    ]
    if not osmotic_records and not miac_records:
        raise InputError("pure_ion regression requires osmotic and/or mean-ionic activity records.")
    _require_record_value(density_records, PURE_ION_MODE, "P")
    _require_record_value(osmotic_records, PURE_ION_MODE, "P")
    _require_record_value(miac_records, PURE_ION_MODE, "P")
    terms: list[FitTerm] = []
    if density_records:
        terms.append(_term_summary(density_records, TERM_DENSITY, 1.0, len(density_records)))
    if osmotic_records:
        terms.append(_term_summary(osmotic_records, TERM_OSMOTIC, 1.0, len(osmotic_records)))
    if miac_records:
        terms.append(_term_summary(miac_records, TERM_MIAC, 1.0, len(miac_records)))
    return tuple(terms)


def _build_binary_terms(records: Sequence[dict[str, Any]]) -> tuple[FitTerm, ...]:
    vle_records = [record for record in records if _prefixed_species_values(record, "y_")]
    lle_records = [
        record
        for record in records
        if _prefixed_species_values(record, "x_alpha_") and _prefixed_species_values(record, "x_beta_")
    ]
    if lle_records:
        raise InputError("binary_pair V1 supports only VLE x/y records; LLE records are not supported yet.")
    if not vle_records:
        raise InputError("binary_pair regression requires VLE records with y_* columns.")
    _require_record_value(vle_records, BINARY_PAIR_MODE, "P")
    return (_term_summary(vle_records, TERM_BINARY_VLE, 1.0, 2 * len(vle_records)),)


def _native_pure_neutral_solver_payload(
    normalized_records: Sequence[dict[str, Any]],
    normalized_component: str,
    assoc_scheme: str,
    dataset: str | Path | None,
    pure_set: str | None,
    normalized_fit_targets: Sequence[str],
    fixed_parameters: Mapping[str, Any] | None,
    initial_guess: Mapping[str, float] | None,
    bounds: FitBounds | Mapping[str, tuple[float | None, float | None]] | None,
) -> tuple[dict[str, Any], dict[str, float], FitBounds, tuple[FitTerm, ...], tuple[str, ...], str | None]:
    if _assoc_is_enabled(assoc_scheme):
        raise InputError("The native pure_neutral workflow currently supports only nonassociating neutral components.")

    bounds_obj = _coerce_bounds(bounds)
    T_ref = float(np.mean([_float_from_record(record, "T", required=True) for record in normalized_records]))
    seed_payload, source_key = _pure_seed_payload(normalized_component, T_ref, assoc_scheme, dataset, pure_set)
    fixed_payload = seed_payload.copy()
    fixed_payload.update(_copy_mapping(fixed_parameters))
    fixed_payload["assoc_scheme"] = str(assoc_scheme or fixed_payload.get("assoc_scheme", ""))
    if "MW" not in fixed_payload or fixed_payload["MW"] in (None, ""):
        raise InputError(
            "pure_neutral regression requires a fixed MW value, either from the dataset or fixed_parameters."
        )
    if "z" not in fixed_payload or fixed_payload["z"] in (None, ""):
        fixed_payload["z"] = 0.0

    initial = _copy_mapping(initial_guess)
    initial_map = {target: _seed_value(target, initial, fixed_payload) for target in normalized_fit_targets}
    optimization_names = _optimization_parameter_names(PURE_NEUTRAL_MODE, normalized_fit_targets, "constant")
    for name in ("m", "s", "e"):
        if name in initial_map:
            fixed_payload[name] = float(initial_map[name])
    terms = _build_pure_neutral_terms(normalized_records)
    pure_file_hint = (
        f"{source_key}.csv" if source_key is not None else _infer_pure_template_name([normalized_component])
    )
    return (
        _ensure_native_vector_payload(fixed_payload),
        initial_map,
        bounds_obj,
        terms,
        optimization_names,
        pure_file_hint,
    )


def _family_metrics(raw_residuals_by_term: Mapping[str, Sequence[float]]) -> dict[str, float]:
    metrics: dict[str, float] = {}
    for family, values in raw_residuals_by_term.items():
        arr = np.asarray(values, dtype=float)
        metrics[family] = float(np.sqrt(np.mean(arr**2))) if arr.size else 0.0
    return metrics


def _json_like_regression(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return _json_like_regression(value.tolist())
    if isinstance(value, Mapping):
        return {str(key): _json_like_regression(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_like_regression(item) for item in value]
    if isinstance(value, (bool, str)) or value is None:
        return value
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        out = float(value)
        return out if math.isfinite(out) else str(out)
    return str(value)


def _source_dataset_label(dataset: str | Path | None) -> str | None:
    if dataset is None:
        return None
    return str(dataset)


def _ensure_native_vector_payload(payload: Mapping[str, Any]) -> dict[str, Any]:
    normalized = dict(payload)
    for key in ("m", "s", "e", "e_assoc", "vol_a", "z", "dielc", "d_born", "f_solv", "MW"):
        if key not in normalized:
            continue
        value = normalized[key]
        if isinstance(value, np.ndarray):
            normalized[key] = np.asarray(value, dtype=float).reshape(-1)
        elif isinstance(value, (list, tuple)):
            normalized[key] = np.asarray(value, dtype=float).reshape(-1)
        else:
            normalized[key] = np.asarray([value], dtype=float)
    return normalized


def _normalize_species_list(species: Iterable[str] | None) -> tuple[str, ...] | None:
    if species is None:
        return None
    normalized = tuple(_normalize_component(str(name)) for name in species)
    if not normalized:
        raise InputError("species must include at least one component.")
    return normalized


def _normalize_pair(pair: Sequence[str]) -> tuple[str, str]:
    names = tuple(_normalize_component(str(name)) for name in pair)
    if len(names) != 2:
        raise InputError("binary_pair regression requires exactly two pair components.")
    if names[0] == names[1]:
        raise InputError("binary_pair regression requires two distinct pair components.")
    return names


def _record_has_molality(record: Mapping[str, Any]) -> bool:
    return _value_from_record(record, "molality", "m_salt", "salt_molality", required=False) is not None


def _ion_species_from_records(records: Sequence[Mapping[str, Any]], species: Iterable[str] | None) -> tuple[str, ...]:
    normalized = _normalize_species_list(species)
    if normalized is not None:
        return normalized
    if any(_prefixed_species_values(record, "x_") for record in records):
        return _infer_species_union(records, ("x_",))
    if any(_record_has_molality(record) for record in records):
        raise InputError("molality-driven pure_ion records require explicit species and solvent arguments.")
    raise InputError(
        "pure_ion records require composition columns with prefix 'x_' or molality with species and solvent."
    )


def _binary_species_from_records(
    records: Sequence[Mapping[str, Any]],
    pair: tuple[str, str],
    species: Iterable[str] | None,
) -> tuple[str, ...]:
    normalized = _normalize_species_list(species)
    inferred = _infer_species_union(records, ("x_", "y_"))
    if normalized is None:
        normalized = inferred
    missing_pair = [name for name in pair if name not in normalized]
    if missing_pair:
        raise InputError(f"binary_pair species must include fitted pair components: {', '.join(missing_pair)}.")
    missing_records = [name for name in normalized if name not in inferred]
    if missing_records:
        raise InputError(f"binary_pair VLE records are missing x_/y_ columns for: {', '.join(missing_records)}.")
    return normalized


def _ion_composition_from_record(record: Mapping[str, Any], species: Sequence[str], solvent: str | None) -> np.ndarray:
    if _prefixed_species_values(record, "x_"):
        return _composition_from_record(record, "x_", species)
    molality = _float_from_record(record, "molality", "m_salt", "salt_molality", required=False)
    if molality is None:
        raise InputError("pure_ion records require composition columns with prefix 'x_' or a molality value.")
    if solvent is None:
        raise InputError("molality-driven pure_ion records require a solvent argument.")
    try:
        return np.asarray(molality_to_molefraction(molality, species=species, solvent=solvent), dtype=float)
    except ValueError as exc:
        raise InputError(str(exc)) from exc


def _params_for_native_record(
    dataset: str | Path,
    species: Sequence[str],
    x: np.ndarray,
    T: float,
    user_options: Mapping[str, Any] | None,
) -> dict[str, Any]:
    return get_prop_dict(dataset, species, x, T, user_options=_copy_mapping(user_options))


def _native_target_payload(
    optimization_names: Sequence[str],
    species: Sequence[str],
    *,
    component: str | None = None,
    pair: tuple[str, str] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    kinds: list[int] = []
    indices: list[int] = []
    indices_2: list[int] = []
    for name in optimization_names:
        if name not in NATIVE_TARGET_KINDS:
            raise InputError(f"Native regression does not support optimization parameter '{name}'.")
        kinds.append(NATIVE_TARGET_KINDS[name])
        if name in {"k_ij", "l_ij", "k_hb_ij"}:
            if pair is None:
                raise InputError(f"Native {name} regression requires a fitted pair.")
            indices.append(species.index(pair[0]))
            indices_2.append(species.index(pair[1]))
        else:
            if component is None:
                raise InputError(f"Native pure-parameter regression requires a component for '{name}'.")
            indices.append(species.index(component))
            indices_2.append(-1)
    return (
        np.asarray(kinds, dtype=int),
        np.asarray(indices, dtype=int),
        np.asarray(indices_2, dtype=int),
    )


def _native_density_record(
    term_name: str,
    record: Mapping[str, Any],
    x: np.ndarray,
    scale: float,
    *,
    phase: str | None = None,
) -> dict[str, Any] | None:
    rho_molar = _float_from_record(record, *PURE_DENSITY_KEYS_MOLAR, required=False)
    rho_mass = _float_from_record(record, *PURE_DENSITY_KEYS_MASS, required=False)
    if rho_molar is None and rho_mass is None:
        return None
    target = rho_molar if rho_molar is not None else rho_mass
    return {
        "term_name": term_name,
        "term": NATIVE_TERM_KINDS[term_name],
        "T": _float_from_record(record, "T", required=True),
        "P": _float_from_record(record, "P", required=True),
        "phase": phase_to_int(phase or _value_from_record(record, "phase", required=False) or "liq"),
        "x": np.asarray(x, dtype=float).tolist(),
        "target": float(target),
        "density_kind": 0 if rho_molar is not None else 1,
        "scale": scale,
    }


def _native_miac_pair_indices(record: Mapping[str, Any], species: Sequence[str]) -> tuple[int, int]:
    explicit = _value_from_record(record, "pair_label", "mean_ionic_label", "salt_label", required=False)
    if explicit is None:
        return -1, -1
    label = str(explicit)
    for i, left in enumerate(species):
        for j, right in enumerate(species):
            if f"{left}{right}" == label:
                return i, j
    raise InputError(f"Requested mean-ionic label '{label}' is not present in the fitted species list.")


def _run_native_generic_least_squares(
    fixed_payloads: Sequence[Mapping[str, Any]],
    native_records: Sequence[Mapping[str, Any]],
    optimization_names: Sequence[str],
    species: Sequence[str],
    theta0: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    *,
    component: str | None = None,
    pair: tuple[str, str] | None = None,
    multistart: int = 0,
    max_nfev: int = 200,
) -> dict[str, Any]:
    target_kinds, target_indices, target_indices_2 = _native_target_payload(
        optimization_names,
        species,
        component=component,
        pair=pair,
    )
    return _fit_generic_native_least_squares(
        [dict(payload) for payload in fixed_payloads],
        [dict(record) for record in native_records],
        target_kinds,
        target_indices,
        target_indices_2,
        theta0,
        lower,
        upper,
        multistart=int(multistart),
        max_nfev=int(max_nfev),
    )


def _run_native_generic_ceres(
    fixed_payloads: Sequence[Mapping[str, Any]],
    native_records: Sequence[Mapping[str, Any]],
    optimization_names: Sequence[str],
    species: Sequence[str],
    theta0: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    *,
    component: str | None = None,
    pair: tuple[str, str] | None = None,
    multistart: int = 0,
    max_nfev: int = 200,
) -> dict[str, Any]:
    target_kinds, target_indices, target_indices_2 = _native_target_payload(
        optimization_names,
        species,
        component=component,
        pair=pair,
    )
    return _fit_generic_native_ceres(
        [dict(payload) for payload in fixed_payloads],
        [dict(record) for record in native_records],
        target_kinds,
        target_indices,
        target_indices_2,
        theta0,
        lower,
        upper,
        multistart=int(multistart),
        max_nfev=int(max_nfev),
    )


def _run_native_generic_derivative_debug(
    fixed_payloads: Sequence[Mapping[str, Any]],
    native_records: Sequence[Mapping[str, Any]],
    optimization_names: Sequence[str],
    species: Sequence[str],
    x: Sequence[float],
    *,
    component: str | None = None,
    pair: tuple[str, str] | None = None,
) -> dict[str, Any]:
    target_kinds, target_indices, target_indices_2 = _native_target_payload(
        optimization_names,
        species,
        component=component,
        pair=pair,
    )
    return _evaluate_generic_native_debug(
        [dict(payload) for payload in fixed_payloads],
        [dict(record) for record in native_records],
        target_kinds,
        target_indices,
        target_indices_2,
        np.asarray(x, dtype=float),
    )


def evaluate_generic_regression_derivatives(
    fixed_payloads: Sequence[Mapping[str, Any]],
    native_records: Sequence[Mapping[str, Any]],
    optimization_names: Sequence[str],
    species: Sequence[str],
    x: Sequence[float],
    *,
    component: str | None = None,
    pair: Sequence[str] | None = None,
    jacobian_backend: str = "auto",
) -> dict[str, Any]:
    """Evaluate native generic-regression residual derivative payloads where available."""

    backend = str(jacobian_backend).strip().lower()
    if backend not in {"auto", "autodiff", "analytic"}:
        raise InputError("jacobian_backend must be 'auto', 'autodiff', or 'analytic'.")
    raise InputError(
        "backend_unavailable: generic regression sensitivities are not implemented for this residual path."
    )

    normalized_species = tuple(_normalize_component(str(name)) for name in species)
    normalized_pair = None if pair is None else _normalize_pair(pair)
    normalized_component = None if component is None else _normalize_component(component)
    normalized_names = tuple(str(name) for name in optimization_names)
    theta = np.asarray(x, dtype=float)
    if theta.size != len(normalized_names):
        raise InputError(f"Expected {len(normalized_names)} values for x, got {theta.size}.")
    return _run_native_generic_derivative_debug(
        fixed_payloads,
        native_records,
        normalized_names,
        normalized_species,
        theta,
        component=normalized_component,
        pair=normalized_pair,
    )


def _terms_support_ion_activity(terms: Sequence[FitTerm]) -> bool:
    return any(term.term_type in {TERM_OSMOTIC, TERM_MIAC} for term in terms)


def _terms_support_dielectric(terms: Sequence[FitTerm]) -> bool:
    return any(term.term_type == TERM_RELATIVE_PERMITTIVITY for term in terms)


def _source_for_pure_target(target: str, terms: Sequence[FitTerm]) -> str:
    if target == "d_born":
        if _terms_support_dielectric(terms):
            return "relative_permittivity"
        if _terms_support_ion_activity(terms):
            return "ion_activity"
        return "mixed_reactive_vle"
    if _terms_support_ion_activity(terms):
        return "ion_activity"
    return "direct_data"


def _pure_parameter_declarations(
    component: str,
    targets: Sequence[str],
    terms: Sequence[FitTerm],
    *,
    provenance: Sequence[FitParameter] | Mapping[str, str] | None,
    allow_unsupported_parameters: bool,
) -> list[FitParameter]:
    declared: dict[str, FitParameter] = {}
    if isinstance(provenance, Mapping):
        for target, source in provenance.items():
            declared[str(target)] = FitParameter(
                component,
                str(target),
                source=str(source),
                allow_without_direct_data=allow_unsupported_parameters,
            )
    elif provenance is not None:
        for item in provenance:
            if not isinstance(item, FitParameter):
                raise InputError("pure-parameter provenance entries must be FitParameter instances.")
            declared[item.parameter] = item
    out: list[FitParameter] = []
    for target in targets:
        out.append(
            declared.get(
                target,
                FitParameter(
                    component,
                    target,
                    source=_source_for_pure_target(target, terms),
                    allow_without_direct_data=allow_unsupported_parameters,
                ),
            )
        )
    return out


def _binary_interaction_declarations(
    pair: tuple[str, str],
    targets: Sequence[str],
    *,
    provenance: Sequence[BinaryInteraction] | Mapping[str, str] | None,
    allow_unsupported_parameters: bool,
) -> list[BinaryInteraction]:
    declared: dict[str, BinaryInteraction] = {}
    if isinstance(provenance, Mapping):
        for target, source in provenance.items():
            declared[str(target)] = BinaryInteraction(
                pair,
                parameter=str(target),
                source=str(source),
                allow_without_direct_data=allow_unsupported_parameters,
            )
    elif provenance is not None:
        for item in provenance:
            if not isinstance(item, BinaryInteraction):
                raise InputError("binary-interaction provenance entries must be BinaryInteraction instances.")
            declared[item.parameter] = item
    out: list[BinaryInteraction] = []
    for target in targets:
        out.append(
            declared.get(
                target,
                BinaryInteraction(
                    pair,
                    parameter=target,
                    source="direct_binary_vle",
                    allow_without_direct_data=allow_unsupported_parameters,
                ),
            )
        )
    return out


def _fit_pure_ion_internal(
    records: Any,
    component: str,
    *,
    dataset: str | Path,
    solvent: str | None = None,
    species: Iterable[str] | None = None,
    fit_targets: Iterable[str] | None = None,
    initial_guess: Mapping[str, float] | None = None,
    bounds: FitBounds | Mapping[str, tuple[float | None, float | None]] | None = None,
    user_options: Mapping[str, Any] | None = None,
    provenance: Sequence[FitParameter] | Mapping[str, str] | None = None,
    allow_unsupported_parameters: bool = False,
    multistart: int = 0,
) -> FitResult:
    normalized_component = _normalize_component(component)
    normalized_solvent = None if solvent is None else _normalize_component(solvent)
    normalized_records = _normalize_records(records)
    normalized_species = _ion_species_from_records(normalized_records, species)
    if normalized_component not in normalized_species:
        raise InputError(f"pure_ion species must include fitted component '{normalized_component}'.")
    if normalized_solvent is not None and normalized_solvent not in normalized_species:
        raise InputError(f"solvent '{normalized_solvent}' is not present in species.")

    normalized_fit_targets = _normalize_fit_targets(PURE_ION_MODE, fit_targets)
    unsupported = [target for target in normalized_fit_targets if target not in {"s", "e", "d_born"}]
    if unsupported:
        raise InputError("pure_ion V1 supports only the targets 's', 'e', and 'd_born'.")
    terms = _build_pure_ion_terms(normalized_records)
    provenance_report = validate_regression_provenance(
        _pure_parameter_declarations(
            normalized_component,
            normalized_fit_targets,
            terms,
            provenance=provenance,
            allow_unsupported_parameters=allow_unsupported_parameters,
        ),
        terms=terms,
        species=normalized_species,
        strict=True,
    )

    T_ref = float(np.mean([_float_from_record(record, "T", required=True) for record in normalized_records]))
    seed_payload, source_key = _pure_seed_payload(normalized_component, T_ref, "", dataset, None)
    initial = _copy_mapping(initial_guess)
    initial_map = {target: _seed_value(target, initial, seed_payload) for target in normalized_fit_targets}
    optimization_names = _optimization_parameter_names(PURE_ION_MODE, normalized_fit_targets, "constant")
    bounds_obj = _coerce_bounds(bounds)
    lower, upper = bounds_obj.arrays_for(optimization_names)
    theta0 = np.asarray([initial_map[name] for name in optimization_names], dtype=float)

    native_records: list[dict[str, Any]] = []
    fixed_payloads: list[dict[str, Any]] = []
    solvent_index = -1 if normalized_solvent is None else normalized_species.index(normalized_solvent)
    for term in terms:
        scale = _family_scale(term)
        for record in term.records:
            T = _float_from_record(record, "T", required=True)
            P = _float_from_record(record, "P", required=True)
            assert T is not None and P is not None
            x = _ion_composition_from_record(record, normalized_species, normalized_solvent)
            if term.term_type == TERM_DENSITY:
                native_record = _native_density_record(term.term_type, record, x, scale)
                if native_record is None:
                    continue
            elif term.term_type == TERM_OSMOTIC:
                native_record = {
                    "term_name": term.term_type,
                    "term": NATIVE_TERM_KINDS[term.term_type],
                    "T": T,
                    "P": P,
                    "phase": phase_to_int(_value_from_record(record, "phase", required=False) or "liq"),
                    "x": x.tolist(),
                    "target": _float_from_record(record, "osmotic_coefficient", "osmotic", required=True),
                    "scale": scale,
                }
            elif term.term_type == TERM_MIAC:
                basis = (
                    str(_value_from_record(record, "activity_basis", "miac_basis", required=False) or "molality")
                    .strip()
                    .lower()
                )
                miac_i, miac_j = _native_miac_pair_indices(record, normalized_species)
                native_record = {
                    "term_name": term.term_type,
                    "term": NATIVE_TERM_KINDS[term.term_type],
                    "T": T,
                    "P": P,
                    "phase": phase_to_int(_value_from_record(record, "phase", required=False) or "liq"),
                    "x": x.tolist(),
                    "target": _float_from_record(
                        record,
                        "mean_ionic_activity",
                        "mean_ionic_activity_coefficient",
                        "miac",
                        required=True,
                    ),
                    "target_index": miac_i,
                    "target_index_2": miac_j,
                    "activity_basis": 1 if basis in {"molality", "m"} else 0,
                    "solvent_index": solvent_index,
                    "scale": scale,
                }
            else:
                continue
            fixed_payloads.append(_params_for_native_record(dataset, normalized_species, x, T, user_options))
            native_records.append(native_record)

    result = _run_native_generic_least_squares(
        fixed_payloads,
        native_records,
        optimization_names,
        normalized_species,
        theta0,
        lower,
        upper,
        component=normalized_component,
        multistart=int(multistart),
    )
    vector_map = _normalize_vector_map(optimization_names, result["x"])
    rendered = {name: float(vector_map[name]) for name in normalized_fit_targets}
    problem = FitProblem(
        mode=PURE_ION_MODE,
        records=tuple(normalized_records),
        component=normalized_component,
        solvent=normalized_solvent,
        dataset=_source_dataset_label(dataset),
        fit_targets=normalized_fit_targets,
        optimization_parameters=optimization_names,
        fixed_parameters=seed_payload,
        initial_guess=initial_map,
        terms=terms,
        pure_file_hint=(
            f"{source_key}.csv" if source_key is not None else _infer_pure_template_name([normalized_component])
        ),
    )
    return FitResult(
        problem=problem,
        fitted_values=vector_map,
        rendered_values=rendered,
        metrics_by_term=result["metrics_by_term"],
        cost=float(result["cost"]),
        residual_norm=float(result["residual_norm"]),
        success=bool(result["success"]),
        status=int(result["status"]),
        message=str(result["message"]),
        nfev=int(result["nfev"]),
        backend=str(result["backend"]),
        **_fit_derivative_metadata(result),
        provenance_report=provenance_report,
    )


def _fit_binary_pair_internal(
    records: Any,
    pair: Sequence[str],
    *,
    dataset: str | Path,
    species: Iterable[str] | None = None,
    fit_targets: Iterable[str] | None = None,
    temperature_model: str = "constant",
    initial_guess: Mapping[str, float] | None = None,
    bounds: FitBounds | Mapping[str, tuple[float | None, float | None]] | None = None,
    user_options: Mapping[str, Any] | None = None,
    provenance: Sequence[BinaryInteraction] | Mapping[str, str] | None = None,
    allow_unsupported_parameters: bool = False,
    multistart: int = 0,
    optimizer_backend: str = "least_squares_native",
) -> FitResult:
    optimizer_backend = _optimizer_backend_from_options(
        {"optimizer_backend": optimizer_backend}, "least_squares_native"
    )
    normalized_pair = _normalize_pair(pair)
    normalized_records = _normalize_records(records)
    normalized_species = _binary_species_from_records(normalized_records, normalized_pair, species)
    normalized_temperature_model = _normalize_temperature_model(temperature_model)
    if normalized_temperature_model != "constant":
        raise InputError("binary_pair V1 supports only temperature_model='constant'.")
    normalized_fit_targets = _normalize_fit_targets(BINARY_PAIR_MODE, fit_targets)
    unsupported_targets = [target for target in normalized_fit_targets if target not in {"k_ij", "l_ij", "k_hb_ij"}]
    if unsupported_targets:
        raise InputError("binary_pair regression supports only the targets 'k_ij', 'l_ij', and 'k_hb_ij'.")
    terms = _build_binary_terms(normalized_records)

    T_ref = float(np.mean([_float_from_record(record, "T", required=True) for record in normalized_records]))
    sample_x = _composition_from_record(normalized_records[0], "x_", normalized_species)
    sample_payload = _params_for_native_record(dataset, normalized_species, sample_x, T_ref, user_options)
    charges = np.asarray(sample_payload.get("z", np.zeros(len(normalized_species))), dtype=float).flatten()
    if charges.size == 0:
        charges = np.zeros(len(normalized_species), dtype=float)
    provenance_report = validate_regression_provenance(
        _binary_interaction_declarations(
            normalized_pair,
            normalized_fit_targets,
            provenance=provenance,
            allow_unsupported_parameters=allow_unsupported_parameters,
        ),
        terms=terms,
        species=normalized_species,
        charges=charges,
        strict=True,
    )
    dataset_obj = _load_dataset(dataset)
    current = {
        target: _matrix_value(
            dataset_obj, "k_hb" if target == "k_hb_ij" else target, normalized_pair[0], normalized_pair[1], T_ref
        )
        for target in normalized_fit_targets
    }
    initial_map = {}
    for target in normalized_fit_targets:
        initial_map.update(
            _binary_seed_value(target, normalized_temperature_model, _copy_mapping(initial_guess), current)
        )
    optimization_names = _optimization_parameter_names(
        BINARY_PAIR_MODE, normalized_fit_targets, normalized_temperature_model
    )
    bounds_obj = _coerce_bounds(bounds)
    lower, upper = bounds_obj.arrays_for(optimization_names)
    theta0 = np.asarray([initial_map[name] for name in optimization_names], dtype=float)
    pair_indices = tuple(normalized_species.index(name) for name in normalized_pair)
    problem = FitProblem(
        mode=BINARY_PAIR_MODE,
        records=tuple(normalized_records),
        pair=normalized_pair,
        dataset=_source_dataset_label(dataset),
        fit_targets=normalized_fit_targets,
        optimization_parameters=optimization_names,
        initial_guess=initial_map,
        temperature_model=normalized_temperature_model,
        terms=terms,
    )

    native_records: list[dict[str, Any]] = []
    fixed_payloads: list[dict[str, Any]] = []
    for term in terms:
        scale = _family_scale(term)
        for record in term.records:
            T = _float_from_record(record, "T", required=True)
            P = _float_from_record(record, "P", required=True)
            assert T is not None and P is not None
            x_liq = _composition_from_record(record, "x_", normalized_species)
            y_vap = _composition_from_record(record, "y_", normalized_species)
            for index in pair_indices:
                _safe_log_fraction(float(x_liq[index]))
                _safe_log_fraction(float(y_vap[index]))
            native_records.append(
                {
                    "term_name": term.term_type,
                    "term": NATIVE_TERM_KINDS[term.term_type],
                    "T": T,
                    "P": P,
                    "x": x_liq.tolist(),
                    "y": y_vap.tolist(),
                    "target_index": pair_indices[0],
                    "target_index_2": pair_indices[1],
                    "scale": scale,
                }
            )
            fixed_payloads.append(_params_for_native_record(dataset, normalized_species, x_liq, T, user_options))

    if optimizer_backend == "ceres":
        if normalized_fit_targets != ("k_ij",):
            raise InputError("backend_unavailable: binary_pair Ceres currently supports only the target 'k_ij'.")
        result = _run_native_generic_ceres(
            fixed_payloads,
            native_records,
            optimization_names,
            normalized_species,
            theta0,
            lower,
            upper,
            pair=normalized_pair,
            multistart=int(multistart),
        )
    else:
        result = _run_native_generic_least_squares(
            fixed_payloads,
            native_records,
            optimization_names,
            normalized_species,
            theta0,
            lower,
            upper,
            pair=normalized_pair,
            multistart=int(multistart),
        )
    vector_map = _normalize_vector_map(optimization_names, result["x"])
    return FitResult(
        problem=problem,
        fitted_values=vector_map,
        rendered_values=_render_binary_values(vector_map, normalized_fit_targets, normalized_temperature_model),
        metrics_by_term=result["metrics_by_term"],
        cost=float(result["cost"]),
        residual_norm=float(result["residual_norm"]),
        success=bool(result["success"]),
        status=int(result["status"]),
        message=str(result["message"]),
        nfev=int(result["nfev"]),
        backend=str(result["backend"]),
        **_fit_derivative_metadata(result),
        provenance_report=provenance_report,
    )


def _native_pure_neutral_runner_args(
    normalized_records: Sequence[dict[str, Any]],
    normalized_component: str,
    assoc_scheme: str,
    dataset: str | Path | None,
    pure_set: str | None,
    normalized_fit_targets: Sequence[str],
    fixed_parameters: Mapping[str, Any] | None,
    initial_guess: Mapping[str, float] | None,
    bounds: FitBounds | Mapping[str, tuple[float | None, float | None]] | None,
):
    fixed_payload, initial_map, bounds_obj, terms, optimization_names, pure_file_hint = (
        _native_pure_neutral_solver_payload(
            normalized_records,
            normalized_component,
            assoc_scheme,
            dataset,
            pure_set,
            normalized_fit_targets,
            fixed_parameters,
            initial_guess,
            bounds,
        )
    )
    lower, upper = bounds_obj.arrays_for(optimization_names)
    density_term = next(term for term in terms if term.term_type == TERM_DENSITY)
    pure_vle_term = next(term for term in terms if term.term_type == TERM_PURE_VLE)
    mw_value = float(np.asarray(fixed_payload["MW"], dtype=float)[0])
    return {
        "fixed_payload": fixed_payload,
        "initial_map": initial_map,
        "terms": terms,
        "optimization_names": optimization_names,
        "pure_file_hint": pure_file_hint,
        "lower": lower,
        "upper": upper,
        "density_T": np.asarray(
            [_float_from_record(record, "T", required=True) for record in density_term.records], dtype=float
        ),
        "density_P": np.asarray(
            [_float_from_record(record, "P", required=True) for record in density_term.records], dtype=float
        ),
        "density_rho_exp": np.asarray(
            [_pure_neutral_density_molar(record, mw_value) for record in density_term.records],
            dtype=float,
        ),
        "density_phase": np.asarray(
            [
                phase_to_int(_value_from_record(record, "phase", required=False) or "liq")
                for record in density_term.records
            ],
            dtype=int,
        ),
        "density_scale": float(_family_scale(density_term)),
        "vle_T": np.asarray(
            [_float_from_record(record, "T", required=True) for record in pure_vle_term.records], dtype=float
        ),
        "vle_P": np.asarray(
            [_float_from_record(record, "P", required=True) for record in pure_vle_term.records], dtype=float
        ),
        "pure_vle_scale": float(_family_scale(pure_vle_term)),
    }


def _fit_pure_neutral_native_ceres(
    fixed_payload,
    density_T,
    density_P,
    density_rho_exp,
    density_phase,
    density_scale,
    vle_T,
    vle_P,
    pure_vle_scale,
    x0,
    lower,
    upper,
    multistart=0,
) -> dict[str, Any]:
    params = check_association(dict(fixed_payload))
    cppargs = create_struct(params)
    result = _core._fit_pure_neutral_native_ceres(
        cppargs,
        np.asarray(density_T, dtype=float),
        np.asarray(density_P, dtype=float),
        np.asarray(density_rho_exp, dtype=float),
        np.asarray(density_phase, dtype=int),
        float(density_scale),
        np.asarray(vle_T, dtype=float),
        np.asarray(vle_P, dtype=float),
        float(pure_vle_scale),
        np.asarray(x0, dtype=float),
        np.asarray(lower, dtype=float),
        np.asarray(upper, dtype=float),
        int(multistart),
    )
    return {
        "x": vector_to_array(result["x"]),
        "cost": float(result["cost"]),
        "residual_norm": float(result["residual_norm"]),
        "density_metric": float(result["density_metric"]),
        "pure_vle_metric": float(result["pure_vle_metric"]),
        "initial_cost": float(result["initial_cost"]),
        "initial_density_metric": float(result["initial_density_metric"]),
        "initial_pure_vle_metric": float(result["initial_pure_vle_metric"]),
        "success": bool(result["success"]),
        "status": int(result["status"]),
        "nfev": int(result["nfev"]),
        "iterations": int(result["iterations"]),
        "starts_tried": int(result["starts_tried"]),
        "objective_evaluations": int(result["objective_evaluations"]),
        "gradient_evaluations": int(result["gradient_evaluations"]),
        "residual_evaluations": int(result["residual_evaluations"]),
        "density_solves": int(result["density_solves"]),
        "fused_state_evaluations": int(result["fused_state_evaluations"]),
        "callback_wall_time_s": float(result["callback_wall_time_s"]),
        "solve_wall_time_s": float(result["solve_wall_time_s"]),
        "message": str(result["message"]),
        "backend": str(result["backend"]),
        "optimizer_backend": str(result.get("optimizer_backend", result["backend"])),
        "derivative_backend": str(result.get("derivative_backend", result.get("jacobian_backend", "not_available"))),
        "objective_initial": float(result.get("objective_initial", result["initial_cost"])),
        "objective_final": float(result.get("objective_final", result["cost"])),
        "n_residual_evaluations": int(result.get("n_residual_evaluations", result["objective_evaluations"])),
        "n_jacobian_evaluations": int(result.get("n_jacobian_evaluations", result["gradient_evaluations"])),
        "gradient_norm": float(result.get("gradient_norm", float("nan"))),
        "step_norm": float(result.get("step_norm", float("nan"))),
        "python_objective_used": bool(result.get("python_objective_used", False)),
        "jacobian_available": bool(result.get("jacobian_available", True)),
        "jacobian_backend": str(result.get("jacobian_backend", result.get("derivative_backend", "not_available"))),
        "jacobian_fallback_used": bool(result.get("jacobian_fallback_used", False)),
        "jacobian_fallback_reason": str(result.get("jacobian_fallback_reason", "")),
        "hessian_available": bool(result.get("hessian_available", False)),
        "hessian_backend": str(result.get("hessian_backend", "not_implemented")),
        "hessian_fallback_used": bool(result.get("hessian_fallback_used", False)),
        "hessian_fallback_reason": str(result.get("hessian_fallback_reason", HESSIAN_NOT_IMPLEMENTED_REASON)),
    }


def _fit_pure_neutral_internal(
    records: Any,
    component: str,
    *,
    assoc_scheme: str = "",
    dataset: str | Path | None = None,
    pure_set: str | None = None,
    fit_targets: Iterable[str] | None = None,
    fixed_parameters: Mapping[str, Any] | None = None,
    initial_guess: Mapping[str, float] | None = None,
    bounds: FitBounds | Mapping[str, tuple[float | None, float | None]] | None = None,
    multistart: int = 0,
    optimizer_backend: str = "least_squares_native",
) -> FitResult:
    fit_result, _ = _fit_pure_neutral_internal_with_native(
        records,
        component,
        assoc_scheme=assoc_scheme,
        dataset=dataset,
        pure_set=pure_set,
        fit_targets=fit_targets,
        fixed_parameters=fixed_parameters,
        initial_guess=initial_guess,
        bounds=bounds,
        multistart=multistart,
        optimizer_backend=optimizer_backend,
    )
    return fit_result


def _fit_pure_neutral_internal_with_native(
    records: Any,
    component: str,
    *,
    assoc_scheme: str = "",
    dataset: str | Path | None = None,
    pure_set: str | None = None,
    fit_targets: Iterable[str] | None = None,
    fixed_parameters: Mapping[str, Any] | None = None,
    initial_guess: Mapping[str, float] | None = None,
    bounds: FitBounds | Mapping[str, tuple[float | None, float | None]] | None = None,
    multistart: int = 0,
    optimizer_backend: str = "least_squares_native",
) -> tuple[FitResult, dict[str, Any]]:
    normalized_component = _normalize_component(component)
    normalized_records = _normalize_records(records)
    normalized_fit_targets = _normalize_fit_targets(PURE_NEUTRAL_MODE, fit_targets, assoc_scheme=assoc_scheme)
    for target in normalized_fit_targets:
        if target not in {"m", "s", "e"}:
            raise InputError("Phase-1 pure_neutral regression supports only the targets 'm', 's', and 'e'.")

    payload = _native_pure_neutral_runner_args(
        normalized_records,
        normalized_component,
        assoc_scheme,
        dataset,
        pure_set,
        normalized_fit_targets,
        fixed_parameters,
        initial_guess,
        bounds,
    )
    theta0 = np.asarray([payload["initial_map"][name] for name in payload["optimization_names"]], dtype=float)
    if optimizer_backend == "least_squares_native":
        native_result = _fit_pure_neutral_native_least_squares(
            payload["fixed_payload"],
            payload["density_T"],
            payload["density_P"],
            payload["density_rho_exp"],
            payload["density_phase"],
            payload["density_scale"],
            payload["vle_T"],
            payload["vle_P"],
            payload["pure_vle_scale"],
            theta0,
            payload["lower"],
            payload["upper"],
            multistart=int(multistart),
        )
    elif optimizer_backend == "ceres":
        native_result = _fit_pure_neutral_native_ceres(
            payload["fixed_payload"],
            payload["density_T"],
            payload["density_P"],
            payload["density_rho_exp"],
            payload["density_phase"],
            payload["density_scale"],
            payload["vle_T"],
            payload["vle_P"],
            payload["pure_vle_scale"],
            theta0,
            payload["lower"],
            payload["upper"],
            multistart=int(multistart),
        )
    else:
        raise InputError(f"Unsupported pure-neutral optimizer backend: {optimizer_backend}")

    vector_map = _normalize_vector_map(payload["optimization_names"], native_result["x"])
    problem = FitProblem(
        mode=PURE_NEUTRAL_MODE,
        records=tuple(normalized_records),
        component=normalized_component,
        dataset=_source_dataset_label(dataset),
        fit_targets=normalized_fit_targets,
        optimization_parameters=payload["optimization_names"],
        fixed_parameters=payload["fixed_payload"],
        initial_guess=payload["initial_map"],
        assoc_scheme=str(payload["fixed_payload"]["assoc_scheme"]),
        terms=payload["terms"],
        pure_file_hint=payload["pure_file_hint"],
    )
    rendered = {name: float(vector_map[name]) for name in normalized_fit_targets}
    row_diagnostics = [
        {"row_family": TERM_DENSITY, "metric": float(native_result["density_metric"])},
        {"row_family": TERM_PURE_VLE, "metric": float(native_result["pure_vle_metric"])},
    ]
    fit_result = FitResult(
        problem=problem,
        fitted_values=vector_map,
        rendered_values=rendered,
        metrics_by_term={
            TERM_DENSITY: float(native_result["density_metric"]),
            TERM_PURE_VLE: float(native_result["pure_vle_metric"]),
        },
        cost=float(native_result["cost"]),
        residual_norm=float(native_result["residual_norm"]),
        success=bool(native_result["success"]),
        status=int(native_result["status"]),
        message=str(native_result["message"]),
        nfev=int(native_result["nfev"]),
        backend=str(native_result["backend"]),
        parameter_map=vector_map,
        row_diagnostics=row_diagnostics,
        active_bounds=[
            name
            for name in payload["optimization_names"]
            if math.isclose(vector_map[name], float(payload["lower"][payload["optimization_names"].index(name)]))
            or math.isclose(vector_map[name], float(payload["upper"][payload["optimization_names"].index(name)]))
        ],
        **_fit_derivative_metadata(native_result),
    )
    return fit_result, native_result


def fit_pure_neutral(
    records: Any,
    component: str,
    *,
    assoc_scheme: str = "",
    dataset: str | Path | None = None,
    pure_set: str | None = None,
    fit_targets: Iterable[str] | None = None,
    fixed_parameters: Mapping[str, Any] | None = None,
    initial_guess: Mapping[str, float] | None = None,
    bounds: FitBounds | Mapping[str, tuple[float | None, float | None]] | None = None,
    user_options: Mapping[str, Any] | None = None,
    multistart: int = 0,
    optimizer_backend: str = "least_squares_native",
) -> FitResult:
    """Fit neutral pure-component m, s, and e parameters."""
    _ = user_options
    return _fit_pure_neutral_internal(
        records,
        component,
        assoc_scheme=assoc_scheme,
        dataset=dataset,
        pure_set=pure_set,
        fit_targets=fit_targets,
        fixed_parameters=fixed_parameters,
        initial_guess=initial_guess,
        bounds=bounds,
        multistart=multistart,
        optimizer_backend=_optimizer_backend_from_options(
            {"optimizer_backend": optimizer_backend}, "least_squares_native"
        ),
    )


def _fit_pure_neutral_least_squares_internal(
    records: Any,
    component: str,
    *,
    assoc_scheme: str = "",
    dataset: str | Path | None = None,
    pure_set: str | None = None,
    fit_targets: Iterable[str] | None = None,
    fixed_parameters: Mapping[str, Any] | None = None,
    initial_guess: Mapping[str, float] | None = None,
    bounds: FitBounds | Mapping[str, tuple[float | None, float | None]] | None = None,
    multistart: int = 0,
) -> FitResult:
    """Internal native least-squares comparison hook for pure-neutral regression."""
    return _fit_pure_neutral_internal(
        records,
        component,
        assoc_scheme=assoc_scheme,
        dataset=dataset,
        pure_set=pure_set,
        fit_targets=fit_targets,
        fixed_parameters=fixed_parameters,
        initial_guess=initial_guess,
        bounds=bounds,
        multistart=multistart,
        optimizer_backend="least_squares_native",
    )


def _associating_pure_payload(
    component: str,
    T_ref: float,
    assoc_scheme: str,
    fixed_parameters: Mapping[str, Any] | None,
    initial_guess: Mapping[str, float] | None,
) -> dict[str, Any]:
    payload, _ = _pure_seed_payload(component, T_ref, assoc_scheme, None, None)
    payload.update(_copy_mapping(fixed_parameters))
    initial = _copy_mapping(initial_guess)
    for field in PURE_REQUIRED_FIELDS:
        if field in payload and payload[field] not in (None, ""):
            continue
        if field in initial:
            payload[field] = initial[field]
    payload["assoc_scheme"] = assoc_scheme
    payload.setdefault("z", 0.0)
    payload.setdefault("dielc", 8.0)
    payload.setdefault("d_born", 0.0)
    payload.setdefault("f_solv", 1.0)
    missing = [field for field in PURE_REQUIRED_FIELDS if field not in payload or payload[field] in (None, "")]
    if missing:
        raise InputError(f"Associating pure-neutral regression is missing fixed values for: {', '.join(missing)}.")
    return payload


def _fit_pure_neutral_associating_native(
    records: Any,
    component: str,
    *,
    assoc_scheme: str,
    fit_targets: Iterable[str] | None = None,
    fixed_parameters: Mapping[str, Any] | None = None,
    initial_guess: Mapping[str, float] | None = None,
    bounds: FitBounds | Mapping[str, tuple[float | None, float | None]] | None = None,
    multistart: int = 0,
    max_nfev: int = 1,
) -> FitResult:
    """Internal native associating pure-neutral regression path for analysis checks."""

    if not _assoc_is_enabled(assoc_scheme):
        raise InputError("Associating pure-neutral regression requires an association scheme.")
    normalized_component = _normalize_component(component)
    normalized_records = _normalize_records(records)
    normalized_fit_targets = (
        tuple(DEFAULT_TARGETS[PURE_NEUTRAL_MODE]["associating"])
        if fit_targets is None
        else _normalize_fit_targets(PURE_NEUTRAL_MODE, fit_targets, assoc_scheme=assoc_scheme)
    )
    unsupported = [target for target in normalized_fit_targets if target not in {"m", "s", "e", "e_assoc", "vol_a"}]
    if unsupported:
        raise InputError("Associating pure-neutral regression supports only m, s, e, e_assoc, and vol_a.")
    terms = _build_pure_neutral_terms(normalized_records)
    T_ref = float(np.mean([_float_from_record(record, "T", required=True) for record in normalized_records]))
    seed_payload = _associating_pure_payload(
        normalized_component,
        T_ref,
        str(assoc_scheme),
        fixed_parameters,
        initial_guess,
    )
    initial = _copy_mapping(initial_guess)
    initial_map = {target: _seed_value(target, initial, seed_payload) for target in normalized_fit_targets}
    lower, upper = _coerce_bounds(bounds).arrays_for(normalized_fit_targets)
    theta0 = np.asarray([initial_map[name] for name in normalized_fit_targets], dtype=float)

    params = _build_single_component_params(normalized_component, seed_payload, str(assoc_scheme))
    fixed_payloads: list[dict[str, Any]] = []
    native_records: list[dict[str, Any]] = []
    x_single = np.asarray([1.0], dtype=float)
    for term in terms:
        scale = _family_scale(term)
        for record in term.records:
            if term.term_type == TERM_DENSITY:
                native_record = _native_density_record(term.term_type, record, x_single, scale)
                if native_record is None:
                    continue
            elif term.term_type == TERM_PURE_VLE:
                native_record = {
                    "term_name": term.term_type,
                    "term": NATIVE_TERM_KINDS[term.term_type],
                    "T": _float_from_record(record, "T", required=True),
                    "P": _float_from_record(record, "P", required=True),
                    "x": x_single.tolist(),
                    "scale": scale,
                }
            else:
                continue
            fixed_payloads.append(params)
            native_records.append(native_record)

    result = _run_native_generic_least_squares(
        fixed_payloads,
        native_records,
        normalized_fit_targets,
        (normalized_component,),
        theta0,
        lower,
        upper,
        component=normalized_component,
        multistart=int(multistart),
        max_nfev=int(max_nfev),
    )
    vector_map = _normalize_vector_map(normalized_fit_targets, result["x"])
    metrics = dict(result["metrics_by_term"])
    metrics["initial_residual_norm"] = float(result["initial_residual_norm"])
    problem = FitProblem(
        mode=PURE_NEUTRAL_MODE,
        records=tuple(normalized_records),
        component=normalized_component,
        fit_targets=normalized_fit_targets,
        optimization_parameters=normalized_fit_targets,
        fixed_parameters=seed_payload,
        initial_guess=initial_map,
        assoc_scheme=str(assoc_scheme),
        terms=tuple(terms),
        pure_file_hint=_infer_pure_template_name([normalized_component]),
    )
    return FitResult(
        problem=problem,
        fitted_values=vector_map,
        rendered_values={name: float(vector_map[name]) for name in normalized_fit_targets},
        metrics_by_term=metrics,
        cost=float(result["cost"]),
        residual_norm=float(result["residual_norm"]),
        success=bool(result["success"]),
        status=int(result["status"]),
        message=str(result["message"]),
        nfev=int(result["nfev"]),
        backend=str(result["backend"]),
        **_fit_derivative_metadata(result),
    )


def _debug_native_pure_neutral_objective(
    records: Any,
    component: str,
    *,
    assoc_scheme: str = "",
    dataset: str | Path | None = None,
    pure_set: str | None = None,
    fit_targets: Iterable[str] | None = None,
    fixed_parameters: Mapping[str, Any] | None = None,
    initial_guess: Mapping[str, float] | None = None,
    bounds: FitBounds | Mapping[str, tuple[float | None, float | None]] | None = None,
    x: Mapping[str, float] | Sequence[float] | None = None,
) -> dict[str, Any]:
    """Internal debug hook for validating the native pure-neutral objective and gradient."""

    normalized_component = _normalize_component(component)
    normalized_records = _normalize_records(records)
    normalized_fit_targets = _normalize_fit_targets(PURE_NEUTRAL_MODE, fit_targets, assoc_scheme=assoc_scheme)
    for target in normalized_fit_targets:
        if target not in {"m", "s", "e"}:
            raise InputError("Phase-1 pure_neutral regression supports only the targets 'm', 's', and 'e'.")
    payload = _native_pure_neutral_runner_args(
        normalized_records,
        normalized_component,
        assoc_scheme,
        dataset,
        pure_set,
        normalized_fit_targets,
        fixed_parameters,
        initial_guess,
        bounds,
    )
    if x is None:
        theta = np.asarray([payload["initial_map"][name] for name in payload["optimization_names"]], dtype=float)
    elif isinstance(x, Mapping):
        theta = np.asarray([float(x[name]) for name in payload["optimization_names"]], dtype=float)
    else:
        theta = np.asarray(tuple(x), dtype=float)
    if theta.size != len(payload["optimization_names"]):
        raise InputError(
            f"Expected {len(payload['optimization_names'])} optimization values for native debug evaluation, got {theta.size}."
        )
    return _fit_pure_neutral_native_debug(
        payload["fixed_payload"],
        payload["density_T"],
        payload["density_P"],
        payload["density_rho_exp"],
        payload["density_phase"],
        payload["density_scale"],
        payload["vle_T"],
        payload["vle_P"],
        payload["pure_vle_scale"],
        theta,
    )


def evaluate_pure_neutral_derivatives(
    records: Any,
    component: str,
    *,
    assoc_scheme: str = "",
    dataset: str | Path | None = None,
    pure_set: str | None = None,
    fit_targets: Iterable[str] | None = None,
    fixed_parameters: Mapping[str, Any] | None = None,
    initial_guess: Mapping[str, float] | None = None,
    bounds: FitBounds | Mapping[str, tuple[float | None, float | None]] | None = None,
    x: Mapping[str, float] | Sequence[float] | None = None,
) -> dict[str, Any]:
    """Evaluate pure-neutral residuals, gradient, and Jacobian diagnostics at a parameter point."""

    return _debug_native_pure_neutral_objective(
        records,
        component,
        assoc_scheme=assoc_scheme,
        dataset=dataset,
        pure_set=pure_set,
        fit_targets=fit_targets,
        fixed_parameters=fixed_parameters,
        initial_guess=initial_guess,
        bounds=bounds,
        x=x,
    )


_USER_TARGET_ALIASES = {
    "sigma": "s",
    "epsilon": "e",
    "epsilon_k": "e",
    "epsilon_over_k": "e",
    "association_energy": "e_assoc",
    "association_volume": "vol_a",
    "born_diameter": "d_born",
    "solvation_factor": "f_solv",
}


def _reject_finite_difference_options(options: Any) -> None:
    tokens = ("finite_difference", "finite-difference", "finite difference")

    def visit(value: Any) -> bool:
        if isinstance(value, str):
            return any(token in value.lower() for token in tokens)
        if isinstance(value, Mapping):
            return any(visit(key) or visit(item) for key, item in value.items())
        if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
            return any(visit(item) for item in value)
        return False

    if visit(options):
        raise InputError("finite-difference derivative options are not part of the public regression API.")


def _optimizer_backend_from_options(options: Mapping[str, Any] | None, default: str) -> str:
    backend = str((options or {}).get("optimizer_backend", default)).strip().lower()
    aliases = {
        "native": "least_squares_native",
        "least_squares": "least_squares_native",
        "least_squares_native": "least_squares_native",
        "ceres": "ceres",
    }
    if backend not in aliases:
        raise InputError(f"Unsupported optimizer_backend: {backend}")
    return aliases[backend]


def _normalize_user_targets(parameters_to_fit: Iterable[str] | None, default: Sequence[str]) -> tuple[str, ...]:
    raw_targets = tuple(default if parameters_to_fit is None else parameters_to_fit)
    return tuple(_USER_TARGET_ALIASES.get(str(target), str(target)) for target in raw_targets)


def _fit_species_tuple(species: str | Sequence[str]) -> tuple[str, ...]:
    if isinstance(species, str):
        return (species,)
    labels = tuple(str(label) for label in species)
    if not labels:
        raise InputError("species must contain at least one component label.")
    return labels


def _annotate_contract_problem(
    result: FitResult,
    *,
    bounds: FitBounds | Mapping[str, tuple[float | None, float | None]] | None,
    weights: Mapping[str, float] | None,
    loss: str,
    solver_options: Mapping[str, Any] | None,
    output_report: bool,
    fixed_parameters: Mapping[str, Any] | None = None,
) -> FitResult:
    result.problem.bounds = _copy_bounds_contract(bounds)
    result.problem.weights = {str(k): float(v) for k, v in (weights or {}).items()}
    result.problem.loss = str(loss)
    result.problem.solver_options = _copy_mapping(solver_options)
    result.problem.output_report = bool(output_report)
    if fixed_parameters is not None:
        result.problem.fixed_parameters = _copy_mapping(fixed_parameters)
    return result


def fit_pure_parameters(
    *,
    species: str | Sequence[str],
    data_rows: Any,
    parameter_set: str | Path | None = None,
    parameters_to_fit: Iterable[str] | None = None,
    fixed_parameters: Mapping[str, Any] | None = None,
    bounds: FitBounds | Mapping[str, tuple[float | None, float | None]] | None = None,
    weights: Mapping[str, float] | None = None,
    loss: str = "linear",
    solver_options: Mapping[str, Any] | None = None,
    output_report: bool = False,
    assoc_scheme: str = "",
    initial_guess: Mapping[str, float] | None = None,
    pure_set: str | None = None,
    user_options: Mapping[str, Any] | None = None,
) -> FitResult:
    """Fit pure-component parameters through the easy regression API contract."""

    _reject_finite_difference_options(solver_options)
    labels = _fit_species_tuple(species)
    if len(labels) != 1:
        raise InputError("fit_pure_parameters expects exactly one species label.")
    targets = _normalize_user_targets(parameters_to_fit, DEFAULT_TARGETS[PURE_NEUTRAL_MODE]["nonassociating"])
    result = fit_pure_neutral(
        data_rows,
        labels[0],
        assoc_scheme=assoc_scheme,
        dataset=parameter_set,
        pure_set=pure_set,
        fit_targets=targets,
        fixed_parameters=fixed_parameters,
        initial_guess=initial_guess,
        bounds=bounds,
        user_options=user_options,
        multistart=int((solver_options or {}).get("multistart", 0)),
    )
    return _annotate_contract_problem(
        result,
        bounds=bounds,
        weights=weights,
        loss=loss,
        solver_options=solver_options,
        output_report=output_report,
        fixed_parameters=fixed_parameters,
    )


def fit_pure_ion(
    records: Any,
    component: str,
    *,
    dataset: str | Path,
    solvent: str | None = None,
    species: Iterable[str] | None = None,
    fit_targets: Iterable[str] | None = None,
    initial_guess: Mapping[str, float] | None = None,
    bounds: FitBounds | Mapping[str, tuple[float | None, float | None]] | None = None,
    user_options: Mapping[str, Any] | None = None,
    provenance: Sequence[FitParameter] | Mapping[str, str] | None = None,
    allow_unsupported_parameters: bool = False,
    multistart: int = 0,
) -> FitResult:
    """Fit ion pure-component parameters against electrolyte records."""
    return _fit_pure_ion_internal(
        records,
        component,
        dataset=dataset,
        solvent=solvent,
        species=species,
        fit_targets=fit_targets,
        initial_guess=initial_guess,
        bounds=bounds,
        user_options=user_options,
        provenance=provenance,
        allow_unsupported_parameters=allow_unsupported_parameters,
        multistart=multistart,
    )


def fit_binary_parameters(
    *,
    species: Sequence[str],
    data_rows: Any,
    parameter_set: str | Path,
    parameters_to_fit: Iterable[str] | None = None,
    fixed_parameters: Mapping[str, Any] | None = None,
    bounds: FitBounds | Mapping[str, tuple[float | None, float | None]] | None = None,
    weights: Mapping[str, float] | None = None,
    loss: str = "linear",
    solver_options: Mapping[str, Any] | None = None,
    output_report: bool = False,
    temperature_model: str = "constant",
    initial_guess: Mapping[str, float] | None = None,
    user_options: Mapping[str, Any] | None = None,
    provenance: Sequence[BinaryInteraction] | Mapping[str, str] | None = None,
    allow_unsupported_parameters: bool = False,
) -> FitResult:
    """Fit binary interaction parameters through the easy regression API contract."""

    _reject_finite_difference_options(solver_options)
    labels = _fit_species_tuple(species)
    if len(labels) != 2:
        raise InputError("fit_binary_parameters expects exactly two species labels.")
    targets = _normalize_user_targets(parameters_to_fit, DEFAULT_TARGETS[BINARY_PAIR_MODE])
    result = fit_binary_pair(
        data_rows,
        labels,
        dataset=parameter_set,
        species=labels,
        fit_targets=targets,
        temperature_model=temperature_model,
        initial_guess=initial_guess,
        bounds=bounds,
        user_options=user_options,
        provenance=provenance,
        allow_unsupported_parameters=allow_unsupported_parameters,
        multistart=int((solver_options or {}).get("multistart", 0)),
        optimizer_backend=_optimizer_backend_from_options(solver_options, "least_squares_native"),
    )
    return _annotate_contract_problem(
        result,
        bounds=bounds,
        weights=weights,
        loss=loss,
        solver_options=solver_options,
        output_report=output_report,
        fixed_parameters=fixed_parameters,
    )


def fit_binary_pair(
    records: Any,
    pair: Sequence[str],
    *,
    dataset: str | Path,
    species: Iterable[str] | None = None,
    fit_targets: Iterable[str] | None = None,
    temperature_model: str = "constant",
    initial_guess: Mapping[str, float] | None = None,
    bounds: FitBounds | Mapping[str, tuple[float | None, float | None]] | None = None,
    user_options: Mapping[str, Any] | None = None,
    provenance: Sequence[BinaryInteraction] | Mapping[str, str] | None = None,
    allow_unsupported_parameters: bool = False,
    multistart: int = 0,
    optimizer_backend: str = "least_squares_native",
) -> FitResult:
    """Fit V1 binary interaction parameters against VLE x/y records."""
    return _fit_binary_pair_internal(
        records,
        pair,
        dataset=dataset,
        species=species,
        fit_targets=fit_targets,
        temperature_model=temperature_model,
        initial_guess=initial_guess,
        bounds=bounds,
        user_options=user_options,
        provenance=provenance,
        allow_unsupported_parameters=allow_unsupported_parameters,
        multistart=multistart,
        optimizer_backend=optimizer_backend,
    )


def _build_liquid_electrolyte_terms(
    records: Sequence[dict[str, Any]], weights: Mapping[str, float] | None
) -> tuple[FitTerm, ...]:
    weight_map = {str(k): float(v) for k, v in (weights or {}).items()}
    families: list[tuple[str, str, tuple[str, ...]]] = [
        (TERM_DENSITY, "density", (*PURE_DENSITY_KEYS_MOLAR, *PURE_DENSITY_KEYS_MASS)),
        (TERM_RELATIVE_PERMITTIVITY, "relative_permittivity", ("epsilon_r_exp", "relative_permittivity")),
        (TERM_OSMOTIC, "osmotic_coefficient", ("osmotic_coefficient", "osmotic")),
        (TERM_MIAC, "mean_ionic_activity", ("mean_ionic_activity", "mean_ionic_activity_coefficient", "miac")),
        (TERM_BINARY_VLE, "binary_vle", ("y_",)),
        (TERM_BINARY_LLE, "binary_lle", ("x_phase2_",)),
    ]
    terms: list[FitTerm] = []
    for term_type, weight_key, keys in families:
        selected: list[dict[str, Any]] = []
        for record in records:
            if any(key.endswith("_") and any(str(name).startswith(key) for name in record) for key in keys):
                selected.append(record)
            elif any(key in record for key in keys if not key.endswith("_")):
                selected.append(record)
        if selected:
            terms.append(
                FitTerm(
                    term_type=term_type,
                    records=tuple(selected),
                    weight=weight_map.get(weight_key, weight_map.get(term_type, 1.0)),
                    residual_count=len(selected),
                )
            )
    if not terms:
        raise InputError("fit_liquid_electrolyte_parameters requires at least one supported data row family.")
    return tuple(terms)


def fit_liquid_electrolyte_parameters(
    *,
    species: Sequence[str],
    data_rows: Any,
    parameter_set: str | Path,
    parameters_to_fit: Iterable[str],
    fixed_parameters: Mapping[str, Any] | None = None,
    bounds: FitBounds | Mapping[str, tuple[float | None, float | None]] | None = None,
    weights: Mapping[str, float] | None = None,
    loss: str = "linear",
    solver_options: Mapping[str, Any] | None = None,
    output_report: bool = False,
    initial_guess: Mapping[str, float] | None = None,
    user_options: Mapping[str, Any] | None = None,
    provenance: Sequence[FitParameter | BinaryInteraction] | Mapping[str, str] | None = None,
) -> FitResult:
    """Declare a liquid-electrolyte regression problem without running unavailable optimizer internals."""

    _reject_finite_difference_options(solver_options)
    labels = _fit_species_tuple(species)
    if len(labels) < 2:
        raise InputError("fit_liquid_electrolyte_parameters expects at least two species labels.")
    records = _normalize_records(data_rows)
    targets = _normalize_user_targets(parameters_to_fit, ())
    problem = FitProblem(
        mode=LIQUID_ELECTROLYTE_MODE,
        records=tuple(records),
        dataset=_source_dataset_label(parameter_set),
        fit_targets=targets,
        optimization_parameters=targets,
        fixed_parameters=fixed_parameters or {},
        initial_guess=initial_guess or {},
        bounds=_copy_bounds_contract(bounds),
        weights=weights or {},
        loss=loss,
        solver_options=solver_options or {},
        output_report=output_report,
        terms=_build_liquid_electrolyte_terms(records, weights),
    )
    provenance_report: dict[str, Any] = {}
    if provenance:
        provenance_report = validate_regression_provenance(provenance, species=labels, strict=False)
    if user_options:
        problem.solver_options = {**problem.solver_options, "user_options": _copy_mapping(user_options)}
    optimizer_backend = _optimizer_backend_from_options(solver_options, "ceres")
    return FitResult(
        problem=problem,
        success=False,
        status=-1,
        message=LIQUID_ELECTROLYTE_BACKEND_UNAVAILABLE_REASON,
        backend="backend_unavailable",
        optimizer_backend=optimizer_backend,
        derivative_backend="backend_unavailable",
        jacobian_backend="backend_unavailable",
        jacobian_available=False,
        jacobian_fallback_used=False,
        python_objective_used=False,
        backend_unavailable_reason=LIQUID_ELECTROLYTE_BACKEND_UNAVAILABLE_REASON,
        row_diagnostics=[
            {
                "row_family": term.term_type,
                "supported": False,
                "backend_unavailable_reason": LIQUID_ELECTROLYTE_BACKEND_UNAVAILABLE_REASON,
            }
            for term in problem.terms
        ],
        provenance_report=provenance_report,
    )


def _find_matching_pure_files(dataset_root: Path, component: str) -> list[Path]:
    pure_dir = dataset_root / "pure"
    if not pure_dir.exists():
        raise FileNotFoundError(f"Dataset folder '{dataset_root}' does not contain a pure/ directory.")
    matches: list[Path] = []
    for path in sorted(pure_dir.glob("*.csv")):
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                if _normalize_component(str(row.get("component", "")).strip()) == component:
                    matches.append(path)
                    break
    return matches


def _choose_pure_file(problem: FitProblem, dataset_root: Path, pure_file: str | Path | None) -> Path:
    if pure_file is not None:
        path = Path(pure_file)
        return path if path.is_absolute() else dataset_root / "pure" / path
    matches = _find_matching_pure_files(dataset_root, str(problem.component))
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise InputError("Multiple pure parameter files contain the target component; pass pure_file explicitly.")
    pure_dir = dataset_root / "pure"
    candidates = sorted(pure_dir.glob("*.csv"))
    if len(candidates) == 1:
        return candidates[0]
    if problem.pure_file_hint is not None:
        hinted = pure_dir / problem.pure_file_hint
        if hinted.exists():
            return hinted
    inferred = pure_dir / _infer_pure_template_name([str(problem.component)])
    if inferred.exists():
        return inferred
    raise InputError("Could not determine which pure parameter CSV should receive the fitted values.")


def _update_csv_row(path: Path, component: str, updates: Mapping[str, str | float], overwrite: bool) -> None:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        rows = [dict(row) for row in reader]

    if "component" not in fieldnames:
        raise InputError(f"CSV file '{path}' does not include a 'component' column.")

    target_row = None
    for row in rows:
        if _normalize_component(str(row.get("component", "")).strip()) == component:
            target_row = row
            break

    if target_row is None:
        target_row = {name: "" for name in fieldnames}
        target_row["component"] = component
        rows.append(target_row)

    for key, value in updates.items():
        if key not in fieldnames:
            raise InputError(f"CSV file '{path}' does not include column '{key}'.")
        existing = str(target_row.get(key, "") or "").strip()
        if existing and not overwrite:
            raise InputError(f"Refusing to overwrite existing value for '{component}' column '{key}' in '{path}'.")
        target_row[key] = value

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _update_matrix_file(
    path: Path,
    pair: tuple[str, str],
    value: str | float,
    overwrite: bool,
) -> None:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        rows = [row for row in reader]

    if not rows:
        raise InputError(f"Matrix file '{path}' is empty.")
    header = rows[0]
    if not header or header[0] != "component":
        raise InputError(f"Matrix file '{path}' must start with a 'component' header column.")
    columns = header[1:]
    normalized_columns = [_normalize_component(name) for name in columns]
    row_lookup = {_normalize_component(row[0]): index for index, row in enumerate(rows[1:], start=1) if row}
    if pair[0] not in row_lookup or pair[1] not in row_lookup:
        raise InputError(f"Matrix file '{path}' does not contain both fitted components.")
    try:
        i_col = normalized_columns.index(pair[0]) + 1
        j_col = normalized_columns.index(pair[1]) + 1
    except ValueError as exc:
        raise InputError(f"Matrix file '{path}' is missing one of the fitted columns.") from exc

    i_row = row_lookup[pair[0]]
    j_row = row_lookup[pair[1]]
    string_value = f"{float(value):.12g}" if isinstance(value, (int, float, np.integer, np.floating)) else str(value)

    for row_index, col_index in ((i_row, j_col), (j_row, i_col)):
        existing = str(rows[row_index][col_index] or "").strip()
        if existing and not overwrite:
            raise InputError(f"Refusing to overwrite existing matrix value in '{path}' without overwrite=True.")
        rows[row_index][col_index] = string_value

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerows(rows)


def _mea_co2_h2o_species_from_records(
    records: Sequence[Mapping[str, Any]], species: Iterable[str] | None
) -> tuple[str, ...]:
    normalized = _normalize_species_list(species)
    inferred = _infer_species_union(records, ("x_",))
    if normalized is None:
        normalized = inferred
    required = ("H2O", "MEA", "CO2", "MEAH+", "MEACOO-", "HCO3-")
    missing_required = [name for name in required if name not in normalized]
    if missing_required:
        raise InputError(f"MEA-CO2-H2O benchmark species are missing: {', '.join(missing_required)}.")
    missing_records = [name for name in normalized if name not in inferred]
    if missing_records:
        raise InputError(f"MEA-CO2-H2O benchmark records are missing x_ columns for: {', '.join(missing_records)}.")
    return normalized


def _build_mea_co2_h2o_terms(records: Sequence[dict[str, Any]]) -> tuple[FitTerm, ...]:
    density_records = [
        record
        for record in records
        if _value_from_record(record, *PURE_DENSITY_KEYS_MOLAR, *PURE_DENSITY_KEYS_MASS, required=False) is not None
    ]
    co2_records = [
        record
        for record in records
        if _value_from_record(record, "lnphi_CO2", "ln_fugacity_coefficient_CO2", required=False) is not None
    ]
    osmotic_records = [
        record
        for record in records
        if _value_from_record(record, "osmotic_coefficient", "osmotic", required=False) is not None
    ]
    if not density_records and not co2_records and not osmotic_records:
        raise InputError("MEA-CO2-H2O benchmark records require density, CO2 fugacity, and/or osmotic targets.")
    _require_record_value(density_records, "MEA-CO2-H2O benchmark", "P")
    _require_record_value(co2_records, "MEA-CO2-H2O benchmark", "P")
    _require_record_value(osmotic_records, "MEA-CO2-H2O benchmark", "P")
    terms: list[FitTerm] = []
    if density_records:
        terms.append(_term_summary(density_records, TERM_MEA_CO2_H2O_DENSITY, 1.0, len(density_records)))
    if co2_records:
        terms.append(_term_summary(co2_records, TERM_MEA_CO2_H2O_CO2_FUGACITY, 1.0, len(co2_records)))
    if osmotic_records:
        terms.append(_term_summary(osmotic_records, TERM_MEA_CO2_H2O_OSMOTIC, 1.0, len(osmotic_records)))
    return tuple(terms)


def _target_bounds(
    targets: Sequence[str], bounds: FitBounds | Mapping[str, tuple[float | None, float | None]] | None
) -> tuple[np.ndarray, np.ndarray]:
    return _coerce_bounds(bounds).arrays_for(targets)


def _benchmark_seed_payloads(
    dataset: str | Path,
    species: Sequence[str],
    T_ref: float,
    components: Sequence[str],
) -> tuple[dict[str, dict[str, Any]], dict[str, str | None]]:
    payloads: dict[str, dict[str, Any]] = {}
    hints: dict[str, str | None] = {}
    for component in components:
        payload, source_key = _pure_seed_payload(component, T_ref, "", dataset, None)
        payloads[component] = payload
        hints[component] = f"{source_key}.csv" if source_key is not None else _infer_pure_template_name(list(species))
    return payloads, hints


def _benchmark_vector_map(targets: Sequence[str], theta: Sequence[float]) -> dict[str, float]:
    return {name: float(value) for name, value in zip(targets, theta)}


def _fit_mea_co2_h2o_component(
    records: Sequence[dict[str, Any]],
    component: str,
    *,
    dataset: str | Path,
    species: Sequence[str],
    fit_targets: Sequence[str],
    seed_payload: Mapping[str, Any],
    pure_file_hint: str | None,
    terms: Sequence[FitTerm],
    initial_guess: Mapping[str, float] | None,
    bounds: FitBounds | Mapping[str, tuple[float | None, float | None]] | None,
    user_options: Mapping[str, Any] | None,
    multistart: int,
    max_nfev: int,
) -> FitResult:
    initial = _copy_mapping(initial_guess)
    initial_map = {target: _seed_value(target, initial, seed_payload) for target in fit_targets}
    lower, upper = _target_bounds(fit_targets, bounds)
    theta0 = np.asarray([initial_map[name] for name in fit_targets], dtype=float)

    native_records: list[dict[str, Any]] = []
    fixed_payloads: list[dict[str, Any]] = []
    co2_index = species.index("CO2")
    for term in terms:
        scale = _family_scale(term)
        for record in term.records:
            T = _float_from_record(record, "T", required=True)
            P = _float_from_record(record, "P", required=True)
            assert T is not None and P is not None
            x = _composition_from_record(record, "x_", species)
            if term.term_type == TERM_MEA_CO2_H2O_DENSITY:
                native_record = _native_density_record(term.term_type, record, x, scale)
                if native_record is None:
                    continue
            elif term.term_type == TERM_MEA_CO2_H2O_CO2_FUGACITY:
                native_record = {
                    "term_name": term.term_type,
                    "term": NATIVE_TERM_KINDS[term.term_type],
                    "T": T,
                    "P": P,
                    "phase": phase_to_int(_value_from_record(record, "phase", required=False) or "liq"),
                    "x": x.tolist(),
                    "target": _float_from_record(record, "lnphi_CO2", "ln_fugacity_coefficient_CO2", required=True),
                    "target_index": co2_index,
                    "scale": scale,
                }
            elif term.term_type == TERM_MEA_CO2_H2O_OSMOTIC:
                native_record = {
                    "term_name": term.term_type,
                    "term": NATIVE_TERM_KINDS[term.term_type],
                    "T": T,
                    "P": P,
                    "phase": phase_to_int(_value_from_record(record, "phase", required=False) or "liq"),
                    "x": x.tolist(),
                    "target": _float_from_record(record, "osmotic_coefficient", "osmotic", required=True),
                    "scale": scale,
                }
            else:
                continue
            fixed_payloads.append(_params_for_native_record(dataset, species, x, T, user_options))
            native_records.append(native_record)

    result = _run_native_generic_least_squares(
        fixed_payloads,
        native_records,
        fit_targets,
        species,
        theta0,
        lower,
        upper,
        component=component,
        multistart=int(multistart),
        max_nfev=int(max_nfev),
    )
    vector_map = _benchmark_vector_map(fit_targets, result["x"])
    mode = PURE_ION_MODE if abs(float(seed_payload.get("z", 0.0))) > 1.0e-12 else PURE_NEUTRAL_MODE
    problem = FitProblem(
        mode=mode,
        records=tuple(records),
        component=component,
        dataset=_source_dataset_label(dataset),
        fit_targets=tuple(fit_targets),
        optimization_parameters=tuple(fit_targets),
        fixed_parameters=seed_payload,
        initial_guess=initial_map,
        assoc_scheme=str(seed_payload.get("assoc_scheme", "")),
        terms=tuple(terms),
        pure_file_hint=pure_file_hint,
    )
    metrics = dict(result["metrics_by_term"])
    metrics["initial_residual_norm"] = float(result["initial_residual_norm"])
    return FitResult(
        problem=problem,
        fitted_values=vector_map,
        rendered_values={name: float(vector_map[name]) for name in fit_targets},
        metrics_by_term=metrics,
        cost=float(result["cost"]),
        residual_norm=float(result["residual_norm"]),
        success=bool(result["success"]),
        status=int(result["status"]),
        message=str(result["message"]),
        nfev=int(result["nfev"]),
        backend=str(result["backend"]),
        **_fit_derivative_metadata(result),
    )


def _fit_mea_co2_h2o_pure_parameter_benchmark(
    records: Any,
    *,
    dataset: str | Path,
    species: Iterable[str] | None = None,
    user_options: Mapping[str, Any] | None = None,
    initial_guess: Mapping[str, Mapping[str, float]] | None = None,
    bounds: FitBounds | Mapping[str, tuple[float | None, float | None]] | None = None,
    multistart: int = 0,
    max_nfev: int = 1,
) -> dict[str, FitResult]:
    """Internal opt-in benchmark hook for MEA-CO2-H2O pure-parameter fitting."""

    normalized_records = _normalize_records(records)
    normalized_species = _mea_co2_h2o_species_from_records(normalized_records, species)
    terms = _build_mea_co2_h2o_terms(normalized_records)
    components = ("MEA", "MEAH+", "MEACOO-", "HCO3-")
    T_ref = float(np.mean([_float_from_record(record, "T", required=True) for record in normalized_records]))
    seed_payloads, pure_file_hints = _benchmark_seed_payloads(dataset, normalized_species, T_ref, components)
    fit_targets = {
        "MEA": ("m", "s", "e", "e_assoc", "vol_a"),
        "MEAH+": ("s", "e", "d_born"),
        "MEACOO-": ("s", "e", "d_born"),
        "HCO3-": ("s", "e", "d_born"),
    }
    guesses = initial_guess or {}
    return {
        component: _fit_mea_co2_h2o_component(
            normalized_records,
            component,
            dataset=dataset,
            species=normalized_species,
            fit_targets=fit_targets[component],
            seed_payload=seed_payloads[component],
            pure_file_hint=pure_file_hints[component],
            terms=terms,
            initial_guess=guesses.get(component, {}),
            bounds=bounds,
            user_options=user_options,
            multistart=multistart,
            max_nfev=max_nfev,
        )
        for component in components
    }


def fit_mea_co2_h2o_electrolyte(
    records: Any,
    *,
    dataset: str | Path,
    species: Iterable[str] | None = None,
    user_options: Mapping[str, Any] | None = None,
    initial_guess: Mapping[str, Mapping[str, float]] | None = None,
    bounds: FitBounds | Mapping[str, tuple[float | None, float | None]] | None = None,
    multistart: int = 0,
    max_nfev: int = 1,
) -> dict[str, FitResult]:
    """Fit the dataset-driven MEA-CO2-H2O electrolyte pure-parameter benchmark.

    This helper exposes the opt-in benchmark path for MEA solvent systems while
    keeping records, chemistry, and dataset ownership outside the core runtime.
    It returns one :class:`FitResult` per fitted pure component and includes
    Born-diameter targets for the ionic species.
    """

    return _fit_mea_co2_h2o_pure_parameter_benchmark(
        records,
        dataset=dataset,
        species=species,
        user_options=user_options,
        initial_guess=initial_guess,
        bounds=bounds,
        multistart=multistart,
        max_nfev=max_nfev,
    )


def evaluate_reactive_electrolyte_bubble_residuals(
    records: Any,
    *,
    species: Sequence[str],
    mixture_factory: Any,
    balances: Mapping[str, Mapping[str, float]],
    reactions: Any,
    vapor_species: Sequence[str],
    volatile_species: Sequence[str] | None = None,
    nonvolatile_species: Sequence[str] | None = None,
    options: Any = None,
    continuation: str = "auto",
    pressure_species: Sequence[str] | None = None,
    speciation_species: Sequence[str] | None = None,
    reaction_names: Sequence[str] | None = None,
    pressure_weight: float = 1.0,
    speciation_weight: float = 1.0,
    reaction_weight: float = 0.0,
    penalty_value: float = 8.0,
) -> ReactiveElectrolyteRegressionResult:
    """Evaluate a fixed-shape coupled reactive-electrolyte residual vector.

    This is a regression objective helper, not an optimizer. Downstream projects
    own target selection, parameter updates, artifact promotion, and optimizer
    policy. The package owns the coupled native speciation plus fixed-liquid
    electrolyte bubble call sequence, continuation seeds, result-mode failure
    shaping, and compact per-record diagnostics.
    """

    from .reactive_electrolyte import ReactiveElectrolyteBubbleOptions
    from .reactive_regression import (
        ReactiveElectrolyteBatch,
        ReactiveElectrolyteBatchOptions,
        ReactiveElectrolyteRegressionContext,
        ReactiveElectrolyteRow,
        build_reactive_regression_objective,
    )

    normalized_records = _normalize_records(records)
    species_labels = tuple(str(label) for label in species)
    if not species_labels:
        raise InputError("species must include at least one label.")
    pressure_labels = _infer_reactive_target_labels(
        normalized_records,
        pressure_species,
        target_keys=("target_partial_pressures", "partial_pressures", "pressure_targets"),
        scalar_keys=("target_pressure", "target_pressure_Pa", "target_partial_pressure", "target_partial_pressure_Pa"),
    )
    speciation_labels = _infer_reactive_target_labels(
        normalized_records,
        speciation_species,
        target_keys=("target_x", "target_composition", "speciation_targets"),
        scalar_keys=(),
    )
    reaction_labels = tuple(str(name) for name in (reaction_names or ()))
    penalty = _positive_regression_weight(penalty_value, "penalty_value")
    continuation_mode = _normalize_regression_continuation(continuation)
    row_weight_scale = {
        "partial_pressure": float(pressure_weight),
        "speciation": float(speciation_weight),
        "reaction": float(reaction_weight if reaction_labels else 0.0),
    }
    rows = []
    for idx, record in enumerate(normalized_records):
        if "T" not in record or "totals" not in record:
            raise InputError("Each reactive electrolyte regression record requires T and totals.")
        pressure_targets = _reactive_target_mapping(
            record,
            pressure_labels,
            target_keys=("target_partial_pressures", "partial_pressures", "pressure_targets"),
            scalar_keys=(
                "target_pressure",
                "target_pressure_Pa",
                "target_partial_pressure",
                "target_partial_pressure_Pa",
            ),
        )
        speciation_targets = _reactive_target_mapping(
            record,
            speciation_labels,
            target_keys=("target_x", "target_composition", "speciation_targets"),
            scalar_keys=(),
        )
        rows.append(
            ReactiveElectrolyteRow(
                row_id=str(record.get("row_id", record.get("id", f"record_{idx}"))),
                T=float(record["T"]),
                P=float(record["P"]) if record.get("P") not in (None, "") else None,
                P_seed=float(record["P_seed"]) if record.get("P_seed") not in (None, "") else None,
                initial_x=record.get("initial_x"),
                balances=record.get("balances", balances),
                totals=record["totals"],
                reactions=_record_reactions(reactions, record),
                vapor_species=vapor_species,
                target_pressure=None,
                target_speciation=speciation_targets,
                target_partial_pressures=pressure_targets,
                source=str(record.get("source", "")),
                split=str(record.get("split", "")),
                metadata=record.get("metadata", {}),
                mode="bubble",
            )
        )
    batch = ReactiveElectrolyteBatch(
        species=species_labels,
        rows=tuple(rows),
        balances=balances,
        reactions=reactions,
        vapor_species=vapor_species,
        volatile_species=volatile_species,
        nonvolatile_species=nonvolatile_species,
        mixture_factory=mixture_factory,
        reactive_bubble_options=(
            options
            if isinstance(options, ReactiveElectrolyteBubbleOptions)
            else ReactiveElectrolyteBubbleOptions(error_mode="result")
        ),
        options=ReactiveElectrolyteBatchOptions(
            warm_start_rows=continuation_mode == "auto",
            warm_start_objective=continuation_mode == "auto",
            reuse_pressure_seeds=True,
            penalty_value=penalty,
            include_state_outputs=False,
        ),
    )
    objective = build_reactive_regression_objective(
        batch,
        residual_weights=row_weight_scale,
        failure_penalty=penalty,
    )
    context = ReactiveElectrolyteRegressionContext.from_batch(
        species=batch.species,
        rows=batch.rows,
        balances=batch.balances,
        reactions=batch.reactions,
        options=batch.options,
        objective=objective,
        vapor_species=batch.vapor_species,
        volatile_species=batch.volatile_species,
        nonvolatile_species=batch.nonvolatile_species,
        mixture_factory=batch.mixture_factory,
        reactive_bubble_options=batch.reactive_bubble_options,
    )
    if reaction_labels:
        context.reaction_names = reaction_labels
    batch_result = context.evaluate({})
    record_results = []
    for row_result in batch_result.row_results:
        record_results.append(
            {
                "row_id": row_result.row_id,
                "success": bool(row_result.success),
                "message": row_result.message,
                "P_total": None if row_result.pressure is None else float(row_result.pressure),
                "partial_pressures": _json_like_regression(row_result.partial_pressures),
                "x_liq": _json_like_regression(row_result.composition),
                "y_vap": _json_like_regression(row_result.y_vap),
                "named_reaction_residuals": _json_like_regression(row_result.named_reaction_residuals),
                "fugacity_residual_norm": float(
                    row_result.failure_diagnostics.get("bubble", {}).get("fugacity_residual_norm", 0.0)
                    if isinstance(row_result.failure_diagnostics.get("bubble"), Mapping)
                    else 0.0
                ),
                "state_failure_count": int(
                    row_result.failure_diagnostics.get("state_failure_count", 0)
                    if isinstance(row_result.failure_diagnostics, Mapping)
                    else 0
                ),
                "diagnostics": _json_like_regression(row_result.failure_diagnostics),
                "residual_count": len(row_result.residual_names),
            }
        )
    return ReactiveElectrolyteRegressionResult(
        residuals=batch_result.residuals,
        residual_names=batch_result.residual_names,
        record_results=tuple(record_results),
        success_count=batch_result.success_count,
        failure_count=batch_result.failure_count,
        diagnostics={
            "record_count": len(normalized_records),
            "pressure_species": pressure_labels,
            "speciation_species": speciation_labels,
            "reaction_names": reaction_labels,
            "continuation": continuation_mode,
            "fixed_shape": True,
            "batch_cache_stats": _json_like_regression(batch_result.cache_stats),
            "timing_summary": _json_like_regression(batch_result.timing_summary),
        },
    )


def _infer_reactive_target_labels(
    records: Sequence[Mapping[str, Any]],
    explicit: Sequence[str] | None,
    *,
    target_keys: Sequence[str],
    scalar_keys: Sequence[str],
) -> tuple[str, ...]:
    if explicit is not None:
        return tuple(str(label) for label in explicit)
    labels: list[str] = []
    for record in records:
        for key in target_keys:
            value = record.get(key)
            if isinstance(value, Mapping):
                labels.extend(str(label) for label in value)
        if any(key in record for key in scalar_keys):
            label = str(record.get("target_species", record.get("target_pressure_species", "")))
            if label:
                labels.append(label)
    return tuple(dict.fromkeys(labels))


def _reactive_target_mapping(
    record: Mapping[str, Any],
    labels: Sequence[str],
    *,
    target_keys: Sequence[str],
    scalar_keys: Sequence[str],
) -> dict[str, float]:
    mapping: dict[str, float] = {}
    for key in target_keys:
        value = record.get(key)
        if isinstance(value, Mapping):
            mapping.update({str(label): float(target) for label, target in value.items()})
    scalar_value = None
    for key in scalar_keys:
        if key in record:
            scalar_value = record[key]
            break
    if scalar_value is not None:
        label = str(record.get("target_species", record.get("target_pressure_species", "")))
        if not label and len(labels) == 1:
            label = str(labels[0])
        if label:
            mapping[label] = float(scalar_value)
    missing = [label for label in labels if label not in mapping]
    if missing:
        raise InputError(f"Reactive electrolyte regression record is missing target(s): {', '.join(missing)}.")
    return {str(label): float(mapping[str(label)]) for label in labels}


def _positive_regression_weight(value: float, label: str) -> float:
    out = float(value)
    if not math.isfinite(out) or out < 0.0:
        raise InputError(f"{label} must be finite and non-negative.")
    return out


def _safe_log10_ratio(predicted: float, target: float, penalty: float) -> float:
    if not (math.isfinite(predicted) and math.isfinite(target)) or predicted <= 0.0 or target <= 0.0:
        return penalty
    return math.log10(max(predicted, 1.0e-300) / max(target, 1.0e-300))


def _normalize_regression_continuation(value: Any) -> str:
    if isinstance(value, bool):
        return "auto" if value else "none"
    token = str(value).strip().lower()
    if token in {"auto", "on", "true", "1"}:
        return "auto"
    if token in {"none", "off", "false", "0", "disabled"}:
        return "none"
    raise InputError("continuation must be 'auto' or 'none'.")


def _reactive_regression_options(options: Any, penalty: float, options_type: Any) -> Any:
    if options is None:
        return options_type(error_mode="result", penalty_value=penalty)
    if not isinstance(options, options_type):
        raise InputError("options must be a ReactiveElectrolyteBubbleOptions instance.")
    if options.error_mode != "result" or float(options.penalty_value) != penalty:
        return options_type(
            speciation_options=options.speciation_options,
            bubble_options=options.bubble_options,
            phase_handoff_mass_tolerance=options.phase_handoff_mass_tolerance,
            phase_handoff_charge_tolerance=options.phase_handoff_charge_tolerance,
            phase_handoff_reaction_tolerance=options.phase_handoff_reaction_tolerance,
            error_mode="result",
            penalty_value=penalty,
        )
    return options


def _reactive_regression_options_with_seed(
    options: Any,
    bubble_options_type: Any,
    *,
    initial_pressure: float,
    initial_y_vap: Mapping[str, float] | None,
    fallback_pressure: float,
) -> Any:
    base = options.bubble_options or bubble_options_type(
        initial_pressure=float(fallback_pressure),
        return_best_effort=True,
    )
    bubble_options = bubble_options_type(
        initial_pressure=float(initial_pressure),
        min_pressure=base.min_pressure,
        max_pressure=base.max_pressure,
        max_iterations=base.max_iterations,
        max_vapor_iterations=base.max_vapor_iterations,
        max_bracket_expansions=base.max_bracket_expansions,
        tolerance=base.tolerance,
        vapor_tolerance=base.vapor_tolerance,
        pressure_factor=base.pressure_factor,
        min_composition=base.min_composition,
        charge_tolerance=base.charge_tolerance,
        return_best_effort=base.return_best_effort,
        initial_y_vap=initial_y_vap or base.initial_y_vap,
    )
    return type(options)(
        speciation_options=options.speciation_options,
        bubble_options=bubble_options,
        phase_handoff_mass_tolerance=options.phase_handoff_mass_tolerance,
        phase_handoff_charge_tolerance=options.phase_handoff_charge_tolerance,
        phase_handoff_reaction_tolerance=options.phase_handoff_reaction_tolerance,
        error_mode=options.error_mode,
        penalty_value=options.penalty_value,
    )


def _record_reactions(reactions: Any, record: Mapping[str, Any]) -> Any:
    if not callable(reactions):
        return reactions
    try:
        return reactions(record)
    except TypeError:
        return reactions(float(record["T"]))


def write_fit_result(
    result: FitResult,
    dataset_root: str | Path,
    *,
    overwrite: bool = False,
    pure_file: str | Path | None = None,
) -> list[Path]:
    """Write fitted values into a user-owned dataset folder."""

    root = Path(dataset_root).expanduser()
    if not root.exists():
        raise FileNotFoundError(f"Dataset folder not found: {root}")
    written: list[Path] = []
    problem = result.problem

    if problem.mode in {PURE_NEUTRAL_MODE, PURE_ION_MODE}:
        path = _choose_pure_file(problem, root, pure_file)
        updates = {target: result.rendered_values[target] for target in problem.fit_targets}
        _update_csv_row(path, str(problem.component), updates, overwrite=overwrite)
        _invalidate_dataset_cache(root)
        written.append(path)
        return written

    if problem.mode == BINARY_PAIR_MODE:
        bi_dir = root / "mixed" / "binary_interaction"
        if not bi_dir.exists():
            raise FileNotFoundError(f"Dataset folder '{root}' does not contain mixed/binary_interaction/.")
        if problem.pair is None:
            raise InputError("binary_pair fit results require problem.pair before writing.")
        for target in problem.fit_targets:
            path = bi_dir / MATRIX_FILE_NAMES[target]
            if not path.exists():
                raise FileNotFoundError(f"Expected matrix file '{path}' to exist.")
            _update_matrix_file(path, problem.pair, result.rendered_values[target], overwrite=overwrite)
            written.append(path)
        _invalidate_dataset_cache(root)
        return written

    raise InputError(f"Unsupported fit result mode '{problem.mode}'.")
