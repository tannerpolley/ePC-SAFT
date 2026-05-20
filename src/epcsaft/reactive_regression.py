"""Batched reactive-electrolyte regression helpers and diagnostics."""

from __future__ import annotations

import csv
import json
import math
import time
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from . import reactive_electrolyte as reactive_electrolyte_module
from . import reactive_speciation as reactive_speciation_module
from ._types import InputError
from .epcsaft import ePCSAFTMixture
from .parameter_schema import ParameterSet, ParameterSource
from .reactive_electrolyte import ReactiveElectrolyteBubbleOptions, ReactiveElectrolyteBubbleResult
from .reactive_speciation import ReactionDefinition, ReactiveSpeciationOptions, ReactiveSpeciationResult
from .regression import _compile_target_family_summaries

_PRESSURE_FAMILY_ALIASES = {
    "pressure_linear": "pressure_linear",
    "pressure_log": "pressure_log",
    "partial_pressure_log": "partial_pressure_log",
}
_SPECIATION_FAMILY_ALIASES = {
    "speciation_mole_fraction": "speciation_mole_fraction",
    "speciation_log_mole_fraction": "speciation_log_mole_fraction",
}
_ACTIVITY_FAMILY_ALIASES = {
    "activity_coefficient": "activity_coefficient",
    "log_activity_coefficient": "log_activity_coefficient",
}
_FUGACITY_FAMILY_ALIASES = {
    "fugacity_coefficient": "fugacity_coefficient",
    "ln_fugacity_coefficient": "ln_fugacity_coefficient",
}
_SCALAR_FAMILY_ALIASES = {
    "density": "density",
    "relative_permittivity": "relative_permittivity",
}

_PARAMETER_FIELD_ALIASES = {
    "m": "m",
    "segment_number": "m",
    "sigma": "s",
    "s": "s",
    "epsilon_k": "e",
    "epsilon_over_k": "e",
    "e": "e",
    "mw": "MW",
    "molecular_weight": "MW",
    "z": "z",
    "dielc": "dielc",
    "relative_permittivity": "dielc",
    "d_born": "d_born",
    "born_diameter": "d_born",
    "f_solv": "f_solv",
    "solvation_factor": "f_solv",
    "e_assoc": "e_assoc",
    "association_energy": "e_assoc",
    "vol_a": "vol_a",
    "association_volume": "vol_a",
}
_BINARY_FIELD_ALIASES = {
    "k_ij": "k_ij",
    "l_ij": "l_ij",
    "k_hb_ij": "k_hb_ij",
}


def _json_like(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return [_json_like(item) for item in value.tolist()]
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, Mapping):
        return {str(k): _json_like(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_like(item) for item in value]
    return value


def _positive(value: float, label: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise InputError(f"{label} must be finite and positive.")
    return number


def _normalize_log_ratio(predicted: float, target: float, penalty: float) -> float:
    if predicted <= 0.0 or target <= 0.0 or not math.isfinite(predicted) or not math.isfinite(target):
        return float(penalty)
    return math.log10(predicted / target)


def _mapping_float(values: Mapping[str, Any] | None) -> dict[str, float]:
    if values is None:
        return {}
    return {str(k): float(v) for k, v in values.items()}


def _weights_float(values: Mapping[str, Any] | None, *, allow_zero: bool = False) -> dict[str, float]:
    if values is None:
        return {}
    out: dict[str, float] = {}
    for key, raw in values.items():
        number = float(raw)
        if not math.isfinite(number) or number < 0.0 or (not allow_zero and number <= 0.0):
            adjective = "non-negative" if allow_zero else "positive"
            raise InputError(f"weights[{key!r}] must be finite and {adjective}.")
        out[str(key)] = number
    return out


def _normalize_reaction_defs(reactions: Sequence[Any] | None) -> tuple[ReactionDefinition, ...]:
    if not reactions:
        return ()
    out: list[ReactionDefinition] = []
    for reaction in reactions:
        if isinstance(reaction, ReactionDefinition):
            out.append(reaction)
            continue
        if isinstance(reaction, Mapping):
            out.append(
                ReactionDefinition(
                    stoichiometry=reaction["stoichiometry"],
                    log_equilibrium_constant=float(reaction["log_equilibrium_constant"]),
                    name=str(reaction.get("name", "")),
                    standard_state=str(reaction.get("standard_state", "mole_fraction_activity")),
                )
            )
            continue
        raise InputError("reactions must contain ReactionDefinition instances or mapping-like definitions.")
    return tuple(out)


def _sequence_str(values: Sequence[Any] | None) -> tuple[str, ...]:
    if values is None:
        return ()
    return tuple(str(value) for value in values)


@dataclass(frozen=True, slots=True)
class ReactiveElectrolyteRow:
    """One reactive-electrolyte regression row."""

    row_id: str
    T: float
    P: float
    initial_x: Sequence[float] | Mapping[str, float] | None = None
    balances: Mapping[str, Mapping[str, float]] | None = None
    totals: Mapping[str, float] | None = None
    reactions: Sequence[ReactionDefinition] | Sequence[Mapping[str, Any]] | None = None
    vapor_species: Sequence[str] | None = None
    target_pressure: float | None = None
    target_speciation: Mapping[str, float] | None = None
    target_activity: Mapping[str, float] | None = None
    target_fugacity: Mapping[str, float] | None = None
    target_density: float | None = None
    target_relative_permittivity: float | None = None
    target_partial_pressures: Mapping[str, float] | None = None
    weights: Mapping[str, float] = field(default_factory=dict)
    source: str = ""
    split: str = ""
    metadata: Mapping[str, Any] = field(default_factory=dict)
    mode: str = "bubble"

    def __post_init__(self) -> None:
        object.__setattr__(self, "row_id", str(self.row_id))
        object.__setattr__(self, "T", float(self.T))
        if not math.isfinite(self.T) or self.T <= 0.0:
            raise InputError("ReactiveElectrolyteRow.T must be finite and positive.")
        object.__setattr__(self, "P", float(self.P))
        if not math.isfinite(self.P) or self.P <= 0.0:
            raise InputError("ReactiveElectrolyteRow.P must be finite and positive.")
        object.__setattr__(self, "balances", dict(self.balances or {}))
        object.__setattr__(self, "totals", _mapping_float(self.totals))
        object.__setattr__(self, "reactions", _normalize_reaction_defs(self.reactions))
        object.__setattr__(self, "vapor_species", _sequence_str(self.vapor_species))
        object.__setattr__(self, "target_speciation", _mapping_float(self.target_speciation))
        object.__setattr__(self, "target_activity", _mapping_float(self.target_activity))
        object.__setattr__(self, "target_fugacity", _mapping_float(self.target_fugacity))
        object.__setattr__(self, "target_partial_pressures", _mapping_float(self.target_partial_pressures))
        object.__setattr__(self, "weights", _weights_float(self.weights))
        object.__setattr__(self, "source", str(self.source))
        object.__setattr__(self, "split", str(self.split))
        object.__setattr__(self, "metadata", dict(self.metadata))
        mode = str(self.mode).strip().lower()
        if mode not in {"bubble", "speciation"}:
            raise InputError("ReactiveElectrolyteRow.mode must be 'bubble' or 'speciation'.")
        object.__setattr__(self, "mode", mode)


@dataclass(frozen=True, slots=True)
class ReactiveElectrolyteBatchOptions:
    """Batch-evaluation policy controls."""

    penalty_value: float = 8.0
    failure_residual_mode: str = "penalty"
    include_state_outputs: bool = True

    def __post_init__(self) -> None:
        mode = str(self.failure_residual_mode).strip().lower()
        if mode not in {"penalty", "drop"}:
            raise InputError("ReactiveElectrolyteBatchOptions.failure_residual_mode must be 'penalty' or 'drop'.")
        object.__setattr__(self, "failure_residual_mode", mode)
        object.__setattr__(self, "penalty_value", _positive(float(self.penalty_value), "penalty_value"))


@dataclass(frozen=True, slots=True)
class ReactiveRegressionObjective:
    """Controls residual packing for reactive-electrolyte regression."""

    pressure_family: str = "pressure_log"
    speciation_family: str = "speciation_log_mole_fraction"
    activity_family: str = "log_activity_coefficient"
    fugacity_family: str = "ln_fugacity_coefficient"
    scalar_families: tuple[str, ...] = ("density", "relative_permittivity")
    row_weights: Mapping[str, float] = field(default_factory=dict)
    source_weights: Mapping[str, float] = field(default_factory=dict)
    split_weights: Mapping[str, float] = field(default_factory=dict)
    residual_weights: Mapping[str, float] = field(default_factory=dict)
    failure_penalty: float = 8.0
    residual_clip: float | None = None

    def __post_init__(self) -> None:
        pressure_family = _PRESSURE_FAMILY_ALIASES.get(str(self.pressure_family).strip().lower())
        if pressure_family is None:
            raise InputError("ReactiveRegressionObjective.pressure_family is not supported.")
        speciation_family = _SPECIATION_FAMILY_ALIASES.get(str(self.speciation_family).strip().lower())
        if speciation_family is None:
            raise InputError("ReactiveRegressionObjective.speciation_family is not supported.")
        activity_family = _ACTIVITY_FAMILY_ALIASES.get(str(self.activity_family).strip().lower())
        if activity_family is None:
            raise InputError("ReactiveRegressionObjective.activity_family is not supported.")
        fugacity_family = _FUGACITY_FAMILY_ALIASES.get(str(self.fugacity_family).strip().lower())
        if fugacity_family is None:
            raise InputError("ReactiveRegressionObjective.fugacity_family is not supported.")
        scalar = []
        for name in self.scalar_families:
            token = _SCALAR_FAMILY_ALIASES.get(str(name).strip().lower())
            if token is None:
                raise InputError(f"ReactiveRegressionObjective scalar family {name!r} is not supported.")
            scalar.append(token)
        object.__setattr__(self, "pressure_family", pressure_family)
        object.__setattr__(self, "speciation_family", speciation_family)
        object.__setattr__(self, "activity_family", activity_family)
        object.__setattr__(self, "fugacity_family", fugacity_family)
        object.__setattr__(self, "scalar_families", tuple(scalar))
        object.__setattr__(self, "row_weights", _weights_float(self.row_weights))
        object.__setattr__(self, "source_weights", _weights_float(self.source_weights))
        object.__setattr__(self, "split_weights", _weights_float(self.split_weights))
        object.__setattr__(self, "residual_weights", _weights_float(self.residual_weights, allow_zero=True))
        object.__setattr__(self, "failure_penalty", _positive(float(self.failure_penalty), "failure_penalty"))
        if self.residual_clip is not None:
            object.__setattr__(self, "residual_clip", _positive(float(self.residual_clip), "residual_clip"))


@dataclass(frozen=True, slots=True)
class ReactiveElectrolyteBatch:
    """Batch input for repeated reactive-electrolyte evaluations."""

    species: Sequence[str]
    rows: Sequence[ReactiveElectrolyteRow]
    balances: Mapping[str, Mapping[str, float]]
    reactions: Sequence[ReactionDefinition] | Sequence[Mapping[str, Any]]
    vapor_species: Sequence[str] | None = None
    volatile_species: Sequence[str] | None = None
    nonvolatile_species: Sequence[str] | None = None
    base_parameters: Mapping[str, Any] | None = None
    user_options: Mapping[str, Any] | None = None
    mixture_factory: Any | None = None
    mixture_factory_builder: Any | None = None
    reactive_speciation_options: ReactiveSpeciationOptions | None = None
    reactive_bubble_options: ReactiveElectrolyteBubbleOptions | None = None
    options: ReactiveElectrolyteBatchOptions = field(default_factory=ReactiveElectrolyteBatchOptions)

    def __post_init__(self) -> None:
        species = _sequence_str(self.species)
        if not species:
            raise InputError("ReactiveElectrolyteBatch.species must include at least one species.")
        rows = tuple(self.rows)
        if not rows:
            raise InputError("ReactiveElectrolyteBatch.rows must include at least one row.")
        if self.base_parameters is None and self.mixture_factory is None and self.mixture_factory_builder is None:
            raise InputError(
                "ReactiveElectrolyteBatch requires base_parameters, mixture_factory, or mixture_factory_builder."
            )
        object.__setattr__(self, "species", species)
        object.__setattr__(self, "rows", rows)
        object.__setattr__(self, "balances", dict(self.balances))
        object.__setattr__(self, "reactions", _normalize_reaction_defs(self.reactions))
        object.__setattr__(self, "vapor_species", _sequence_str(self.vapor_species))
        object.__setattr__(self, "volatile_species", _sequence_str(self.volatile_species))
        object.__setattr__(self, "nonvolatile_species", _sequence_str(self.nonvolatile_species))
        user_options = _copy_parameter_payload(self.user_options or {})
        if self.base_parameters is not None:
            object.__setattr__(
                self,
                "base_parameters",
                ParameterSource(self.base_parameters, species=species).to_runtime_dict(user_options=user_options),
            )
        object.__setattr__(self, "user_options", user_options)


@dataclass(frozen=True, slots=True)
class ReactiveElectrolyteRowResult:
    """Structured result for one regression row."""

    row_id: str
    success: bool
    message: str
    composition: Mapping[str, float]
    pressure: float | None
    ln_fugacity: Mapping[str, float]
    activity_coefficients: Mapping[str, float]
    density: float | None
    relative_permittivity: float | None
    residuals: Mapping[str, float]
    residual_names: tuple[str, ...]
    failure_diagnostics: Mapping[str, Any]
    active_bounds: Mapping[str, bool]
    elapsed_seconds: float
    partial_pressures: Mapping[str, float] = field(default_factory=dict)
    y_vap: Mapping[str, float] = field(default_factory=dict)
    named_reaction_residuals: Mapping[str, float] = field(default_factory=dict)
    source: str = ""
    split: str = ""
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "row_id": self.row_id,
            "success": bool(self.success),
            "message": self.message,
            "composition": _json_like(self.composition),
            "pressure": None if self.pressure is None else float(self.pressure),
            "ln_fugacity": _json_like(self.ln_fugacity),
            "activity_coefficients": _json_like(self.activity_coefficients),
            "density": None if self.density is None else float(self.density),
            "relative_permittivity": (
                None if self.relative_permittivity is None else float(self.relative_permittivity)
            ),
            "residuals": _json_like(self.residuals),
            "residual_names": list(self.residual_names),
            "failure_diagnostics": _json_like(self.failure_diagnostics),
            "active_bounds": _json_like(self.active_bounds),
            "elapsed_seconds": float(self.elapsed_seconds),
            "partial_pressures": _json_like(self.partial_pressures),
            "y_vap": _json_like(self.y_vap),
            "named_reaction_residuals": _json_like(self.named_reaction_residuals),
            "source": self.source,
            "split": self.split,
            "metadata": _json_like(self.metadata),
        }


@dataclass(frozen=True, slots=True)
class ReactiveElectrolyteBatchResult:
    """Structured result for one batch evaluation."""

    success_count: int
    failure_count: int
    row_results: tuple[ReactiveElectrolyteRowResult, ...]
    residuals: np.ndarray
    residual_names: tuple[str, ...]
    residual_row_map: tuple[str, ...]
    diagnostics: Mapping[str, Any]
    cache_stats: Mapping[str, Any]
    timing_summary: Mapping[str, Any]

    def __post_init__(self) -> None:
        object.__setattr__(self, "residuals", np.asarray(self.residuals, dtype=float).reshape(-1))
        if len(self.residual_names) != int(self.residuals.size):
            raise InputError("ReactiveElectrolyteBatchResult residual_names must match residual length.")
        if len(self.residual_row_map) != int(self.residuals.size):
            raise InputError("ReactiveElectrolyteBatchResult residual_row_map must match residual length.")

    def to_dict(self) -> dict[str, Any]:
        return {
            "success_count": int(self.success_count),
            "failure_count": int(self.failure_count),
            "row_results": [row.to_dict() for row in self.row_results],
            "residuals": [float(value) for value in self.residuals],
            "residual_names": list(self.residual_names),
            "residual_row_map": list(self.residual_row_map),
            "diagnostics": _json_like(self.diagnostics),
            "cache_stats": _json_like(self.cache_stats),
            "timing_summary": _json_like(self.timing_summary),
        }

    def flatten_rows(self) -> list[dict[str, Any]]:
        rows: list[dict[str, Any]] = []
        for row in self.row_results:
            item = row.to_dict()
            item["residual_norm"] = float(
                np.linalg.norm(np.asarray(list(row.residuals.values()), dtype=float)) if row.residuals else 0.0
            )
            rows.append(item)
        return rows

    def flatten_residuals(self) -> list[dict[str, Any]]:
        rows: list[dict[str, Any]] = []
        for name, row_id, value in zip(self.residual_names, self.residual_row_map, self.residuals):
            rows.append({"row_id": row_id, "residual_name": name, "residual_value": float(value)})
        return rows


@dataclass(frozen=True, slots=True)
class ReactiveRegressionObjectiveResult:
    """Objective-level view over a batch evaluation."""

    batch_result: ReactiveElectrolyteBatchResult
    residuals: np.ndarray
    residual_names: tuple[str, ...]
    residual_row_map: tuple[str, ...]
    metrics: Mapping[str, Any]
    diagnostics: Mapping[str, Any]

    def __post_init__(self) -> None:
        object.__setattr__(self, "residuals", np.asarray(self.residuals, dtype=float).reshape(-1))

    @property
    def objective(self) -> float:
        return 0.5 * float(np.dot(self.residuals, self.residuals))

    def to_dict(self) -> dict[str, Any]:
        return {
            "objective": self.objective,
            "residuals": [float(value) for value in self.residuals],
            "residual_names": list(self.residual_names),
            "residual_row_map": list(self.residual_row_map),
            "metrics": _json_like(self.metrics),
            "diagnostics": _json_like(self.diagnostics),
            "batch_result": self.batch_result.to_dict(),
        }


@dataclass(frozen=True, slots=True)
class ReactiveRegressionJacobianResult:
    """Derivative payload for reactive regression paths with implemented sensitivities."""

    jacobian: np.ndarray
    gradient: np.ndarray
    residuals: np.ndarray
    residual_names: tuple[str, ...]
    residual_row_map: tuple[str, ...]
    parameter_names: tuple[str, ...]
    diagnostics: Mapping[str, Any]

    def __post_init__(self) -> None:
        object.__setattr__(self, "jacobian", np.asarray(self.jacobian, dtype=float))
        object.__setattr__(self, "gradient", np.asarray(self.gradient, dtype=float).reshape(-1))
        object.__setattr__(self, "residuals", np.asarray(self.residuals, dtype=float).reshape(-1))

    def to_dict(self) -> dict[str, Any]:
        return {
            "jacobian": _json_like(self.jacobian),
            "jacobian_shape": list(self.jacobian.shape),
            "gradient": _json_like(self.gradient),
            "residuals": _json_like(self.residuals),
            "residual_names": list(self.residual_names),
            "residual_row_map": list(self.residual_row_map),
            "parameter_names": list(self.parameter_names),
            "diagnostics": _json_like(self.diagnostics),
        }


@dataclass(slots=True)
class ReactiveElectrolyteRegressionContext:
    """Compiled reusable context for repeated reactive regression evaluations."""

    batch: ReactiveElectrolyteBatch
    objective: ReactiveRegressionObjective
    species: tuple[str, ...]
    species_index: dict[str, int]
    charge_vector: np.ndarray
    neutral_mask: np.ndarray
    ion_mask: np.ndarray
    reaction_names: tuple[str, ...]
    balance_names: tuple[str, ...]
    residual_name_template: tuple[str, ...]
    row_static_metadata: tuple[dict[str, Any], ...]
    context_diagnostics: dict[str, Any]
    _stats: Counter = field(default_factory=Counter, repr=False)

    @classmethod
    def from_batch(
        cls,
        *,
        species: Sequence[str],
        rows: Sequence[ReactiveElectrolyteRow],
        balances: Mapping[str, Mapping[str, float]],
        reactions: Sequence[ReactionDefinition] | Sequence[Mapping[str, Any]],
        options: ReactiveElectrolyteBatchOptions | None = None,
        objective: ReactiveRegressionObjective | None = None,
        vapor_species: Sequence[str] | None = None,
        volatile_species: Sequence[str] | None = None,
        nonvolatile_species: Sequence[str] | None = None,
        base_parameters: Mapping[str, Any] | None = None,
        user_options: Mapping[str, Any] | None = None,
        mixture_factory: Any | None = None,
        mixture_factory_builder: Any | None = None,
        reactive_speciation_options: ReactiveSpeciationOptions | None = None,
        reactive_bubble_options: ReactiveElectrolyteBubbleOptions | None = None,
    ) -> ReactiveElectrolyteRegressionContext:
        batch = ReactiveElectrolyteBatch(
            species=species,
            rows=rows,
            balances=balances,
            reactions=reactions,
            vapor_species=vapor_species,
            volatile_species=volatile_species,
            nonvolatile_species=nonvolatile_species,
            base_parameters=base_parameters,
            user_options=user_options,
            mixture_factory=mixture_factory,
            mixture_factory_builder=mixture_factory_builder,
            reactive_speciation_options=reactive_speciation_options,
            reactive_bubble_options=reactive_bubble_options,
            options=options or ReactiveElectrolyteBatchOptions(),
        )
        resolved_objective = objective or build_reactive_regression_objective(batch)
        species_labels = tuple(batch.species)
        species_index = {label: idx for idx, label in enumerate(species_labels)}
        reference_params = None
        if batch.base_parameters is not None:
            reference_params = batch.base_parameters
        else:
            row0 = batch.rows[0]
            mixture = _build_row_mixture(
                batch, row0, {}, x_override=_row_initial_x(row0, species_labels), P_override=row0.P or 101325.0
            )
            reference_params = getattr(mixture, "_params", None)
        if reference_params is None:
            charge_vector = np.zeros(len(species_labels), dtype=float)
        else:
            charge_vector = np.asarray(reference_params.get("z", np.zeros(len(species_labels))), dtype=float).reshape(
                -1
            )
            if charge_vector.size != len(species_labels):
                charge_vector = np.zeros(len(species_labels), dtype=float)
        reaction_names = _reaction_names(batch.reactions)
        residual_template = _compile_residual_name_template(batch.rows, reaction_names)
        target_family_schema = _compile_target_family_summaries(
            Counter(_residual_family_from_name(name) for name in residual_template),
            count_label="residual_count",
        )
        row_static = [
            {
                "row_id": row.row_id,
                "source": row.source,
                "split": row.split,
                "mode": row.mode,
                "target_keys": sorted(
                    set(row.target_speciation)
                    | set(row.target_activity)
                    | set(row.target_fugacity)
                    | set(row.target_partial_pressures)
                ),
            }
            for row in batch.rows
        ]
        return cls(
            batch=batch,
            objective=resolved_objective,
            species=species_labels,
            species_index=species_index,
            charge_vector=charge_vector,
            neutral_mask=np.abs(charge_vector) <= 1.0e-12,
            ion_mask=np.abs(charge_vector) > 1.0e-12,
            reaction_names=reaction_names,
            balance_names=tuple(str(name) for name in batch.balances),
            residual_name_template=residual_template,
            row_static_metadata=tuple(row_static),
            context_diagnostics={
                "species_order": list(species_labels),
                "species_index": dict(species_index),
                "charge_vector": _json_like(charge_vector),
                "neutral_mask": _json_like(np.abs(charge_vector) <= 1.0e-12),
                "ion_mask": _json_like(np.abs(charge_vector) > 1.0e-12),
                "reaction_names": list(reaction_names),
                "balance_names": [str(name) for name in batch.balances],
                "row_count": len(batch.rows),
                "target_schema_size": len(residual_template),
                "target_family_schema": target_family_schema,
            },
        )

    def evaluate(
        self,
        parameter_map: Mapping[str, float] | None = None,
    ) -> ReactiveElectrolyteBatchResult:
        parameter_map = {str(k): float(v) for k, v in (parameter_map or {}).items()}
        start_ns = time.perf_counter_ns()
        row_results: list[ReactiveElectrolyteRowResult] = []
        residuals: list[float] = []
        residual_names: list[str] = []
        residual_row_map: list[str] = []
        success_count = 0
        failure_count = 0
        target_counter: Counter[str] = Counter()
        residual_values_by_family: dict[str, list[float]] = {}
        solve_counter: Counter[str] = Counter()
        local_stats: Counter[str] = Counter()

        for row in self.batch.rows:
            row_start_ns = time.perf_counter_ns()
            try:
                pressure = float(row.P)
                if row.mode == "speciation":
                    solve_counter["speciation_solves"] += 1
                    raw_result = reactive_speciation_module.solve_reactive_speciation(
                        species=self.species,
                        mixture_factory=_mixture_factory_from_batch(self.batch, parameter_map),
                        T=float(row.T),
                        P=pressure,
                        balances=row.balances or self.batch.balances,
                        totals=row.totals,
                        reactions=row.reactions or self.batch.reactions,
                        initial_x=_row_initial_x(row, self.species),
                        options=self.batch.reactive_speciation_options,
                    )
                    row_result = _row_result_from_speciation(
                        context=self,
                        row=row,
                        raw_result=raw_result,
                        elapsed_seconds=(time.perf_counter_ns() - row_start_ns) / 1.0e9,
                    )
                else:
                    solve_counter["bubble_solves"] += 1
                    raw_result = reactive_electrolyte_module.solve_reactive_electrolyte_bubble(
                        species=self.species,
                        mixture_factory=_mixture_factory_from_batch(self.batch, parameter_map),
                        T=float(row.T),
                        P=pressure,
                        balances=row.balances or self.batch.balances,
                        totals=row.totals,
                        reactions=row.reactions or self.batch.reactions,
                        initial_x=_row_initial_x(row, self.species),
                        vapor_species=row.vapor_species or self.batch.vapor_species,
                        volatile_species=self.batch.volatile_species,
                        nonvolatile_species=self.batch.nonvolatile_species,
                        options=self.batch.reactive_bubble_options,
                    )
                    row_result = _row_result_from_bubble(
                        context=self,
                        row=row,
                        raw_result=raw_result,
                        parameter_map=parameter_map,
                        elapsed_seconds=(time.perf_counter_ns() - row_start_ns) / 1.0e9,
                    )
                if row_result.success:
                    success_count += 1
                else:
                    failure_count += 1
                row_results.append(row_result)
            except Exception as exc:
                failure_count += 1
                row_results.append(
                    _failed_row_result(
                        row=row,
                        exc=exc,
                        elapsed_seconds=(time.perf_counter_ns() - row_start_ns) / 1.0e9,
                    )
                )

            packed_names, packed_values = _pack_row_residuals(
                row=row,
                row_result=row_results[-1],
                objective=self.objective,
                reaction_names=self.reaction_names,
                penalty=self.batch.options.penalty_value,
            )
            for name, value in zip(packed_names, packed_values):
                family = _residual_family_from_name(name)
                residual_names.append(name)
                residual_row_map.append(row.row_id)
                residual_value = float(value)
                residuals.append(residual_value)
                target_counter[family] += 1
                residual_values_by_family.setdefault(family, []).append(residual_value)
            local_stats["row_evaluations"] += 1

        self._stats["context_evaluations"] += 1
        self._stats.update(local_stats)
        total_elapsed = (time.perf_counter_ns() - start_ns) / 1.0e9
        timing_summary = {
            "elapsed_seconds": total_elapsed,
            "mean_row_seconds": total_elapsed / max(len(self.batch.rows), 1),
        }
        cache_stats = {
            "context_evaluations": int(self._stats["context_evaluations"]),
            "row_evaluations": int(self._stats["row_evaluations"]),
        }
        target_family_summaries = _compile_target_family_summaries(
            target_counter,
            {
                family: float(np.linalg.norm(np.asarray(values, dtype=float)))
                for family, values in residual_values_by_family.items()
            },
            count_label="residual_count",
        )
        diagnostics = {
            "parameter_map_keys": sorted(parameter_map),
            "reaction_names": list(self.reaction_names),
            "target_family_counts": dict(target_counter),
            "target_family_summaries": target_family_summaries,
            "solve_counts": dict(solve_counter),
            "context": _json_like(self.context_diagnostics),
        }
        return ReactiveElectrolyteBatchResult(
            success_count=success_count,
            failure_count=failure_count,
            row_results=tuple(row_results),
            residuals=np.asarray(residuals, dtype=float),
            residual_names=tuple(residual_names),
            residual_row_map=tuple(residual_row_map),
            diagnostics=diagnostics,
            cache_stats=cache_stats,
            timing_summary=timing_summary,
        )

    def evaluate_objective(
        self,
        parameter_map: Mapping[str, float] | None = None,
    ) -> ReactiveRegressionObjectiveResult:
        batch_result = self.evaluate(parameter_map)
        metrics = _objective_metrics(batch_result)
        return ReactiveRegressionObjectiveResult(
            batch_result=batch_result,
            residuals=batch_result.residuals,
            residual_names=batch_result.residual_names,
            residual_row_map=batch_result.residual_row_map,
            metrics=metrics,
            diagnostics={
                "objective": asdict(self.objective),
                "cache_stats": _json_like(batch_result.cache_stats),
            },
        )

    def evaluate_derivatives(
        self,
        parameter_map: Mapping[str, float],
        *,
        parameters: Sequence[str],
    ) -> ReactiveRegressionJacobianResult:
        _ = parameter_map, parameters
        raise InputError("reactive-regression sensitivities require native Ceres derivative coverage.")


# AlgID: reactive_electrolyte_batch_residual_context
def build_reactive_regression_objective(
    batch: ReactiveElectrolyteBatch,
    *,
    pressure_family: str = "pressure_log",
    speciation_family: str = "speciation_log_mole_fraction",
    activity_family: str = "log_activity_coefficient",
    fugacity_family: str = "ln_fugacity_coefficient",
    scalar_families: Sequence[str] = ("density", "relative_permittivity"),
    row_weights: Mapping[str, float] | None = None,
    source_weights: Mapping[str, float] | None = None,
    split_weights: Mapping[str, float] | None = None,
    residual_weights: Mapping[str, float] | None = None,
    failure_penalty: float | None = None,
    residual_clip: float | None = None,
) -> ReactiveRegressionObjective:
    return ReactiveRegressionObjective(
        pressure_family=pressure_family,
        speciation_family=speciation_family,
        activity_family=activity_family,
        fugacity_family=fugacity_family,
        scalar_families=tuple(scalar_families),
        row_weights=dict(row_weights or {}),
        source_weights=dict(source_weights or {}),
        split_weights=dict(split_weights or {}),
        residual_weights=dict(residual_weights or {}),
        failure_penalty=float(failure_penalty if failure_penalty is not None else batch.options.penalty_value),
        residual_clip=residual_clip,
    )


# AlgID: reactive_electrolyte_batch_residual_context
def evaluate_reactive_regression_objective(
    batch_or_context: ReactiveElectrolyteBatch | ReactiveElectrolyteRegressionContext,
    *,
    parameter_map: Mapping[str, float] | None = None,
    objective: ReactiveRegressionObjective | None = None,
) -> ReactiveRegressionObjectiveResult:
    if isinstance(batch_or_context, ReactiveElectrolyteRegressionContext):
        context = batch_or_context
    else:
        context = ReactiveElectrolyteRegressionContext.from_batch(
            species=batch_or_context.species,
            rows=batch_or_context.rows,
            balances=batch_or_context.balances,
            reactions=batch_or_context.reactions,
            options=batch_or_context.options,
            objective=objective,
            vapor_species=batch_or_context.vapor_species,
            volatile_species=batch_or_context.volatile_species,
            nonvolatile_species=batch_or_context.nonvolatile_species,
            base_parameters=batch_or_context.base_parameters,
            user_options=batch_or_context.user_options,
            mixture_factory=batch_or_context.mixture_factory,
            mixture_factory_builder=batch_or_context.mixture_factory_builder,
            reactive_speciation_options=batch_or_context.reactive_speciation_options,
            reactive_bubble_options=batch_or_context.reactive_bubble_options,
        )
    return context.evaluate_objective(parameter_map)


def summarize_regression_result(result: ReactiveRegressionObjectiveResult) -> dict[str, Any]:
    objective_result = result
    batch = objective_result.batch_result
    row_norms: list[tuple[str, float]] = []
    source_counter: dict[str, list[float]] = {}
    split_counter: dict[str, list[float]] = {}
    target_counter: dict[str, list[float]] = {}
    species_counter: dict[str, list[float]] = {}
    for row in batch.row_results:
        norm = float(np.linalg.norm(np.asarray(list(row.residuals.values()), dtype=float)) if row.residuals else 0.0)
        row_norms.append((row.row_id, norm))
        source_counter.setdefault(row.source or "unspecified", []).append(norm)
        split_counter.setdefault(row.split or "unspecified", []).append(norm)
        for name, value in row.residuals.items():
            target_counter.setdefault(_residual_family_from_name(name), []).append(float(abs(value)))
            species = name.rsplit(".", 1)[-1]
            species_counter.setdefault(species, []).append(float(abs(value)))
    row_norms.sort(key=lambda item: item[1], reverse=True)
    payload = {
        "objective": objective_result.objective,
        "success_count": batch.success_count,
        "failure_count": batch.failure_count,
        "residual_norm": float(np.linalg.norm(objective_result.residuals)),
        "top_failing_rows": [{"row_id": row_id, "residual_norm": norm} for row_id, norm in row_norms[:5]],
        "by_source": {
            source: {
                "count": len(values),
                "mean_row_residual_norm": float(np.mean(values)) if values else 0.0,
            }
            for source, values in source_counter.items()
        },
        "by_split": {
            split: {
                "count": len(values),
                "mean_row_residual_norm": float(np.mean(values)) if values else 0.0,
            }
            for split, values in split_counter.items()
        },
        "by_target_type": {
            target: {
                "count": len(values),
                "mean_abs_residual": float(np.mean(values)) if values else 0.0,
            }
            for target, values in target_counter.items()
        },
        "by_species": {
            species: {
                "count": len(values),
                "mean_abs_residual": float(np.mean(values)) if values else 0.0,
            }
            for species, values in species_counter.items()
        },
        "train_validation": {
            split: {
                "count": len(values),
                "residual_norm": float(np.linalg.norm(np.asarray(values, dtype=float))) if values else 0.0,
            }
            for split, values in split_counter.items()
        },
        "target_family_summaries": _json_like(
            batch.diagnostics.get("target_family_summaries", {}) if isinstance(batch.diagnostics, Mapping) else {}
        ),
        "residual_block_norms": {
            family: float(summary.get("residual_block_norm", 0.0))
            for family, summary in (
                batch.diagnostics.get("target_family_summaries", {}) if isinstance(batch.diagnostics, Mapping) else {}
            ).items()
        },
        "cache_stats": _json_like(batch.cache_stats),
        "timing_summary": _json_like(batch.timing_summary),
    }
    return payload


def write_regression_summary(result: ReactiveRegressionObjectiveResult, path: str | Path) -> Path:
    target = Path(path)
    payload = summarize_regression_result(result)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.suffix.lower() == ".json":
        target.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
        return target
    lines = [f"{key}: {value}" for key, value in payload.items()]
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target


def write_regression_row_table(result: ReactiveRegressionObjectiveResult, path: str | Path) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    rows = result.batch_result.flatten_rows()
    _write_csv_rows(target, rows)
    return target


def write_regression_parameter_table(
    parameter_map: Mapping[str, float],
    path: str | Path,
    *,
    seed_map: Mapping[str, float] | None = None,
) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    parameter_map = {str(k): float(v) for k, v in parameter_map.items()}
    seed_map = {str(k): float(v) for k, v in (seed_map or {}).items()}
    for key, value in sorted((str(k), float(v)) for k, v in parameter_map.items()):
        seed = seed_map.get(key)
        rows.append(
            {
                "parameter": key,
                "seed_value": seed,
                "final_value": value,
                "movement": None if seed is None else float(value - seed),
                "relative_movement": None if seed in (None, 0.0) else float((value - seed) / seed),
            }
        )
    _write_csv_rows(target, rows)
    return target


def write_regression_residual_table(result: ReactiveRegressionObjectiveResult, path: str | Path) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    _write_csv_rows(target, result.batch_result.flatten_residuals())
    return target


def _write_csv_rows(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    fieldnames: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            key_str = str(key)
            if key_str not in seen:
                seen.add(key_str)
                fieldnames.append(key_str)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _csv_scalar(value) for key, value in row.items()})


def _csv_scalar(value: Any) -> Any:
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(_json_like(value), sort_keys=True)
    return value


def _reaction_names(reactions: Sequence[ReactionDefinition]) -> tuple[str, ...]:
    names = []
    for idx, reaction in enumerate(reactions):
        names.append(reaction.name or f"reaction_{idx}")
    return tuple(names)


def _compile_residual_name_template(
    rows: Sequence[ReactiveElectrolyteRow],
    reaction_names: Sequence[str],
) -> tuple[str, ...]:
    names: list[str] = []
    for row in rows:
        for key in row.target_partial_pressures:
            names.append(f"{row.row_id}.partial_pressure.{key}")
        if row.target_pressure is not None:
            names.append(f"{row.row_id}.pressure")
        for key in row.target_speciation:
            names.append(f"{row.row_id}.x.{key}")
        for key in row.target_activity:
            names.append(f"{row.row_id}.activity.{key}")
        for key in row.target_fugacity:
            names.append(f"{row.row_id}.ln_phi.{key}")
        if row.target_density is not None:
            names.append(f"{row.row_id}.density")
        if row.target_relative_permittivity is not None:
            names.append(f"{row.row_id}.relative_permittivity")
        for name in reaction_names:
            names.append(f"{row.row_id}.reaction.{name}")
    return tuple(names)


def _objective_metrics(batch_result: ReactiveElectrolyteBatchResult) -> dict[str, Any]:
    rows = batch_result.flatten_rows()
    return {
        "success_count": batch_result.success_count,
        "failure_count": batch_result.failure_count,
        "residual_count": int(batch_result.residuals.size),
        "residual_norm": float(np.linalg.norm(batch_result.residuals)),
        "mean_row_residual_norm": float(np.mean([row["residual_norm"] for row in rows])) if rows else 0.0,
    }


def _build_row_mixture(
    batch: ReactiveElectrolyteBatch,
    row: ReactiveElectrolyteRow,
    parameter_map: Mapping[str, float],
    *,
    x_override: Sequence[float] | None,
    P_override: float,
) -> ePCSAFTMixture:
    x_array = np.asarray(x_override if x_override is not None else _row_initial_x(row, batch.species), dtype=float)
    if batch.mixture_factory_builder is not None:
        mixture_factory = batch.mixture_factory_builder(parameter_map)
        return mixture_factory(x_array, float(row.T), float(P_override))
    if batch.base_parameters is not None:
        params = _copy_parameter_payload(batch.base_parameters)
        _apply_parameter_map(params, batch.species, parameter_map)
        return ePCSAFTMixture.from_params(params, species=batch.species)
    assert batch.mixture_factory is not None
    return batch.mixture_factory(x_array, float(row.T), float(P_override))


def _mixture_factory_from_batch(batch: ReactiveElectrolyteBatch, parameter_map: Mapping[str, float]) -> Any:
    def factory(x: Sequence[float], _T: float, P: float) -> ePCSAFTMixture:
        return _build_row_mixture(batch, batch.rows[0], parameter_map, x_override=x, P_override=P)

    return factory


def _row_initial_x(
    row: ReactiveElectrolyteRow,
    species: Sequence[str],
) -> np.ndarray:
    if isinstance(row.initial_x, Mapping):
        return np.asarray([float(row.initial_x.get(label, 0.0)) for label in species], dtype=float)
    if row.initial_x is None:
        if row.totals:
            total = sum(float(v) for v in row.totals.values())
            if total > 0.0:
                values = [float(row.totals.get(label, 0.0)) / total for label in species]
                return np.asarray(values, dtype=float)
        return np.full(len(species), 1.0 / len(species), dtype=float)
    return np.asarray(row.initial_x, dtype=float).reshape(-1)


def _row_result_from_speciation(
    *,
    context: ReactiveElectrolyteRegressionContext,
    row: ReactiveElectrolyteRow,
    raw_result: ReactiveSpeciationResult,
    elapsed_seconds: float,
) -> ReactiveElectrolyteRowResult:
    composition = dict(raw_result.x)
    density = None
    relative_permittivity = None
    ln_fugacity: dict[str, float] = {}
    if context.batch.options.include_state_outputs and raw_result.success:
        mixture = _build_row_mixture(
            context.batch,
            row,
            {},
            x_override=[composition[label] for label in context.species],
            P_override=float(row.P),
        )
        state = mixture.state(T=float(row.T), P=float(row.P), x=list(composition.values()))
        density = float(state.molar_density())
        relative_permittivity = float(state.relative_permittivity()[0])
        ln_fugacity = {label: float(value) for label, value in zip(context.species, state.fugacity_coefficient())}
    return ReactiveElectrolyteRowResult(
        row_id=row.row_id,
        success=raw_result.success,
        message=raw_result.message,
        composition=composition,
        pressure=float(row.P),
        ln_fugacity=ln_fugacity,
        activity_coefficients=dict(raw_result.activity_coefficients),
        density=density,
        relative_permittivity=relative_permittivity,
        residuals={},
        residual_names=(),
        failure_diagnostics=dict(raw_result.diagnostics if not raw_result.success else {}),
        active_bounds={},
        elapsed_seconds=float(elapsed_seconds),
        named_reaction_residuals=dict(raw_result.named_reaction_residuals),
        source=row.source,
        split=row.split,
        metadata=dict(row.metadata),
    )


def _row_result_from_bubble(
    *,
    context: ReactiveElectrolyteRegressionContext,
    row: ReactiveElectrolyteRow,
    raw_result: ReactiveElectrolyteBubbleResult,
    parameter_map: Mapping[str, float],
    elapsed_seconds: float,
) -> ReactiveElectrolyteRowResult:
    composition = dict(getattr(raw_result, "x_liq", {}) or {})
    density = None
    relative_permittivity = None
    ln_fugacity: dict[str, float] = {}
    if context.batch.options.include_state_outputs and bool(getattr(raw_result, "success", False)):
        mixture = _build_row_mixture(
            context.batch,
            row,
            parameter_map,
            x_override=[composition[label] for label in context.species],
            P_override=float(getattr(raw_result, "P_total", row.P)),
        )
        state = mixture.state(
            T=float(row.T),
            P=float(getattr(raw_result, "P_total", row.P)),
            x=[composition[label] for label in context.species],
            phase="liq",
        )
        density = float(state.molar_density())
        relative_permittivity = float(state.relative_permittivity()[0])
        ln_fugacity = {label: float(value) for label, value in zip(context.species, state.fugacity_coefficient())}
    return ReactiveElectrolyteRowResult(
        row_id=row.row_id,
        success=bool(getattr(raw_result, "success", False)),
        message=str(getattr(raw_result, "message", "")),
        composition=composition,
        pressure=float(getattr(raw_result, "P_total", row.P)),
        ln_fugacity=ln_fugacity,
        activity_coefficients=dict(getattr(raw_result, "activity_coefficients", {}) or {}),
        density=density,
        relative_permittivity=relative_permittivity,
        residuals={},
        residual_names=(),
        failure_diagnostics=(
            dict(getattr(raw_result, "diagnostics", {}) or {})
            if not bool(getattr(raw_result, "success", False))
            else {}
        ),
        active_bounds={},
        elapsed_seconds=float(elapsed_seconds),
        partial_pressures=dict(getattr(raw_result, "partial_pressures", {}) or {}),
        y_vap=dict(getattr(raw_result, "y_vap", {}) or {}),
        named_reaction_residuals=dict(getattr(raw_result, "named_reaction_residuals", {}) or {}),
        source=row.source,
        split=row.split,
        metadata=dict(row.metadata),
    )


def _failed_row_result(
    *,
    row: ReactiveElectrolyteRow,
    exc: Exception,
    elapsed_seconds: float,
) -> ReactiveElectrolyteRowResult:
    diagnostics = dict(getattr(exc, "diagnostics", {}) or {})
    return ReactiveElectrolyteRowResult(
        row_id=row.row_id,
        success=False,
        message=f"{type(exc).__name__}: {str(exc).splitlines()[0]}",
        composition={},
        pressure=None,
        ln_fugacity={},
        activity_coefficients={},
        density=None,
        relative_permittivity=None,
        residuals={},
        residual_names=(),
        failure_diagnostics=diagnostics,
        active_bounds={},
        elapsed_seconds=float(elapsed_seconds),
        source=row.source,
        split=row.split,
        metadata=dict(row.metadata),
    )


def _pack_row_residuals(
    *,
    row: ReactiveElectrolyteRow,
    row_result: ReactiveElectrolyteRowResult,
    objective: ReactiveRegressionObjective,
    reaction_names: Sequence[str],
    penalty: float,
) -> tuple[tuple[str, ...], np.ndarray]:
    names: list[str] = []
    values: list[float] = []
    row_weight = float(objective.row_weights.get(row.row_id, 1.0))
    source_weight = float(objective.source_weights.get(row.source, 1.0))
    split_weight = float(objective.split_weights.get(row.split, 1.0))
    base_scale = math.sqrt(row_weight * source_weight * split_weight)
    penalty_value = base_scale * float(objective.failure_penalty if not row_result.success else penalty)
    family_scale = {
        "partial_pressure": math.sqrt(float(objective.residual_weights.get("partial_pressure", 1.0))),
        "pressure": math.sqrt(float(objective.residual_weights.get("pressure", 1.0))),
        "speciation": math.sqrt(float(objective.residual_weights.get("speciation", 1.0))),
        "activity": math.sqrt(float(objective.residual_weights.get("activity", 1.0))),
        "fugacity": math.sqrt(float(objective.residual_weights.get("fugacity", 1.0))),
        "density": math.sqrt(float(objective.residual_weights.get("density", 1.0))),
        "relative_permittivity": math.sqrt(float(objective.residual_weights.get("relative_permittivity", 1.0))),
        "reaction": math.sqrt(float(objective.residual_weights.get("reaction", 1.0))),
    }

    for label, target in row.target_partial_pressures.items():
        names.append(f"{row.row_id}.partial_pressure.{label}")
        predicted = float(row_result.partial_pressures.get(label, float("nan")))
        raw = penalty_value if not row_result.success else _normalize_log_ratio(predicted, float(target), penalty_value)
        values.append(base_scale * family_scale["partial_pressure"] * _apply_clip(raw, objective.residual_clip))
    if row.target_pressure is not None:
        names.append(f"{row.row_id}.pressure")
        predicted = float(row_result.pressure if row_result.pressure is not None else float("nan"))
        target = float(row.target_pressure)
        if not row_result.success:
            raw = penalty_value
        elif objective.pressure_family == "pressure_linear":
            raw = penalty_value if not math.isfinite(predicted) else (predicted - target) / max(abs(target), 1.0)
        else:
            raw = _normalize_log_ratio(predicted, target, penalty_value)
        values.append(base_scale * family_scale["pressure"] * _apply_clip(raw, objective.residual_clip))
    for label, target in row.target_speciation.items():
        names.append(f"{row.row_id}.x.{label}")
        predicted = float(row_result.composition.get(label, float("nan")))
        if not row_result.success:
            raw = penalty_value
        elif objective.speciation_family == "speciation_mole_fraction":
            raw = penalty_value if not math.isfinite(predicted) else (predicted - target) / max(abs(target), 1.0e-12)
        else:
            raw = _normalize_log_ratio(predicted, float(target), penalty_value)
        values.append(base_scale * family_scale["speciation"] * _apply_clip(raw, objective.residual_clip))
    for label, target in row.target_activity.items():
        names.append(f"{row.row_id}.activity.{label}")
        predicted = float(row_result.activity_coefficients.get(label, float("nan")))
        if not row_result.success:
            raw = penalty_value
        elif objective.activity_family == "activity_coefficient":
            raw = penalty_value if not math.isfinite(predicted) else (predicted - target) / max(abs(target), 1.0e-12)
        else:
            raw = _normalize_log_ratio(predicted, float(target), penalty_value)
        values.append(base_scale * family_scale["activity"] * _apply_clip(raw, objective.residual_clip))
    for label, target in row.target_fugacity.items():
        names.append(f"{row.row_id}.ln_phi.{label}")
        predicted = float(row_result.ln_fugacity.get(label, float("nan")))
        if not row_result.success:
            raw = penalty_value
        elif objective.fugacity_family == "fugacity_coefficient":
            pred_phi = math.exp(predicted) if math.isfinite(predicted) else float("nan")
            raw = penalty_value if not math.isfinite(pred_phi) else (pred_phi - target) / max(abs(target), 1.0e-12)
        else:
            raw = penalty_value if not math.isfinite(predicted) else predicted - float(target)
        values.append(base_scale * family_scale["fugacity"] * _apply_clip(raw, objective.residual_clip))
    if row.target_density is not None:
        names.append(f"{row.row_id}.density")
        predicted = float(row_result.density if row_result.density is not None else float("nan"))
        raw = (
            penalty_value
            if not row_result.success or not math.isfinite(predicted)
            else (predicted - float(row.target_density)) / max(abs(float(row.target_density)), 1.0)
        )
        values.append(base_scale * family_scale["density"] * _apply_clip(raw, objective.residual_clip))
    if row.target_relative_permittivity is not None:
        names.append(f"{row.row_id}.relative_permittivity")
        predicted = float(
            row_result.relative_permittivity if row_result.relative_permittivity is not None else float("nan")
        )
        raw = (
            penalty_value
            if not row_result.success or not math.isfinite(predicted)
            else (predicted - float(row.target_relative_permittivity))
            / max(abs(float(row.target_relative_permittivity)), 1.0)
        )
        values.append(base_scale * family_scale["relative_permittivity"] * _apply_clip(raw, objective.residual_clip))
    for name in reaction_names:
        names.append(f"{row.row_id}.reaction.{name}")
        predicted = (
            penalty_value
            if not row_result.success
            else float(row_result.named_reaction_residuals.get(name, penalty_value))
        )
        if not math.isfinite(predicted):
            predicted = penalty_value
        values.append(base_scale * family_scale["reaction"] * _apply_clip(predicted, objective.residual_clip))
    residual_map = {name: float(value) for name, value in zip(names, values)}
    object.__setattr__(row_result, "residuals", residual_map)
    object.__setattr__(row_result, "residual_names", tuple(names))
    return tuple(names), np.asarray(values, dtype=float)


def _apply_clip(value: float, clip: float | None) -> float:
    if clip is None:
        return float(value)
    return float(max(-clip, min(clip, value)))


def _residual_family_from_name(name: str) -> str:
    if ".partial_pressure." in name:
        return "partial_pressure"
    if name.endswith(".pressure"):
        return "pressure"
    if ".x." in name:
        return "speciation"
    if ".activity." in name:
        return "activity"
    if ".ln_phi." in name:
        return "fugacity"
    if name.endswith(".density"):
        return "density"
    if name.endswith(".relative_permittivity"):
        return "relative_permittivity"
    if ".reaction." in name:
        return "reaction"
    return "unknown"


def _copy_parameter_payload(payload: Mapping[str, Any]) -> dict[str, Any]:
    copied: dict[str, Any] = {}
    for key, value in payload.items():
        if isinstance(value, np.ndarray):
            copied[str(key)] = np.asarray(value).copy()
        elif isinstance(value, list):
            copied[str(key)] = list(value)
        elif isinstance(value, tuple):
            copied[str(key)] = list(value)
        elif isinstance(value, Mapping):
            copied[str(key)] = _copy_parameter_payload(value)
        else:
            copied[str(key)] = value
    return copied


def _apply_parameter_map(params: dict[str, Any], species: Sequence[str], parameter_map: Mapping[str, float]) -> None:
    if not parameter_map:
        return
    species_index = {label: idx for idx, label in enumerate(species)}
    for label, value in parameter_map.items():
        _apply_parameter_value(params, species_index, label, float(value))


def _apply_parameter_value(
    params: dict[str, Any],
    species_index: Mapping[str, int],
    parameter_label: str,
    value: float,
) -> None:
    if ":" in parameter_label and "." in parameter_label:
        pair, raw_field = parameter_label.split(".", 1)
        left, right = pair.split(":", 1)
        field = _BINARY_FIELD_ALIASES.get(raw_field.strip().lower())
        if field is None:
            raise InputError(f"Unsupported binary parameter label: {parameter_label}")
        i = species_index[str(left)]
        j = species_index[str(right)]
        matrix = np.asarray(
            params.get(field, np.zeros((len(species_index), len(species_index)), dtype=float)), dtype=float
        )
        if matrix.shape != (len(species_index), len(species_index)):
            matrix = np.zeros((len(species_index), len(species_index)), dtype=float)
        matrix[i, j] = float(value)
        matrix[j, i] = float(value)
        params[field] = matrix
        return
    if "." not in parameter_label:
        raise InputError(f"Unsupported parameter label: {parameter_label}")
    component, raw_field = parameter_label.split(".", 1)
    field = _PARAMETER_FIELD_ALIASES.get(raw_field.strip().lower())
    if field is None:
        raise InputError(f"Unsupported pure parameter label: {parameter_label}")
    idx = species_index[str(component)]
    arr = np.asarray(params.get(field), dtype=object if field == "assoc_scheme" else float)
    if arr.size != len(species_index):
        raise InputError(f"Parameter field {field!r} is not aligned with the species list.")
    arr = arr.copy()
    arr[idx] = value
    params[field] = arr


__all__ = [
    "ReactiveElectrolyteBatch",
    "ReactiveElectrolyteBatchOptions",
    "ReactiveElectrolyteBatchResult",
    "ReactiveElectrolyteRegressionContext",
    "ReactiveElectrolyteRow",
    "ReactiveElectrolyteRowResult",
    "ReactiveRegressionJacobianResult",
    "ReactiveRegressionObjective",
    "ReactiveRegressionObjectiveResult",
    "build_reactive_regression_objective",
    "evaluate_reactive_regression_objective",
    "summarize_regression_result",
    "write_regression_parameter_table",
    "write_regression_residual_table",
    "write_regression_row_table",
    "write_regression_summary",
]
