"""Parameter regression helpers for user-owned ePC-SAFT datasets."""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from ._types import InputError
from ._types import phase_to_int
from .epcsaft import _fit_pure_neutral_native_least_squares
from .epcsaft import _fit_pure_neutral_native_debug
from .parameter_templates import _infer_pure_template_name
from .parameters import _deterministic_default
from .parameters import _invalidate_dataset_cache
from .parameters import _load_dataset
from .parameters import _MISSING
from .parameters import _matrix_value
from .parameters import _normalize_component
from .parameters import _resolve_component_field_with_source
from .parameters import _solvent_token_for_component
from .parameters import get_prop_dict


PURE_NEUTRAL_MODE = "pure_neutral"
PURE_ION_MODE = "pure_ion"
BINARY_PAIR_MODE = "binary_pair"

TERM_DENSITY = "density"
TERM_PURE_VLE = "pure_vle_fugacity_balance"
TERM_OSMOTIC = "osmotic_coefficient"
TERM_MIAC = "mean_ionic_activity"
TERM_BINARY_VLE = "binary_vle_fugacity_balance"
TERM_BINARY_LLE = "binary_lle_fugacity_balance"

PURE_DENSITY_KEYS_MOLAR = ("rho",)
PURE_DENSITY_KEYS_MASS = ("rho_kg_m3", "rho_mass_kg_m3", "rho_liq_kg_m3", "rho_sat_liq_kg_m3")

PURE_NUMERIC_FIELDS = (
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
MATRIX_PARAM_NAMES = {
    "k_ij": "k_ij",
    "l_ij": "l_ij",
    "k_hb_ij": "k_hb",
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

SMALL_VALUE_FLOOR = 1.0e-8


def _copy_mapping(mapping: Mapping[str, Any] | None) -> dict[str, Any]:
    return {} if mapping is None else {str(k): v for k, v in mapping.items()}


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
    assoc_scheme: str = ""
    temperature_model: str = "constant"
    terms: tuple[FitTerm, ...] = field(default_factory=tuple)
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
        self.assoc_scheme = str(self.assoc_scheme or "")
        self.temperature_model = str(self.temperature_model)
        self.terms = tuple(self.terms)
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

    def __post_init__(self) -> None:
        self.fitted_values = {str(k): float(v) for k, v in self.fitted_values.items()}
        rendered: dict[str, str | float] = {}
        for key, value in self.rendered_values.items():
            rendered[str(key)] = float(value) if isinstance(value, (int, float, np.integer, np.floating)) else str(value)
        self.rendered_values = rendered
        self.metrics_by_term = {str(k): float(v) for k, v in self.metrics_by_term.items()}
        self.cost = float(self.cost)
        self.residual_norm = float(self.residual_norm)
        self.success = bool(self.success)
        self.status = int(self.status)
        self.message = str(self.message)
        self.nfev = int(self.nfev)
        self.backend = str(self.backend)


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
        name = _normalize_component(str(key)[len(prefix):])
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


def _relative_residual(calc: float, exp: float) -> float:
    denom = max(abs(exp), SMALL_VALUE_FLOOR)
    return (calc - exp) / denom


def _relative_or_absolute_residual(calc: float, exp: float) -> float:
    if abs(exp) <= SMALL_VALUE_FLOOR:
        return calc - exp
    return (calc - exp) / abs(exp)


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


def _apply_component_overrides(params: dict[str, Any], species: Sequence[str], component: str, updates: Mapping[str, float]) -> dict[str, Any]:
    out = {
        key: (np.asarray(value, dtype=float).copy() if isinstance(value, np.ndarray) else value.copy() if isinstance(value, list) else value)
        for key, value in params.items()
    }
    idx = [str(name) for name in species].index(str(component))
    for field, value in updates.items():
        if field not in PURE_NUMERIC_FIELDS:
            continue
        if field not in out:
            continue
        arr = np.asarray(out[field], dtype=float).copy()
        if arr.size == 0:
            if field == "z":
                arr = np.zeros(len(species), dtype=float)
            else:
                raise InputError(f"Cannot override parameter '{field}' because it is missing in the resolved payload.")
        arr[idx] = float(value)
        out[field] = arr
    return out


def _apply_binary_overrides(
    params: dict[str, Any],
    species: Sequence[str],
    pair: tuple[str, str],
    target_values: Mapping[str, float],
) -> dict[str, Any]:
    out = {
        key: (np.asarray(value, dtype=float).copy() if isinstance(value, np.ndarray) else value.copy() if isinstance(value, list) else value)
        for key, value in params.items()
    }
    names = [str(name) for name in species]
    i = names.index(str(pair[0]))
    j = names.index(str(pair[1]))
    for public_name, matrix_key in MATRIX_PARAM_NAMES.items():
        if public_name not in target_values:
            continue
        matrix = np.asarray(out[matrix_key], dtype=float).copy()
        matrix[i, j] = float(target_values[public_name])
        matrix[j, i] = float(target_values[public_name])
        out[matrix_key] = matrix
    return out


def _normalize_vector_map(names: Sequence[str], values: Sequence[float]) -> dict[str, float]:
    return {str(name): float(value) for name, value in zip(names, values)}


def _binary_target_values(vector_map: Mapping[str, float], fit_targets: Sequence[str], temperature_model: str, T: float) -> dict[str, float]:
    values: dict[str, float] = {}
    if temperature_model == "constant":
        for target in fit_targets:
            values[str(target)] = float(vector_map[str(target)])
        return values
    for target in fit_targets:
        slope = float(vector_map[f"{target}_slope"])
        intercept = float(vector_map[f"{target}_intercept"])
        values[str(target)] = slope * float(T) + intercept
    return values


def _render_binary_values(vector_map: Mapping[str, float], fit_targets: Sequence[str], temperature_model: str) -> dict[str, str | float]:
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
        raise InputError("pure_neutral regression requires at least one saturation record with an experimental 'P' value.")
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
    density_records = [record for record in records if _value_from_record(record, "rho", required=False) is not None]
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
    if not density_records:
        raise InputError("pure_ion regression requires density records with a 'rho' value.")
    if not osmotic_records and not miac_records:
        raise InputError("pure_ion regression requires osmotic and/or mean-ionic activity records.")
    _require_record_value(density_records, PURE_ION_MODE, "P")
    _require_record_value(osmotic_records, PURE_ION_MODE, "P")
    _require_record_value(miac_records, PURE_ION_MODE, "P")
    terms = [_term_summary(density_records, TERM_DENSITY, 1.0, len(density_records))]
    if osmotic_records:
        terms.append(_term_summary(osmotic_records, TERM_OSMOTIC, 1.0, len(osmotic_records)))
    if miac_records:
        terms.append(_term_summary(miac_records, TERM_MIAC, 1.0, len(miac_records)))
    return tuple(terms)


def _build_binary_terms(records: Sequence[dict[str, Any]]) -> tuple[FitTerm, ...]:
    vle_records = [record for record in records if _prefixed_species_values(record, "y_")]
    lle_records = [
        record for record in records if _prefixed_species_values(record, "x_alpha_") and _prefixed_species_values(record, "x_beta_")
    ]
    if not vle_records and not lle_records:
        raise InputError("binary_pair regression requires VLE records with y_* columns and/or LLE records with x_alpha_* and x_beta_* columns.")
    _require_record_value(vle_records, BINARY_PAIR_MODE, "P")
    _require_record_value(lle_records, BINARY_PAIR_MODE, "P")
    terms: list[FitTerm] = []
    if vle_records:
        terms.append(_term_summary(vle_records, TERM_BINARY_VLE, 1.0, 2 * len(vle_records)))
    if lle_records:
        terms.append(_term_summary(lle_records, TERM_BINARY_LLE, 1.0, 2 * len(lle_records)))
    return tuple(terms)
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
        raise InputError("pure_neutral regression requires a fixed MW value, either from the dataset or fixed_parameters.")
    if "z" not in fixed_payload or fixed_payload["z"] in (None, ""):
        fixed_payload["z"] = 0.0

    initial = _copy_mapping(initial_guess)
    initial_map = {target: _seed_value(target, initial, fixed_payload) for target in normalized_fit_targets}
    optimization_names = _optimization_parameter_names(PURE_NEUTRAL_MODE, normalized_fit_targets, "constant")
    for name in ("m", "s", "e"):
        if name in initial_map:
            fixed_payload[name] = float(initial_map[name])
    terms = _build_pure_neutral_terms(normalized_records)
    pure_file_hint = f"{source_key}.csv" if source_key is not None else _infer_pure_template_name([normalized_component])
    return _ensure_native_vector_payload(fixed_payload), initial_map, bounds_obj, terms, optimization_names, pure_file_hint


def _family_metrics(raw_residuals_by_term: Mapping[str, Sequence[float]]) -> dict[str, float]:
    metrics: dict[str, float] = {}
    for family, values in raw_residuals_by_term.items():
        arr = np.asarray(values, dtype=float)
        metrics[family] = float(np.sqrt(np.mean(arr**2))) if arr.size else 0.0
    return metrics


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


def _deferred_workflow_error(mode: str) -> NotImplementedError:
    return NotImplementedError(
        f"{mode} regression is deferred. Phase 1 supports only neutral-component fitting of m, s, and e."
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
    fixed_payload, initial_map, bounds_obj, terms, optimization_names, pure_file_hint = _native_pure_neutral_solver_payload(
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
        "density_T": np.asarray([_float_from_record(record, "T", required=True) for record in density_term.records], dtype=float),
        "density_P": np.asarray([_float_from_record(record, "P", required=True) for record in density_term.records], dtype=float),
        "density_rho_exp": np.asarray(
            [_pure_neutral_density_molar(record, mw_value) for record in density_term.records],
            dtype=float,
        ),
        "density_phase": np.asarray(
            [phase_to_int(_value_from_record(record, "phase", required=False) or "liq") for record in density_term.records],
            dtype=int,
        ),
        "density_scale": float(_family_scale(density_term)),
        "vle_T": np.asarray([_float_from_record(record, "T", required=True) for record in pure_vle_term.records], dtype=float),
        "vle_P": np.asarray([_float_from_record(record, "P", required=True) for record in pure_vle_term.records], dtype=float),
        "pure_vle_scale": float(_family_scale(pure_vle_term)),
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
        backend="least_squares_native",
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
) -> FitResult:
    """Fit neutral pure-component m, s, and e parameters."""
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


def fit_pure_ion(
    records: Any,
    component: str,
    *,
    dataset: str | Path,
    solvent: str | None = None,
    fit_targets: Iterable[str] | None = None,
    initial_guess: Mapping[str, float] | None = None,
    bounds: FitBounds | Mapping[str, tuple[float | None, float | None]] | None = None,
    user_options: Mapping[str, Any] | None = None,
    multistart: int = 0,
) -> FitResult:
    """Deferred until a later regression phase."""

    raise _deferred_workflow_error(PURE_ION_MODE)


def fit_binary_pair(
    records: Any,
    pair: Sequence[str],
    *,
    dataset: str | Path,
    fit_targets: Iterable[str] | None = None,
    temperature_model: str = "constant",
    initial_guess: Mapping[str, float] | None = None,
    bounds: FitBounds | Mapping[str, tuple[float | None, float | None]] | None = None,
    user_options: Mapping[str, Any] | None = None,
    multistart: int = 0,
) -> FitResult:
    """Deferred until a later regression phase."""

    raise _deferred_workflow_error(BINARY_PAIR_MODE)


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
        assert problem.pair is not None
        for target in problem.fit_targets:
            path = bi_dir / MATRIX_FILE_NAMES[target]
            if not path.exists():
                raise FileNotFoundError(f"Expected matrix file '{path}' to exist.")
            _update_matrix_file(path, problem.pair, result.rendered_values[target], overwrite=overwrite)
            written.append(path)
        _invalidate_dataset_cache(root)
        return written

    raise InputError(f"Unsupported fit result mode '{problem.mode}'.")
