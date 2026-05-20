"""Canonical parameter records with legacy native-payload conversion."""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from ._types import InputError

_PURE_PARAMETER_KEYS = {
    "MW",
    "molar_mass",
    "m",
    "s",
    "sigma",
    "e",
    "epsilon_k",
    "z",
    "charge",
    "e_assoc",
    "epsilon_k_ab",
    "vol_a",
    "kappa_ab",
    "assoc_scheme",
    "association_scheme",
    "dielc",
    "relative_permittivity",
    "d_born",
    "born_diameter",
    "f_solv",
    "solvation_factor",
}
_BINARY_PARAMETER_KEYS = {"k_ij", "l_ij", "k_hb", "k_hb_ij"}
_GENERATED_PARAMETER_KEYS = {"assoc_num", "assoc_matrix"}
_STRUCTURAL_KEYS = {
    "schema",
    "schema_version",
    "components",
    "pure_records",
    "binary_records",
    "metadata",
    "runtime_options",
}
_PARAMETER_PAYLOAD_KEYS = _PURE_PARAMETER_KEYS | _BINARY_PARAMETER_KEYS | _GENERATED_PARAMETER_KEYS | _STRUCTURAL_KEYS


@dataclass(frozen=True, slots=True)
class ComponentIdentifier:
    """Stable component identifier for canonical parameter records."""

    name: str
    aliases: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", str(self.name))
        object.__setattr__(self, "aliases", tuple(str(alias) for alias in self.aliases))


@dataclass(frozen=True, slots=True)
class AssociationSite:
    """One named association site on a component."""

    label: str
    kind: str = "generic"

    def __post_init__(self) -> None:
        object.__setattr__(self, "label", str(self.label))
        object.__setattr__(self, "kind", str(self.kind))


@dataclass(frozen=True, slots=True)
class PureRecord:
    """Canonical pure-component ePC-SAFT parameters.

    ``molar_mass`` is stored in kg/mol to match the legacy native ``MW`` payload.
    Use :meth:`from_g_per_mol` when source tables report g/mol.
    """

    component: str | ComponentIdentifier
    molar_mass: float
    m: float
    sigma: float
    epsilon_k: float
    charge: float = 0.0
    epsilon_k_ab: float = 0.0
    kappa_ab: float = 0.0
    association_scheme: str | None = None
    association_sites: tuple[AssociationSite, ...] = ()
    relative_permittivity: float = 1.0
    born_diameter: float = 0.0
    solvation_factor: float = 1.0

    def __post_init__(self) -> None:
        component = self.component.name if isinstance(self.component, ComponentIdentifier) else str(self.component)
        object.__setattr__(self, "component", component)
        object.__setattr__(self, "association_sites", tuple(self.association_sites))
        for field_name in (
            "molar_mass",
            "m",
            "sigma",
            "epsilon_k",
            "charge",
            "epsilon_k_ab",
            "kappa_ab",
            "relative_permittivity",
            "born_diameter",
            "solvation_factor",
        ):
            object.__setattr__(self, field_name, float(getattr(self, field_name)))
        if not np.isfinite(float(self.molar_mass)) or float(self.molar_mass) <= 0.0:
            raise InputError(f"{component}.molar_mass must be finite and positive in kg/mol.")
        if float(self.molar_mass) > 1.0:
            raise InputError(
                f"PureRecord.molar_mass is interpreted as kg/mol. Got {float(self.molar_mass)} for {component}. "
                "Use a kg/mol value such as 18.01528e-3, or use PureRecord.from_g_per_mol(...)."
            )

    @classmethod
    def from_g_per_mol(
        cls,
        component: str | ComponentIdentifier,
        *,
        molar_mass_g_per_mol: float,
        m: float,
        sigma: float,
        epsilon_k: float,
        charge: float = 0.0,
        epsilon_k_ab: float = 0.0,
        kappa_ab: float = 0.0,
        association_scheme: str | None = None,
        association_sites: Sequence[AssociationSite] = (),
        relative_permittivity: float = 1.0,
        born_diameter: float = 0.0,
        solvation_factor: float = 1.0,
    ) -> PureRecord:
        """Construct a pure record from a source molar mass reported in g/mol."""
        value = float(molar_mass_g_per_mol)
        if not np.isfinite(value) or value <= 0.0:
            raise InputError("molar_mass_g_per_mol must be finite and positive.")
        return cls(
            component=component,
            molar_mass=value * 1.0e-3,
            m=m,
            sigma=sigma,
            epsilon_k=epsilon_k,
            charge=charge,
            epsilon_k_ab=epsilon_k_ab,
            kappa_ab=kappa_ab,
            association_scheme=association_scheme,
            association_sites=tuple(association_sites),
            relative_permittivity=relative_permittivity,
            born_diameter=born_diameter,
            solvation_factor=solvation_factor,
        )

    @classmethod
    def from_legacy(cls, component: str, payload: Mapping[str, Any]) -> PureRecord:
        return cls(
            component=component,
            molar_mass=_legacy_float(payload, "MW", "molar_mass", default=0.0),
            m=_legacy_float(payload, "m", default=1.0),
            sigma=_legacy_float(payload, "s", "sigma", default=0.0),
            epsilon_k=_legacy_float(payload, "e", "epsilon_k", default=0.0),
            charge=_legacy_float(payload, "z", "charge", default=0.0),
            epsilon_k_ab=_legacy_float(payload, "e_assoc", "epsilon_k_ab", default=0.0),
            kappa_ab=_legacy_float(payload, "vol_a", "kappa_ab", default=0.0),
            association_scheme=payload.get("assoc_scheme", payload.get("association_scheme")),
            relative_permittivity=_legacy_float(payload, "dielc", "relative_permittivity", default=1.0),
            born_diameter=_legacy_float(payload, "d_born", "born_diameter", default=0.0),
            solvation_factor=_legacy_float(payload, "f_solv", "solvation_factor", default=1.0),
        )


@dataclass(frozen=True, slots=True)
class BinaryRecord:
    """Canonical binary interaction record."""

    components: tuple[str, str]
    k_ij: float = 0.0
    l_ij: float = 0.0
    k_hb_ij: float = 0.0

    def __post_init__(self) -> None:
        if len(self.components) != 2:
            raise InputError("BinaryRecord.components must contain exactly two component labels.")
        object.__setattr__(self, "components", (str(self.components[0]), str(self.components[1])))
        object.__setattr__(self, "k_ij", float(self.k_ij))
        object.__setattr__(self, "l_ij", float(self.l_ij))
        object.__setattr__(self, "k_hb_ij", float(self.k_hb_ij))


@dataclass(frozen=True, slots=True)
class PermittivityRecord:
    """Canonical relative-permittivity record for one component."""

    component: str
    relative_permittivity: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "component", str(self.component))
        object.__setattr__(self, "relative_permittivity", float(self.relative_permittivity))


@dataclass(frozen=True, slots=True)
class ParameterSet:
    """Canonical parameter set that can emit the legacy native payload."""

    components: tuple[str, ...]
    pure_records: tuple[PureRecord, ...]
    binary_records: tuple[BinaryRecord, ...] = ()
    metadata: Mapping[str, Any] = field(default_factory=dict)
    runtime_options: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        components = tuple(str(component) for component in self.components)
        object.__setattr__(self, "components", components)
        object.__setattr__(self, "pure_records", tuple(self.pure_records))
        object.__setattr__(self, "binary_records", tuple(self.binary_records))
        object.__setattr__(self, "metadata", dict(self.metadata))
        object.__setattr__(self, "runtime_options", _copy_payload_mapping(self.runtime_options))
        self.validate()

    @classmethod
    def from_records(
        cls,
        pure_records: Sequence[PureRecord],
        binary_records: Sequence[BinaryRecord] | None = None,
        *,
        metadata: Mapping[str, Any] | None = None,
        runtime_options: Mapping[str, Any] | None = None,
    ) -> ParameterSet:
        pure = tuple(pure_records)
        return cls(
            components=tuple(str(record.component) for record in pure),
            pure_records=pure,
            binary_records=tuple(binary_records or ()),
            metadata=dict(metadata or {}),
            runtime_options=_copy_payload_mapping(runtime_options),
        )

    @classmethod
    def from_dict(
        cls,
        payload: Mapping[str, Any],
        *,
        species: Sequence[str] | None = None,
        metadata: Mapping[str, Any] | None = None,
    ) -> ParameterSet:
        arrays = {str(key): value for key, value in payload.items()}
        if "pure_records" in arrays:
            return cls._from_canonical_payload(arrays, species=species, metadata=metadata)
        if species is None:
            if "components" in arrays:
                species = [str(item) for item in arrays["components"]]
            else:
                ncomp = int(np.asarray(arrays["m"], dtype=float).size)
                species = [str(idx) for idx in range(ncomp)]
        labels = tuple(str(item) for item in species)
        pure_records = []
        for idx, label in enumerate(labels):
            pure_records.append(
                PureRecord.from_legacy(
                    label,
                    {
                        key: _array_value(value, idx, default=None)
                        for key, value in arrays.items()
                        if key
                        in {
                            "MW",
                            "molar_mass",
                            "m",
                            "s",
                            "sigma",
                            "e",
                            "epsilon_k",
                            "z",
                            "charge",
                            "e_assoc",
                            "epsilon_k_ab",
                            "vol_a",
                            "kappa_ab",
                            "assoc_scheme",
                            "association_scheme",
                            "dielc",
                            "relative_permittivity",
                            "d_born",
                            "born_diameter",
                            "f_solv",
                            "solvation_factor",
                        }
                    },
                )
            )
        binary_records = _binary_records_from_legacy(labels, arrays)
        return cls.from_records(
            pure_records,
            binary_records,
            metadata=metadata,
            runtime_options=_runtime_options_from_payload(arrays),
        )

    @classmethod
    def _from_canonical_payload(
        cls,
        payload: Mapping[str, Any],
        *,
        species: Sequence[str] | None = None,
        metadata: Mapping[str, Any] | None = None,
    ) -> ParameterSet:
        pure_payloads = list(payload.get("pure_records", ()))
        if not pure_payloads:
            raise InputError("canonical parameter payload must include pure_records.")
        payload_components = tuple(str(item) for item in payload.get("components", ()))
        if not payload_components:
            payload_components = tuple(str(item.get("component")) for item in pure_payloads)
        if species is not None and tuple(str(item) for item in species) != payload_components:
            raise InputError("canonical parameter dataset species must match payload components in order.")
        payload_metadata = _copy_payload_mapping(payload.get("metadata", {}))
        payload_metadata.update(_copy_payload_mapping(metadata))
        return cls(
            components=payload_components,
            pure_records=tuple(_pure_record_from_canonical(item) for item in pure_payloads),
            binary_records=tuple(_binary_record_from_canonical(item) for item in payload.get("binary_records", ())),
            metadata=payload_metadata,
            runtime_options=_copy_payload_mapping(payload.get("runtime_options", {})),
        )

    @classmethod
    def from_json(
        cls,
        path: str | Path,
        *,
        species: Sequence[str] | None = None,
        metadata: Mapping[str, Any] | None = None,
    ) -> ParameterSet:
        """Load a canonical parameter-set JSON file or folder."""
        json_path = Path(path).expanduser()
        if json_path.is_dir():
            json_path = json_path / "parameter_set.json"
        payload = json.loads(json_path.read_text(encoding="utf-8"))
        if not isinstance(payload, Mapping):
            raise InputError(f"Canonical parameter file '{json_path}' must contain a JSON object.")
        return cls.from_dict(payload, species=species, metadata=metadata)

    @classmethod
    def from_dataset(
        cls,
        dataset_name: str | Path,
        species: Sequence[str],
        x: Sequence[float] | None = None,
        T: float = 298.15,
        user_options: Mapping[str, Any] | None = None,
    ) -> ParameterSet:
        from .parameters import get_prop_dict

        dataset_path = Path(dataset_name).expanduser()
        canonical_path = dataset_path / "parameter_set.json"
        if dataset_path.is_dir() and canonical_path.exists():
            params = cls.from_json(
                canonical_path,
                species=species,
                metadata={"dataset": str(dataset_name), "T": float(T)},
            )
            runtime_options = _deep_update_mapping(
                params.runtime_options,
                _load_canonical_user_options(dataset_path),
            )
            runtime_options = _deep_update_mapping(runtime_options, user_options or {})
            return cls(
                components=params.components,
                pure_records=params.pure_records,
                binary_records=params.binary_records,
                metadata=params.metadata,
                runtime_options=runtime_options,
            )

        composition = x if x is not None else [1.0 / len(species)] * len(species)
        payload = get_prop_dict(dataset_name, species, composition, T, user_options=user_options)
        return cls.from_dict(payload, species=species, metadata={"dataset": str(dataset_name), "T": float(T)})

    def validate(self) -> dict[str, Any]:
        errors: list[str] = []
        seen = set()
        for record in self.pure_records:
            label = str(record.component)
            if label in seen:
                errors.append(f"Duplicate component record for {label}.")
            seen.add(label)
            for field_name in ("molar_mass", "m", "sigma", "epsilon_k"):
                value = float(getattr(record, field_name))
                if not np.isfinite(value) or value <= 0.0:
                    errors.append(f"{label}.{field_name} must be finite and positive.")
            if abs(float(record.charge)) > 1.0e-12 and float(record.born_diameter) <= 0.0:
                errors.append(f"{label}.born_diameter must be positive for charged species.")
        missing = [label for label in self.components if label not in seen]
        if missing:
            errors.append(f"Missing pure records for components: {', '.join(missing)}.")
        known = set(self.components)
        for record in self.binary_records:
            for label in record.components:
                if label not in known:
                    errors.append(f"BinaryRecord references unknown component {label}.")
        reserved = sorted(set(self.runtime_options) & _PARAMETER_PAYLOAD_KEYS)
        if reserved:
            errors.append(f"runtime_options cannot override parameter payload keys: {', '.join(reserved)}.")
        if errors:
            raise InputError("; ".join(errors))
        return {
            "valid": True,
            "component_count": len(self.components),
            "binary_count": len(self.binary_records),
            "runtime_option_count": len(self.runtime_options),
        }

    def to_runtime_dict(self, user_options: Mapping[str, Any] | None = None) -> dict[str, Any]:
        records = {str(record.component): record for record in self.pure_records}
        ordered = [records[label] for label in self.components]
        charge_vector = np.asarray([record.charge for record in ordered], dtype=float)
        payload: dict[str, Any] = {
            "MW": np.asarray([record.molar_mass for record in ordered], dtype=float),
            "m": np.asarray([record.m for record in ordered], dtype=float),
            "s": np.asarray([record.sigma for record in ordered], dtype=float),
            "e": np.asarray([record.epsilon_k for record in ordered], dtype=float),
            "e_assoc": np.asarray([record.epsilon_k_ab for record in ordered], dtype=float),
            "vol_a": np.asarray([record.kappa_ab for record in ordered], dtype=float),
            "assoc_scheme": [
                None if record.association_scheme in (None, "") else str(record.association_scheme)
                for record in ordered
            ],
            "z": charge_vector if np.any(np.abs(charge_vector) > 1.0e-12) else np.asarray([], dtype=float),
            "dielc": np.asarray([record.relative_permittivity for record in ordered], dtype=float),
            "d_born": np.asarray([record.born_diameter for record in ordered], dtype=float),
            "f_solv": np.asarray([record.solvation_factor for record in ordered], dtype=float),
        }
        ncomp = len(self.components)
        matrices = {
            "k_ij": np.zeros((ncomp, ncomp), dtype=float),
            "l_ij": np.zeros((ncomp, ncomp), dtype=float),
            "k_hb": np.zeros((ncomp, ncomp), dtype=float),
        }
        index = {label: idx for idx, label in enumerate(self.components)}
        for record in self.binary_records:
            i, j = index[record.components[0]], index[record.components[1]]
            matrices["k_ij"][i, j] = matrices["k_ij"][j, i] = record.k_ij
            matrices["l_ij"][i, j] = matrices["l_ij"][j, i] = record.l_ij
            matrices["k_hb"][i, j] = matrices["k_hb"][j, i] = record.k_hb_ij
        payload.update(matrices)
        runtime_options = _deep_update_mapping(self.runtime_options, user_options or {})
        reserved = sorted(set(runtime_options) & _PARAMETER_PAYLOAD_KEYS)
        if reserved:
            raise InputError(f"runtime_options cannot override parameter payload keys: {', '.join(reserved)}.")
        payload.update(runtime_options)
        return payload

    def to_legacy_dict(self) -> dict[str, Any]:
        return self.to_runtime_dict()

    def to_json(self, path: str | Path | None = None) -> str:
        payload = {
            "components": list(self.components),
            "pure_records": [asdict(record) for record in self.pure_records],
            "binary_records": [asdict(record) for record in self.binary_records],
            "metadata": _json_ready(self.metadata),
            "runtime_options": _json_ready(self.runtime_options),
        }
        text = json.dumps(payload, indent=2, sort_keys=True)
        if path is not None:
            Path(path).write_text(text + "\n", encoding="utf-8")
        return text


@dataclass(frozen=True, slots=True)
class ParameterSource:
    """Canonical resolver for parameter source labels and runtime payloads."""

    source: str | Path | Mapping[str, Any] | ParameterSet
    species: Sequence[str] | None = None

    def __post_init__(self) -> None:
        if self.species is not None:
            object.__setattr__(self, "species", tuple(str(label) for label in self.species))

    @property
    def label(self) -> str:
        """Return a stable human-readable label for reports and diagnostics."""
        if isinstance(self.source, ParameterSet):
            source = self.source.metadata.get("dataset", self.source.metadata.get("source"))
            return str(source) if source not in (None, "") else "ParameterSet"
        if isinstance(self.source, Mapping):
            metadata = self.source.get("metadata")
            if isinstance(metadata, Mapping):
                source = metadata.get("dataset", metadata.get("source"))
                if source not in (None, ""):
                    return str(source)
            source = self.source.get("dataset", self.source.get("source"))
            return str(source) if source not in (None, "") else "parameter-payload"
        return str(self.source)

    def to_parameter_set(
        self,
        species: Sequence[str] | None = None,
        x: Sequence[float] | None = None,
        T: float = 298.15,
        user_options: Mapping[str, Any] | None = None,
    ) -> ParameterSet:
        """Resolve this source into a canonical parameter set."""
        labels = _resolve_parameter_source_species(self.species, species)
        if isinstance(self.source, ParameterSet):
            _require_parameter_source_species(labels, self.source.components, "ParameterSet")
            return _parameter_set_with_runtime_options(self.source, user_options)
        if isinstance(self.source, Mapping):
            params = ParameterSet.from_dict(self.source, species=labels)
            return _parameter_set_with_runtime_options(params, user_options)
        if labels is None:
            raise InputError("ParameterSource species must be provided for dataset sources.")
        composition = x if x is not None else [1.0 / len(labels)] * len(labels)
        return ParameterSet.from_dataset(self.source, labels, composition, T, user_options=user_options)

    def to_runtime_dict(
        self,
        species: Sequence[str] | None = None,
        x: Sequence[float] | None = None,
        T: float = 298.15,
        user_options: Mapping[str, Any] | None = None,
    ) -> dict[str, Any]:
        """Resolve this source into the native runtime payload."""
        return self.to_parameter_set(species=species, x=x, T=T, user_options=user_options).to_runtime_dict()


def _legacy_float(payload: Mapping[str, Any], *keys: str, default: float) -> float:
    for key in keys:
        value = payload.get(key)
        if value not in (None, ""):
            return float(value)
    return float(default)


def _pure_record_from_canonical(payload: Mapping[str, Any]) -> PureRecord:
    component = str(payload.get("component", "")).strip()
    if not component:
        raise InputError("canonical pure record requires component.")
    molar_mass = _required_float(payload, "molar_mass", component)
    units = str(payload.get("molar_mass_units", "kg/mol")).strip().lower()
    if units in {"g/mol", "g per mol", "g mol^-1"}:
        molar_mass *= 1.0e-3
    elif units not in {"kg/mol", "kg per mol", "kg mol^-1"}:
        raise InputError(f"{component}.molar_mass_units must be 'kg/mol' or 'g/mol'.")
    association_sites = tuple(_association_site_from_canonical(item) for item in payload.get("association_sites", ()))
    return PureRecord(
        component=component,
        molar_mass=molar_mass,
        m=_required_float(payload, "m", component),
        sigma=_required_float(payload, "sigma", component),
        epsilon_k=_required_float(payload, "epsilon_k", component),
        charge=_optional_float(payload, "charge", 0.0),
        epsilon_k_ab=_optional_float(payload, "epsilon_k_ab", 0.0),
        kappa_ab=_optional_float(payload, "kappa_ab", 0.0),
        association_scheme=payload.get("association_scheme"),
        association_sites=association_sites,
        relative_permittivity=_optional_float(payload, "relative_permittivity", 1.0),
        born_diameter=_optional_float(payload, "born_diameter", 0.0),
        solvation_factor=_optional_float(payload, "solvation_factor", 1.0),
    )


def _association_site_from_canonical(payload: Any) -> AssociationSite:
    if isinstance(payload, Mapping):
        return AssociationSite(
            label=str(payload.get("label", "")),
            kind=str(payload.get("kind", "generic")),
        )
    return AssociationSite(label=str(payload))


def _binary_record_from_canonical(payload: Mapping[str, Any]) -> BinaryRecord:
    return BinaryRecord(
        components=tuple(str(item) for item in payload.get("components", ())),
        k_ij=_optional_float(payload, "k_ij", 0.0),
        l_ij=_optional_float(payload, "l_ij", 0.0),
        k_hb_ij=_optional_float(payload, "k_hb_ij", 0.0),
    )


def _required_float(payload: Mapping[str, Any], key: str, component: str) -> float:
    value = payload.get(key)
    if value in (None, ""):
        raise InputError(f"{component}.{key} must be filled in canonical parameter data.")
    return float(value)


def _optional_float(payload: Mapping[str, Any], key: str, default: float) -> float:
    value = payload.get(key, default)
    if value in (None, ""):
        return float(default)
    return float(value)


def _resolve_parameter_source_species(
    base: Sequence[str] | None,
    override: Sequence[str] | None,
) -> tuple[str, ...] | None:
    if base is None and override is None:
        return None
    if base is None:
        return tuple(str(label) for label in override or ())
    labels = tuple(str(label) for label in base)
    if override is not None and tuple(str(label) for label in override) != labels:
        raise InputError("ParameterSource species override must match construction species.")
    return labels


def _require_parameter_source_species(
    requested: Sequence[str] | None,
    available: Sequence[str],
    label: str,
) -> None:
    if requested is None:
        return
    if tuple(str(item) for item in requested) != tuple(str(item) for item in available):
        raise InputError(f"{label} species order must match ParameterSource species.")


def _parameter_set_with_runtime_options(
    params: ParameterSet,
    user_options: Mapping[str, Any] | None,
) -> ParameterSet:
    if not user_options:
        return params
    return ParameterSet(
        components=params.components,
        pure_records=params.pure_records,
        binary_records=params.binary_records,
        metadata=params.metadata,
        runtime_options=_deep_update_mapping(params.runtime_options, user_options),
    )


def _load_canonical_user_options(dataset_dir: Path) -> dict[str, Any]:
    path = dataset_dir / "user_options.json"
    if not path.exists():
        return {}
    payload = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(payload, Mapping) and "canonical_user_options" in payload:
        payload = payload.get("canonical_user_options", {})
    if not isinstance(payload, Mapping):
        return {}
    return _copy_payload_mapping(payload)


def _copy_payload_value(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return np.asarray(value).copy()
    if isinstance(value, Mapping):
        return {str(key): _copy_payload_value(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_copy_payload_value(item) for item in value]
    if isinstance(value, tuple):
        return tuple(_copy_payload_value(item) for item in value)
    return value


def _copy_payload_mapping(payload: Mapping[str, Any] | None) -> dict[str, Any]:
    if payload is None:
        return {}
    return {str(key): _copy_payload_value(value) for key, value in payload.items()}


def _deep_update_mapping(base: Mapping[str, Any] | None, updates: Mapping[str, Any] | None) -> dict[str, Any]:
    merged = _copy_payload_mapping(base)
    for key, value in _copy_payload_mapping(updates).items():
        if isinstance(merged.get(key), Mapping) and isinstance(value, Mapping):
            merged[key] = _deep_update_mapping(merged[key], value)
        else:
            merged[key] = value
    return merged


def _runtime_options_from_payload(payload: Mapping[str, Any]) -> dict[str, Any]:
    return {
        str(key): _copy_payload_value(value)
        for key, value in payload.items()
        if str(key) not in _PARAMETER_PAYLOAD_KEYS
    }


def _json_ready(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return _json_ready(value.tolist())
    if isinstance(value, Mapping):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    return value


def _array_value(value: Any, idx: int, *, default: Any) -> Any:
    if value is None:
        return default
    if isinstance(value, (str, bytes)):
        return value
    array = np.asarray(value, dtype=object)
    if array.ndim == 0:
        return array.item()
    if array.size <= idx:
        return default
    item = array.reshape(-1)[idx]
    return item.item() if hasattr(item, "item") else item


def _binary_records_from_legacy(labels: Sequence[str], payload: Mapping[str, Any]) -> tuple[BinaryRecord, ...]:
    out: list[BinaryRecord] = []
    for i, left in enumerate(labels):
        for j in range(i + 1, len(labels)):
            right = labels[j]
            values = {
                "k_ij": _matrix_value(payload.get("k_ij"), i, j),
                "l_ij": _matrix_value(payload.get("l_ij"), i, j),
                "k_hb_ij": _matrix_value(payload.get("k_hb_ij", payload.get("k_hb")), i, j),
            }
            if any(abs(value) > 0.0 for value in values.values()):
                out.append(BinaryRecord((left, right), **values))
    return tuple(out)


def _matrix_value(value: Any, i: int, j: int) -> float:
    if value is None:
        return 0.0
    matrix = np.asarray(value, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] <= i or matrix.shape[1] <= j:
        return 0.0
    return float(matrix[i, j])
