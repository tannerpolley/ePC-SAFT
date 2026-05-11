"""Canonical parameter records with legacy native-payload conversion."""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from ._types import InputError


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
            molar_mass=float(payload.get("MW", payload.get("molar_mass", 0.0))),
            m=float(payload.get("m", 1.0)),
            sigma=float(payload.get("s", payload.get("sigma", 0.0))),
            epsilon_k=float(payload.get("e", payload.get("epsilon_k", 0.0))),
            charge=float(payload.get("z", payload.get("charge", 0.0))),
            epsilon_k_ab=float(payload.get("e_assoc", payload.get("epsilon_k_ab", 0.0))),
            kappa_ab=float(payload.get("vol_a", payload.get("kappa_ab", 0.0))),
            association_scheme=payload.get("assoc_scheme", payload.get("association_scheme")),
            relative_permittivity=float(payload.get("dielc", payload.get("relative_permittivity", 1.0))),
            born_diameter=float(payload.get("d_born", payload.get("born_diameter", 0.0))),
            solvation_factor=float(payload.get("f_solv", payload.get("solvation_factor", 1.0))),
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

    def __post_init__(self) -> None:
        components = tuple(str(component) for component in self.components)
        object.__setattr__(self, "components", components)
        object.__setattr__(self, "pure_records", tuple(self.pure_records))
        object.__setattr__(self, "binary_records", tuple(self.binary_records))
        object.__setattr__(self, "metadata", dict(self.metadata))
        self.validate()

    @classmethod
    def from_records(
        cls,
        pure_records: Sequence[PureRecord],
        binary_records: Sequence[BinaryRecord] | None = None,
        *,
        metadata: Mapping[str, Any] | None = None,
    ) -> ParameterSet:
        pure = tuple(pure_records)
        return cls(
            components=tuple(str(record.component) for record in pure),
            pure_records=pure,
            binary_records=tuple(binary_records or ()),
            metadata=dict(metadata or {}),
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
        return cls.from_records(pure_records, binary_records, metadata=metadata)

    @classmethod
    def from_dataset(
        cls,
        dataset_name: str,
        species: Sequence[str],
        x: Sequence[float] | None = None,
        T: float = 298.15,
        user_options: Mapping[str, Any] | None = None,
    ) -> ParameterSet:
        from .parameters import get_prop_dict

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
        if errors:
            raise InputError("; ".join(errors))
        return {"valid": True, "component_count": len(self.components), "binary_count": len(self.binary_records)}

    def to_legacy_dict(self) -> dict[str, Any]:
        records = {str(record.component): record for record in self.pure_records}
        ordered = [records[label] for label in self.components]
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
            "z": np.asarray([record.charge for record in ordered], dtype=float),
            "dielc": np.asarray([record.relative_permittivity for record in ordered], dtype=float),
            "d_born": np.asarray([record.born_diameter for record in ordered], dtype=float),
            "f_solv": np.asarray([record.solvation_factor for record in ordered], dtype=float),
        }
        ncomp = len(self.components)
        matrices = {
            "k_ij": np.zeros((ncomp, ncomp), dtype=float),
            "l_ij": np.zeros((ncomp, ncomp), dtype=float),
            "k_hb_ij": np.zeros((ncomp, ncomp), dtype=float),
        }
        index = {label: idx for idx, label in enumerate(self.components)}
        for record in self.binary_records:
            i, j = index[record.components[0]], index[record.components[1]]
            matrices["k_ij"][i, j] = matrices["k_ij"][j, i] = record.k_ij
            matrices["l_ij"][i, j] = matrices["l_ij"][j, i] = record.l_ij
            matrices["k_hb_ij"][i, j] = matrices["k_hb_ij"][j, i] = record.k_hb_ij
        payload.update(matrices)
        return payload

    def to_json(self, path: str | Path | None = None) -> str:
        payload = {
            "components": list(self.components),
            "pure_records": [asdict(record) for record in self.pure_records],
            "binary_records": [asdict(record) for record in self.binary_records],
            "metadata": dict(self.metadata),
        }
        text = json.dumps(payload, indent=2, sort_keys=True)
        if path is not None:
            Path(path).write_text(text + "\n", encoding="utf-8")
        return text


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
