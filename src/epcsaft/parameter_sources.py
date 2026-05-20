"""Shared helpers for parameter source identity, copies, and option loading."""

from __future__ import annotations

import copy
import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import numpy as np


def copy_parameter_value(value: Any) -> Any:
    """Copy parameter payload values without aliasing mutable arrays or containers."""

    if isinstance(value, np.ndarray):
        return np.asarray(value).copy()
    if isinstance(value, Mapping):
        return {str(key): copy_parameter_value(item) for key, item in value.items()}
    if isinstance(value, list):
        return [copy_parameter_value(item) for item in value]
    if isinstance(value, tuple):
        return tuple(copy_parameter_value(item) for item in value)
    return copy.deepcopy(value)


def copy_parameter_mapping(payload: Mapping[str, Any] | None) -> dict[str, Any]:
    if payload is None:
        return {}
    return {str(key): copy_parameter_value(value) for key, value in payload.items()}


def deep_update_parameter_mapping(
    base: Mapping[str, Any] | None,
    updates: Mapping[str, Any] | None,
) -> dict[str, Any]:
    merged = copy_parameter_mapping(base)
    for key, value in copy_parameter_mapping(updates).items():
        if isinstance(merged.get(key), Mapping) and isinstance(value, Mapping):
            merged[key] = deep_update_parameter_mapping(merged[key], value)
        else:
            merged[key] = value
    return merged


def load_canonical_user_options(dataset_dir: str | Path) -> dict[str, Any]:
    """Load and sanitize dataset-owned canonical user options."""

    path = Path(dataset_dir) / "user_options.json"
    if not path.exists():
        return {}
    payload = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(payload, Mapping) and "canonical_user_options" in payload:
        canonical = payload.get("canonical_user_options", {})
    else:
        canonical = payload
    if not isinstance(canonical, Mapping):
        return {}

    elec_model = canonical.get("elec_model", {})
    if isinstance(elec_model, Mapping) and "polar_model" in elec_model:
        raise ValueError(
            f"Dataset '{Path(dataset_dir)}' canonical_user_options still contains removed key "
            "'elec_model.polar_model'."
        )

    cleaned = copy_parameter_mapping(canonical)
    elec = cleaned.get("elec_model")
    if isinstance(elec, Mapping):
        elec.pop("preset", None)
        elec.pop("base", None)
    return cleaned
