"""Internal conversion helpers for native equilibrium payloads."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from .._types import SolutionError


def _diagnostics(payload: Mapping[str, Any]) -> dict[str, Any]:
    diagnostics = payload.get("diagnostics", {})
    if isinstance(diagnostics, Mapping):
        return dict(diagnostics)
    return {}


def _phase_payload_to_public(phase: Mapping[str, Any]):
    from ..equilibrium import EquilibriumPhase

    ln_fugacity = phase.get("ln_fugacity_coefficient")
    fugacity = phase.get("fugacity_coefficient")
    return EquilibriumPhase(
        str(phase["label"]),
        composition=np.asarray(phase["composition"], dtype=float),
        density=float(phase["density"]),
        temperature=float(phase["temperature"]),
        pressure=float(phase["pressure"]),
        phase_fraction=float(phase["phase_fraction"]),
        ln_fugacity_coefficient=None if ln_fugacity is None else np.asarray(ln_fugacity, dtype=float),
        fugacity_coefficient=None if fugacity is None else np.asarray(fugacity, dtype=float),
        diagnostics=_diagnostics(phase),
    )


def neutral_two_phase_payload_to_result(payload: Mapping[str, Any]):
    """Convert an accepted native neutral two-phase payload into public dataclasses."""
    from ..equilibrium import EquilibriumResult

    diagnostics = _diagnostics(payload)
    if not bool(payload.get("accepted", False)):
        reason = str(payload.get("rejection_reason", diagnostics.get("rejection_reason", "native_rejected")))
        if "rejection_reason" not in diagnostics:
            diagnostics["rejection_reason"] = reason
        raise SolutionError(f"Native neutral two-phase EOS result was rejected: {reason}", diagnostics)

    phases_raw = payload.get("phases", ())
    if not isinstance(phases_raw, Sequence) or isinstance(phases_raw, (str, bytes)):
        raise SolutionError("Native neutral two-phase EOS result did not contain a phase sequence.", diagnostics)
    phases = tuple(_phase_payload_to_public(phase) for phase in phases_raw)
    if not phases:
        raise SolutionError("Native neutral two-phase EOS result accepted without phase payloads.", diagnostics)

    return EquilibriumResult(
        backend=str(payload.get("backend", "native_equilibrium_nlp")),
        problem_kind=str(payload.get("problem_kind", "neutral_two_phase_eos")),
        phases=phases,
        stable=bool(payload.get("stable", False)),
        split_detected=bool(payload.get("split_detected", True)),
        diagnostics=diagnostics,
    )
