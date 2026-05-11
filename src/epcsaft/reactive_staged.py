"""Staged reactive-equilibrium workflow helpers."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from ._types import InputError
from .reactive_speciation import (
    ReactiveSpeciationOptions,
    ReactiveSpeciationResult,
    _json_like,
    solve_reactive_speciation,
)


@dataclass(frozen=True, slots=True)
class ReactiveStagedEquilibriumResult:
    """Chemical-equilibrium result composed with an existing phase route."""

    success: bool
    message: str
    z: Mapping[str, float]
    chemical: ReactiveSpeciationResult
    phase: Any
    diagnostics: Mapping[str, Any]

    def __post_init__(self) -> None:
        object.__setattr__(self, "success", bool(self.success))
        object.__setattr__(self, "message", str(self.message))
        object.__setattr__(self, "z", {str(k): float(v) for k, v in self.z.items()})
        object.__setattr__(self, "diagnostics", dict(self.diagnostics))

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-like staged workflow payload."""
        phase_payload = self.phase.to_dict() if hasattr(self.phase, "to_dict") else self.phase
        return {
            "success": self.success,
            "message": self.message,
            "z": dict(self.z),
            "chemical": self.chemical.to_dict(),
            "phase": _json_like(phase_payload),
            "diagnostics": _json_like(self.diagnostics),
        }


def solve_reactive_staged_equilibrium(
    *,
    species: Sequence[str],
    mixture_factory: Any,
    T: float,
    P: float,
    balances: Mapping[str, Mapping[str, float]],
    totals: Mapping[str, float],
    reactions: Any,
    initial_x: Any,
    phase_kind: str = "auto",
    speciation_options: ReactiveSpeciationOptions | None = None,
    phase_options: Any = None,
    phase_kwargs: Mapping[str, Any] | None = None,
) -> ReactiveStagedEquilibriumResult:
    """Solve chemical equilibrium first, then pass the staged feed to a phase route.

    This helper is intentionally staged and explicit; it does not claim to be a
    fully coupled reactive flash calculation.
    """
    labels = [str(label) for label in species]
    kind = str(phase_kind).strip()
    if not kind:
        raise InputError("phase_kind must be a non-empty equilibrium route label.")
    if kind == "reactive_flash_tp":
        raise InputError("reactive_flash_tp is not exposed; use an explicit staged phase_kind.")
    extra_phase_kwargs = dict(phase_kwargs or {})
    if "z" in extra_phase_kwargs:
        raise InputError("phase_kwargs must not include z; the staged chemical composition is used as the feed.")
    if "kind" in extra_phase_kwargs:
        raise InputError("phase_kwargs must not include kind; use phase_kind.")

    chemical = solve_reactive_speciation(
        species=labels,
        mixture_factory=mixture_factory,
        T=T,
        P=P,
        balances=balances,
        totals=totals,
        reactions=reactions,
        initial_x=initial_x,
        options=speciation_options,
    )
    z = {label: chemical.x[label] for label in labels}
    z_vector = [z[label] for label in labels]
    mixture = mixture_factory(z_vector, T, P)
    phase = _solve_phase_route(
        mixture,
        kind=kind,
        T=T,
        P=P,
        z=z_vector,
        phase_options=phase_options,
        phase_kwargs=extra_phase_kwargs,
    )
    phase_success = bool(getattr(phase, "success", True))
    success = bool(chemical.success and phase_success)
    diagnostics = {
        "workflow": "chemical_equilibrium_then_phase_equilibrium",
        "reactive_workflow_class": "staged",
        "phase_kind": kind,
        "chemical_success": bool(chemical.success),
        "phase_success": phase_success,
        "phase_problem_kind": str(getattr(phase, "problem_kind", kind)),
        "staged_feed": dict(z),
    }
    return ReactiveStagedEquilibriumResult(
        success=success,
        message="converged" if success else "staged reactive equilibrium did not converge",
        z=z,
        chemical=chemical,
        phase=phase,
        diagnostics=diagnostics,
    )


def _solve_phase_route(
    mixture: Any,
    *,
    kind: str,
    T: float,
    P: float,
    z: Sequence[float],
    phase_options: Any,
    phase_kwargs: Mapping[str, Any],
) -> Any:
    route = kind.strip()
    aliases = {
        "flash_tp": "tp_flash",
        "lle_tp": "lle_flash",
        "electrolyte_lle_tp": "electrolyte_lle",
        "stability_tp": "stability",
        "electrolyte_stability_tp": "electrolyte_stability",
        "electrolyte_bubble_p": "electrolyte_bubble_pressure",
        "electrolyte_bubble": "electrolyte_bubble_pressure",
    }
    route = aliases.get(route, route)
    if route == "auto":
        if not hasattr(mixture, "equilibrium"):
            raise InputError("mixture_factory must return an object with an equilibrium method for auto phase routes.")
        return mixture.equilibrium(kind="auto", T=T, P=P, z=z, options=phase_options, **phase_kwargs)
    if route == "tp_flash":
        return mixture.flash_tp(T=T, P=P, z=z, options=phase_options)
    if route == "lle_flash":
        return mixture.lle_tp(
            T=T,
            P=P,
            z=z,
            options=phase_options,
            initial_phases=phase_kwargs.get("initial_phases"),
        )
    if route == "stability":
        return mixture.stability_tp(
            T=T,
            P=P,
            z=z,
            options=phase_options,
            parent_phase=phase_kwargs.get("parent_phase"),
            trial_phases=phase_kwargs.get("trial_phases"),
        )
    if route == "electrolyte_lle":
        return mixture.electrolyte_lle_tp(
            T=T,
            P=P,
            z=z,
            solvent_feed=phase_kwargs.get("solvent_feed"),
            salt_molality=phase_kwargs.get("salt_molality"),
            initial_phases=phase_kwargs.get("initial_phases"),
            options=phase_options,
        )
    if route == "electrolyte_stability":
        return mixture.electrolyte_stability_tp(
            T=T,
            P=P,
            z=z,
            solvent_feed=phase_kwargs.get("solvent_feed"),
            salt_molality=phase_kwargs.get("salt_molality"),
            options=phase_options,
        )
    if route == "electrolyte_bubble_pressure":
        bubble_options = phase_options
        if bubble_options is None:
            from .electrolyte_bubble import ElectrolyteBubbleOptions

            bubble_options = ElectrolyteBubbleOptions(initial_pressure=float(P))
        return mixture.electrolyte_bubble_p(
            T=T,
            z=z,
            vapor_species=phase_kwargs.get("vapor_species"),
            volatile_species=phase_kwargs.get("volatile_species"),
            nonvolatile_species=phase_kwargs.get("nonvolatile_species"),
            options=bubble_options,
        )
    raise InputError(
        "phase_kind must be one of auto, tp_flash, lle_flash, stability, electrolyte_lle, "
        "electrolyte_stability, or electrolyte_bubble_pressure."
    )
