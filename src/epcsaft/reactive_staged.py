"""Sequential reactive-equilibrium workflow helpers, not a full reactive flash."""

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
        """Return a JSON-like sequential workflow payload."""
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
    workflow_options: Mapping[str, Any] | None = None,
) -> ReactiveStagedEquilibriumResult:
    """Solve chemical equilibrium first, then pass the speciated feed to a phase route.

    This helper is intentionally sequential and explicit; it does not claim to be a
    fully coupled reactive flash calculation.
    """
    labels = [str(label) for label in species]
    kind = str(phase_kind).strip()
    if not kind:
        raise InputError("phase_kind must be a non-empty equilibrium route label.")
    if kind == "reactive_flash_tp":
        raise InputError("reactive_flash_tp is not exposed; use an explicit downstream phase_kind.")
    extra_phase_kwargs = dict(phase_kwargs or {})
    if "z" in extra_phase_kwargs:
        raise InputError("phase_kwargs must not include z; the chemical-equilibrium composition is used as the feed.")
    if "kind" in extra_phase_kwargs:
        raise InputError("phase_kwargs must not include kind; use phase_kind.")
    workflow = _normalize_workflow_options(workflow_options)

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
    phase_audit = _phase_equilibrium_audit(phase)
    chemical_audit = _chemical_equilibrium_audit(chemical)
    diagnostics = {
        "workflow": "chemical_equilibrium_then_phase_equilibrium",
        "reactive_workflow_class": "staged",
        "reactive_phase_method": "chemical_equilibrium_then_phase_equilibrium",
        "coupling_level": "staged_not_full_simultaneous_nlp",
        "reaction_constant_policy": "fixed_literature_constants_first",
        "reaction_constant_fitting_role": workflow["reaction_constant_fitting"],
        "parameter_regression_boundary": "fit_epcsaft_parameters_after_fixed_constant_speciation",
        "full_simultaneous_reactive_nlp": False,
        "reaction_coordinates": chemical_audit["reaction_coordinates"],
        "element_balance_residuals": chemical_audit["element_balance_residuals"],
        "reaction_equilibrium_residuals": chemical_audit["reaction_equilibrium_residuals"],
        "nonnegativity": chemical_audit["nonnegativity"],
        "phase_split": phase_audit["phase_split"],
        "fugacity_equality": phase_audit["fugacity_equality"],
        "material_balance_error": phase_audit["material_balance_error"],
        "ascani_benchmark_attempt": _ascani_benchmark_attempt(kind),
        "derivative_policy": {
            "finite_difference_backend_available": False,
            "accepted_derivative_backends": [
                "auto",
                "analytic",
                "cppad",
                "analytic_implicit",
                "cppad_implicit",
                "eigen_forward",
                "legacy_eigen_forward",
                "backend_unavailable",
            ],
            "unsupported_derivative_status": "backend_unavailable",
        },
        "phase_kind": kind,
        "chemical_success": bool(chemical.success),
        "phase_success": phase_success,
        "phase_problem_kind": str(getattr(phase, "problem_kind", kind)),
        "staged_feed": dict(z),
        "phase_equilibrium_diagnostics": phase_audit["diagnostics"],
    }
    return ReactiveStagedEquilibriumResult(
        success=success,
        message="converged" if success else "sequential reactive equilibrium did not converge",
        z=z,
        chemical=chemical,
        phase=phase,
        diagnostics=diagnostics,
    )


def _chemical_equilibrium_audit(chemical: ReactiveSpeciationResult) -> dict[str, Any]:
    x = dict(chemical.x)
    return {
        "reaction_coordinates": {
            "status": "implicit_in_solved_composition",
            "reaction_count": len(chemical.reaction_residuals),
            "named_reactions": list(chemical.named_reaction_residuals),
        },
        "element_balance_residuals": dict(chemical.mass_balance_residuals),
        "reaction_equilibrium_residuals": dict(chemical.named_reaction_residuals),
        "nonnegativity": {
            "status": "pass" if all(value >= 0.0 for value in x.values()) else "fail",
            "minimum_mole_fraction": min(x.values()) if x else 0.0,
        },
    }


def _phase_equilibrium_audit(phase: Any) -> dict[str, Any]:
    diagnostics = dict(getattr(phase, "diagnostics", {}) or {})
    split_detected = bool(getattr(phase, "split_detected", False))
    phases = tuple(getattr(phase, "phases", ()) or ())
    return {
        "phase_split": {
            "status": "split_detected" if split_detected else "single_phase_or_not_detected",
            "phase_count": len(phases),
            "phase_labels": [str(getattr(item, "label", "")) for item in phases],
            "phase_distance": _float_or_none(diagnostics.get("phase_distance")),
        },
        "fugacity_equality": {
            "status": "reported" if "fugacity_residual_norm" in diagnostics else "not_reported_by_phase_route",
            "fugacity_residual_norm": _float_or_none(diagnostics.get("fugacity_residual_norm")),
        },
        "material_balance_error": _float_or_none(diagnostics.get("material_balance_error")),
        "diagnostics": diagnostics,
    }


def _ascani_benchmark_attempt(kind: str) -> dict[str, Any]:
    route = str(kind).strip()
    if route in {"electrolyte_lle", "electrolyte_lle_flash", "electrolyte_lle_tp"}:
        return {
            "status": "supported_route_available",
            "reference": "Ascani electrolyte LLE transformed-basis route",
            "route": route,
        }
    if route in {"lle_flash", "lle_tp"}:
        return {
            "status": "not_applicable_to_neutral_route",
            "reference": "Ascani electrolyte benchmark requires an electrolyte LLE phase route",
            "route": route,
        }
    return {
        "status": "not_attempted_for_route",
        "reference": "Use phase_kind='electrolyte_lle' for Ascani electrolyte benchmark attempts.",
        "route": route,
    }


def _float_or_none(value: Any) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _normalize_workflow_options(workflow_options: Mapping[str, Any] | None) -> dict[str, str]:
    workflow = dict(workflow_options or {})
    reaction_constant_fitting = str(workflow.get("reaction_constant_fitting", "secondary_optional")).strip().lower()
    aliases = {
        "secondary": "secondary_optional",
        "optional": "secondary_optional",
        "lower_priority": "secondary_optional",
        "fixed": "secondary_optional",
    }
    reaction_constant_fitting = aliases.get(reaction_constant_fitting, reaction_constant_fitting)
    if reaction_constant_fitting in {"primary", "default", "required", "blocking"}:
        raise InputError("reaction-constant fitting is not a default in sequential workflows; keep constants fixed/literature first.")
    if reaction_constant_fitting != "secondary_optional":
        raise InputError(
            "workflow_options.reaction_constant_fitting must be 'secondary_optional' for sequential workflows."
        )
    return {"reaction_constant_fitting": reaction_constant_fitting}


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
