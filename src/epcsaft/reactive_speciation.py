"""Homogeneous reactive speciation helpers using ePC-SAFT activity states."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ._types import InputError, SolutionError
from .implicit_sensitivity import (
    ImplicitSolveResult,
    not_available_implicit_result,
    implicit_backend_for_residual_backend,
)

_SPECIATION_EVALUATION_ERRORS = (
    InputError,
    SolutionError,
    ValueError,
    RuntimeError,
    ArithmeticError,
    np.linalg.LinAlgError,
)
_COMPOSITION_NORMALIZATION_ERRORS = (
    InputError,
    TypeError,
    ValueError,
    ArithmeticError,
)

_REACTION_STANDARD_STATES = {
    "mole_fraction_activity": 0,
    "ideal_mole_fraction": 1,
    "concentration": 2,
    "thermodynamic_activity": 0,
    "molality": None,
    "apparent": None,
}

_REACTION_STANDARD_STATE_DEFAULT_BASIS = {
    "mole_fraction_activity": "mole_fraction",
    "ideal_mole_fraction": "mole_fraction",
    "concentration": "concentration",
    "thermodynamic_activity": "activity",
    "molality": "molality",
    "apparent": "apparent",
}

_REACTION_CONSTANT_KINDS = {"thermodynamic", "apparent", "fitted", "corrected"}
_REACTION_CONSTANT_FITTING_ROLES = {
    "fixed_input",
    "secondary_optional",
    "fitted_parameter",
    "regularized_correction",
}
_PYTHON_FIXED_POINT_STANDARD_STATES = {
    "mole_fraction_activity",
    "thermodynamic_activity",
    "concentration",
    "apparent",
}


@dataclass(frozen=True, slots=True)
class ReactionConstantConvention:
    """Explicit convention for interpreting a reaction equilibrium constant."""

    standard_state: str = "mole_fraction_activity"
    basis: str | None = None
    constant_kind: str = "thermodynamic"
    fitting_role: str = "fixed_input"
    correction_terms: Mapping[str, float] = field(default_factory=dict)
    source: str = ""

    def __post_init__(self) -> None:
        standard_state = str(self.standard_state).strip().lower()
        if standard_state not in _REACTION_STANDARD_STATES:
            supported = "', '".join(_REACTION_STANDARD_STATES)
            raise InputError(f"ReactionConstantConvention.standard_state must be one of '{supported}'.")
        basis = self.basis
        if basis is None:
            basis = _REACTION_STANDARD_STATE_DEFAULT_BASIS[standard_state]
        basis = str(basis).strip().lower()
        expected_basis = _REACTION_STANDARD_STATE_DEFAULT_BASIS[standard_state]
        if basis != expected_basis:
            raise InputError(
                "ReactionConstantConvention.basis is incompatible with standard_state "
                f"'{standard_state}'; expected '{expected_basis}'."
            )
        constant_kind = str(self.constant_kind).strip().lower()
        if constant_kind not in _REACTION_CONSTANT_KINDS:
            supported = "', '".join(sorted(_REACTION_CONSTANT_KINDS))
            raise InputError(f"ReactionConstantConvention.constant_kind must be one of '{supported}'.")
        fitting_role = str(self.fitting_role).strip().lower()
        if fitting_role not in _REACTION_CONSTANT_FITTING_ROLES:
            supported = "', '".join(sorted(_REACTION_CONSTANT_FITTING_ROLES))
            raise InputError(f"ReactionConstantConvention.fitting_role must be one of '{supported}'.")
        if constant_kind in {"fitted", "corrected"} and fitting_role == "fixed_input":
            raise InputError("Fitted or corrected reaction constants require an explicit non-fixed fitting_role.")
        corrections = {str(k): float(v) for k, v in dict(self.correction_terms).items()}
        for name, value in corrections.items():
            if not np.isfinite(value):
                raise InputError(f"ReactionConstantConvention.correction_terms['{name}'] must be finite.")
        object.__setattr__(self, "standard_state", standard_state)
        object.__setattr__(self, "basis", basis)
        object.__setattr__(self, "constant_kind", constant_kind)
        object.__setattr__(self, "fitting_role", fitting_role)
        object.__setattr__(self, "correction_terms", corrections)
        object.__setattr__(self, "source", str(self.source))

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any] | ReactionConstantConvention) -> ReactionConstantConvention:
        """Normalize a mapping or convention instance."""
        if isinstance(value, ReactionConstantConvention):
            return value
        return cls(**dict(value))

    @property
    def native_standard_state_code(self) -> int:
        """Return the native standard-state code or fail loudly when unsupported."""
        code = _REACTION_STANDARD_STATES[self.standard_state]
        if code is None:
            raise InputError(
                "not_available: reaction constant convention "
                f"'{self.standard_state}' is defined but is not supported by the native speciation backend."
            )
        return int(code)

    @property
    def requires_activity_coefficients(self) -> bool:
        """Whether the convention needs activity coefficients in the current native route."""
        return self.native_standard_state_code == _REACTION_STANDARD_STATES["mole_fraction_activity"]

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-like convention payload."""
        return {
            "standard_state": self.standard_state,
            "basis": self.basis,
            "constant_kind": self.constant_kind,
            "fitting_role": self.fitting_role,
            "correction_terms": dict(self.correction_terms),
            "source": self.source,
            "native_standard_state_code": _REACTION_STANDARD_STATES[self.standard_state],
        }


@dataclass(frozen=True, slots=True)
class ReactionDefinition:
    """One reaction residual definition for reactive speciation."""

    stoichiometry: Mapping[str, float]
    log_equilibrium_constant: float
    name: str = ""
    standard_state: str = "mole_fraction_activity"
    metadata: Mapping[str, Any] = field(default_factory=dict)
    convention: ReactionConstantConvention | Mapping[str, Any] | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "stoichiometry", {str(k): float(v) for k, v in self.stoichiometry.items()})
        object.__setattr__(self, "log_equilibrium_constant", float(self.log_equilibrium_constant))
        object.__setattr__(self, "name", str(self.name))
        metadata = {str(k): _json_like(v) for k, v in dict(self.metadata).items()}
        standard_state = str(self.standard_state).strip().lower()
        if self.convention is None:
            if standard_state not in _REACTION_STANDARD_STATES:
                supported = "', '".join(_REACTION_STANDARD_STATES)
                raise InputError(f"ReactionDefinition.standard_state must be one of '{supported}'.")
            convention = ReactionConstantConvention(
                standard_state=standard_state,
                constant_kind=str(metadata.get("constant_kind", "thermodynamic")),
                fitting_role=str(metadata.get("fitting_role", "fixed_input")),
                source=str(metadata.get("source", "")),
            )
        else:
            convention = ReactionConstantConvention.from_mapping(self.convention)
            if standard_state != "mole_fraction_activity" and standard_state != convention.standard_state:
                raise InputError("ReactionDefinition.standard_state is incompatible with convention.standard_state.")
        standard_state = convention.standard_state
        metadata.setdefault("constant_kind", convention.constant_kind)
        metadata.setdefault("fitting_role", convention.fitting_role)
        if convention.source:
            metadata.setdefault("source", convention.source)
        metadata.setdefault("constant_convention", convention.to_dict())
        object.__setattr__(self, "standard_state", standard_state)
        object.__setattr__(self, "metadata", metadata)
        object.__setattr__(self, "convention", convention)

    @classmethod
    def from_literature_constant(
        cls,
        stoichiometry: Mapping[str, float],
        *,
        log_equilibrium_constant: float,
        name: str = "",
        standard_state: str = "mole_fraction_activity",
        convention: ReactionConstantConvention | Mapping[str, Any] | None = None,
        source: str = "",
        metadata: Mapping[str, Any] | None = None,
    ) -> ReactionDefinition:
        """Build a reaction with an explicitly fixed literature equilibrium constant."""
        merged_metadata = dict(metadata or {})
        merged_metadata.setdefault("constant_source", "literature")
        merged_metadata.setdefault("fitting_role", "fixed_input")
        if source:
            merged_metadata.setdefault("source", str(source))
        return cls(
            stoichiometry=stoichiometry,
            log_equilibrium_constant=log_equilibrium_constant,
            name=name,
            standard_state=standard_state,
            metadata=merged_metadata,
            convention=convention,
        )

    @classmethod
    def from_fitted_constant(
        cls,
        stoichiometry: Mapping[str, float],
        *,
        log_equilibrium_constant: float,
        name: str = "",
        convention: ReactionConstantConvention | Mapping[str, Any] | None = None,
        metadata: Mapping[str, Any] | None = None,
    ) -> ReactionDefinition:
        """Build a reaction with an explicitly fitted equilibrium constant."""
        merged_metadata = dict(metadata or {})
        merged_metadata.setdefault("constant_source", "fit")
        if convention is None:
            convention = ReactionConstantConvention(
                constant_kind="fitted",
                fitting_role="fitted_parameter",
            )
        return cls(
            stoichiometry=stoichiometry,
            log_equilibrium_constant=log_equilibrium_constant,
            name=name,
            metadata=merged_metadata,
            convention=convention,
        )


@dataclass(frozen=True, slots=True)
class ReactiveSpeciationOptions:
    """Numerical controls for homogeneous reactive speciation."""

    max_iterations: int = 50
    tolerance: float = 1.0e-8
    damping: float = 0.5
    min_mole_fraction: float = 1.0e-14
    jacobian_backend: str = "auto"
    solver_backend: str = "auto"
    hessian_strategy: str = "gauss_newton"
    phase: str = "liq"
    return_best_effort: bool = False
    error_mode: str = "raise"
    activity_output: str = "auto"
    mass_tolerance: float | None = None
    charge_tolerance: float | None = None
    reaction_tolerance: float | None = None


@dataclass(frozen=True, slots=True)
class ReactiveSpeciationResult:
    """Structured result returned by reactive speciation solves."""

    success: bool
    message: str
    x: dict[str, float]
    activity_coefficients: dict[str, float]
    mass_balance_residuals: dict[str, float]
    charge_residual: float
    reaction_residuals: list[float]
    named_reaction_residuals: dict[str, float]
    state_failure_count: int
    diagnostics: dict[str, Any]
    continuation_state: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "success", bool(self.success))
        object.__setattr__(self, "message", str(self.message))
        object.__setattr__(self, "x", {str(k): float(v) for k, v in self.x.items()})
        object.__setattr__(
            self, "activity_coefficients", {str(k): float(v) for k, v in self.activity_coefficients.items()}
        )
        object.__setattr__(
            self, "mass_balance_residuals", {str(k): float(v) for k, v in self.mass_balance_residuals.items()}
        )
        object.__setattr__(self, "charge_residual", float(self.charge_residual))
        object.__setattr__(self, "reaction_residuals", [float(v) for v in self.reaction_residuals])
        object.__setattr__(
            self,
            "named_reaction_residuals",
            {str(k): float(v) for k, v in self.named_reaction_residuals.items()},
        )
        object.__setattr__(self, "state_failure_count", int(self.state_failure_count))
        object.__setattr__(self, "diagnostics", dict(self.diagnostics))
        object.__setattr__(self, "continuation_state", _json_like(dict(self.continuation_state)))

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-like result payload."""
        return {
            "success": self.success,
            "message": self.message,
            "x": dict(self.x),
            "activity_coefficients": dict(self.activity_coefficients),
            "mass_balance_residuals": dict(self.mass_balance_residuals),
            "charge_residual": self.charge_residual,
            "reaction_residuals": list(self.reaction_residuals),
            "named_reaction_residuals": dict(self.named_reaction_residuals),
            "state_failure_count": self.state_failure_count,
            "diagnostics": _json_like(self.diagnostics),
            "continuation_state": _json_like(self.continuation_state),
        }


def solve_reactive_speciation(
    *,
    species: Any,
    mixture_factory: Any,
    T: float,
    P: float,
    balances: Mapping[str, Mapping[str, float]],
    totals: Mapping[str, float],
    reactions: Any,
    initial_x: Any,
    options: ReactiveSpeciationOptions | None = None,
    warm_start: Any = None,
) -> ReactiveSpeciationResult:
    """Solve homogeneous activity-coupled reactive speciation."""
    opts = _normalize_options(options)
    labels = [str(label) for label in species]
    if not labels:
        raise InputError("species must include at least one label.")
    initial, initial_source = _initial_composition_from_inputs(
        initial_x=initial_x,
        warm_start=warm_start,
        ncomp=len(labels),
        min_value=opts.min_mole_fraction,
    )
    balance_matrix, total_vector, balance_names = _normalize_balances(labels, balances, totals)
    reaction_defs = _normalize_reactions(labels, reactions)
    if opts.solver_backend == "ipopt":
        from .ipopt_backend import solve_reactive_speciation_ipopt

        return solve_reactive_speciation_ipopt(
            species=labels,
            mixture_factory=mixture_factory,
            T=T,
            P=P,
            balance_matrix=balance_matrix,
            total_vector=total_vector,
            balance_names=balance_names,
            reactions=reaction_defs,
            initial_x=initial,
            initial_x_source=initial_source,
            options=opts,
        )
    if _requires_python_activity_fixed_point(reaction_defs):
        return _solve_reactive_speciation_activity_fixed_point(
            species=labels,
            mixture_factory=mixture_factory,
            T=T,
            P=P,
            balance_matrix=balance_matrix,
            total_vector=total_vector,
            balance_names=balance_names,
            reactions=reaction_defs,
            initial_x=initial,
            initial_x_source=initial_source,
            options=opts,
        )
    return _solve_reactive_speciation_native(
        species=labels,
        mixture_factory=mixture_factory,
        T=T,
        P=P,
        balance_matrix=balance_matrix,
        total_vector=total_vector,
        balance_names=balance_names,
        reactions=reaction_defs,
        initial_x=initial,
        initial_x_source=initial_source,
        options=opts,
    )


def _requires_python_activity_fixed_point(reactions: list[ReactionDefinition]) -> bool:
    return any(reaction.standard_state in _PYTHON_FIXED_POINT_STANDARD_STATES for reaction in reactions)


def _validate_mixture_for_species(mixture: Any, species: list[str]) -> None:
    native = getattr(mixture, "_native", None)
    if native is None:
        raise InputError("activity fixed-point speciation requires mixture_factory to return an ePCSAFTMixture.")
    if list(getattr(mixture, "species", species)) != species:
        raise InputError("activity fixed-point speciation requires mixture species to match the species argument.")


def _activity_fixed_point_payload(
    *,
    species: list[str],
    mixture: Any,
    T: float,
    P: float,
    x: np.ndarray,
    reactions: list[ReactionDefinition],
    options: ReactiveSpeciationOptions,
) -> dict[str, Any]:
    _validate_mixture_for_species(mixture, species)
    log_factors: list[np.ndarray] = []
    activity_coefficients = np.ones(len(species), dtype=float)
    activity_model = "ideal_or_apparent"
    activity_evaluation_count = 0
    density_solve_count = 0
    state = None
    for reaction in reactions:
        standard_state = reaction.standard_state
        if standard_state == "ideal_mole_fraction" or standard_state == "apparent":
            log_factors.append(np.zeros(len(species), dtype=float))
            continue
        if standard_state in {"mole_fraction_activity", "thermodynamic_activity"}:
            if state is None:
                state = mixture.state(T=T, P=P, x=x, phase=options.phase)
                density_solve_count += 1
            if _mixture_has_ionic_species(mixture):
                activity_model = "epcsaft_component_activity"
                gamma_map = state.activity_coefficient(species=species)
                activity_coefficients = np.asarray([gamma_map[label] for label in species], dtype=float)
            else:
                activity_model = "epcsaft_neutral_fugacity_activity"
                activity_coefficients = _neutral_fugacity_activity_coefficients(
                    mixture=mixture,
                    T=T,
                    P=P,
                    x=x,
                    phase=options.phase,
                    min_mole_fraction=options.min_mole_fraction,
                )
                density_solve_count += len(species)
            activity_evaluation_count += 1
            log_factors.append(np.log(np.clip(activity_coefficients, options.min_mole_fraction, None)))
            continue
        if standard_state == "concentration":
            if state is None:
                state = mixture.state(T=T, P=P, x=x, phase=options.phase)
                density_solve_count += 1
            molar_density = max(float(state.molar_density()), options.min_mole_fraction)
            activity_model = "concentration"
            log_factors.append(np.full(len(species), np.log(molar_density), dtype=float))
            continue
        raise InputError(
            "not_available: reaction constant convention "
            f"'{standard_state}' is defined but is not supported by activity fixed-point speciation."
        )
    exposes_activity_coefficients = activity_model.startswith("epcsaft_") or options.activity_output == "always"
    activity_coefficients_map = (
        {label: float(value) for label, value in zip(species, activity_coefficients)}
        if exposes_activity_coefficients
        else {}
    )
    return {
        "log_factors": log_factors,
        "activity_coefficients": activity_coefficients.tolist(),
        "activity_coefficients_map": activity_coefficients_map,
        "activity_model": activity_model,
        "activity_evaluation_count": activity_evaluation_count,
        "density_solve_count": density_solve_count,
    }


def _neutral_fugacity_activity_coefficients(
    *,
    mixture: Any,
    T: float,
    P: float,
    x: np.ndarray,
    phase: str,
    min_mole_fraction: float,
) -> np.ndarray:
    state = mixture.state(T=T, P=P, x=x, phase=phase)
    ln_phi = np.asarray(state.fugacity_coefficient(natural_log=True), dtype=float)
    ln_gamma: list[float] = []
    ncomp = len(x)
    for idx in range(ncomp):
        x_ref = np.full(ncomp, max(min_mole_fraction, 1.0e-14), dtype=float)
        x_ref[idx] = 1.0 - x_ref[0] * float(ncomp - 1)
        x_ref = x_ref / float(np.sum(x_ref))
        ref = mixture.state(T=T, P=P, x=x_ref, phase=phase)
        ln_phi_ref = np.asarray(ref.fugacity_coefficient(natural_log=True), dtype=float)
        ln_gamma.append(float(np.clip(ln_phi[idx] - ln_phi_ref[idx], -700.0, 700.0)))
    return np.exp(np.asarray(ln_gamma, dtype=float))


def _mixture_has_ionic_species(mixture: Any) -> bool:
    charges = np.asarray(getattr(mixture, "parameters", {}).get("z", []), dtype=float).flatten()
    return bool(charges.size and np.any(np.abs(charges) > 1.0e-12))


def _idealized_reactions_for_activity_fixed_point(
    species: list[str],
    reactions: list[ReactionDefinition],
    log_factors: list[np.ndarray],
) -> list[ReactionDefinition]:
    out: list[ReactionDefinition] = []
    for reaction, factors in zip(reactions, log_factors):
        correction = _reaction_log_factor_correction(species, reaction, factors)
        metadata = dict(reaction.metadata)
        metadata["activity_fixed_point_source_standard_state"] = reaction.standard_state
        metadata["activity_fixed_point_constant_kind"] = reaction.convention.constant_kind
        convention = ReactionConstantConvention(
            standard_state="ideal_mole_fraction",
            constant_kind=reaction.convention.constant_kind,
            fitting_role=reaction.convention.fitting_role,
            source=reaction.convention.source,
        )
        out.append(
            ReactionDefinition(
                stoichiometry=reaction.stoichiometry,
                log_equilibrium_constant=reaction.log_equilibrium_constant - correction,
                name=reaction.name,
                metadata=metadata,
                convention=convention,
            )
        )
    return out


def _reaction_log_factor_correction(
    species: list[str],
    reaction: ReactionDefinition,
    factors: np.ndarray,
) -> float:
    return float(
        sum(float(coeff) * float(factors[species.index(label)]) for label, coeff in reaction.stoichiometry.items())
    )


def _reaction_residuals_with_log_factors(
    *,
    species: list[str],
    x: np.ndarray,
    reactions: list[ReactionDefinition],
    log_factors: list[np.ndarray],
    min_mole_fraction: float,
) -> list[float]:
    ln_x = np.log(np.clip(np.asarray(x, dtype=float), min_mole_fraction, None))
    residuals: list[float] = []
    for reaction, factors in zip(reactions, log_factors):
        value = -float(reaction.log_equilibrium_constant)
        for label, coeff in reaction.stoichiometry.items():
            value += float(coeff) * float(ln_x[species.index(label)] + factors[species.index(label)])
        residuals.append(float(value))
    return residuals


def _charge_residual_from_mixture(mixture: Any, x: np.ndarray) -> float:
    charges = np.asarray(getattr(mixture, "parameters", {}).get("z", []), dtype=float).flatten()
    if charges.size != np.asarray(x).size:
        return 0.0
    return float(charges @ np.asarray(x, dtype=float))


def _damped_composition_update(
    current: np.ndarray,
    candidate: np.ndarray,
    damping: float,
    min_mole_fraction: float,
) -> np.ndarray:
    updated = (1.0 - float(damping)) * np.asarray(current, dtype=float) + float(damping) * np.asarray(
        candidate, dtype=float
    )
    return _normalize_composition(np.clip(updated, min_mole_fraction, None), len(updated), min_mole_fraction)


def _activity_fixed_point_diagnostics(
    *,
    reactions: list[ReactionDefinition],
    options: ReactiveSpeciationOptions,
    inner: ReactiveSpeciationResult,
    activity_payload: dict[str, Any],
    iteration: int,
    history: list[float],
    initial_x_source: str,
    mass_tolerance: float,
    charge_tolerance: float,
    reaction_tolerance: float,
    mass_residual_norm: float,
    charge_residual: float,
    reaction_residual_norm: float,
    residual_norm: float,
    named_reaction_residuals: dict[str, float],
    x_map: dict[str, float],
    activity_evaluation_count: int,
    density_solve_count: int,
    state_failure_count: int,
) -> dict[str, Any]:
    diagnostics = dict(inner.diagnostics)
    _normalize_reactive_derivative_diagnostics(diagnostics)
    diagnostics.update(
        {
            "solver_language": "python",
            "native_entrypoint": "_solve_chemical_equilibrium_native",
            "problem_class": "homogeneous_chemical_equilibrium",
            "activity_fixed_point": True,
            "activity_fixed_point_updates": int(iteration),
            "activity_model": activity_payload["activity_model"],
            "activity_basis": _reaction_standard_state_summary(reactions),
            "reaction_standard_states": [reaction.standard_state for reaction in reactions],
            "reaction_constant_conventions": _reaction_constant_conventions(reactions),
            "reaction_constant_sources": _reaction_constant_sources(reactions),
            "reaction_constant_policy": "fixed_literature_constants_first",
            "backend": "native_plus_python_activity_fixed_point",
            "thermodynamic_backend": "epcsaft_state_activity_chemical_potential_api",
            "selected_solver_backend": "activity_fixed_point_native_inner",
            "solver_selection_reason": "nonideal_or_apparent_standard_state",
            "requested_solver_backend": str(options.solver_backend),
            "requested_hessian_strategy": str(options.hessian_strategy),
            "requested_jacobian_backend": str(options.jacobian_backend),
            "jacobian_backend": "analytic",
            "derivative_backend": "analytic",
            "derivative_status": "analytic",
            "derivative_available": True,
            "jacobian_available": True,
            "not_available_reason": "",
            "jacobian_fallback_used": False,
            "hessian_fallback_used": False,
            "numerical_derivative_backend_available": False,
            "activity_derivative_in_jacobian": False,
            "activity_derivative_policy": "not_used_by_fixed_point_outer_iteration",
            "residual_norm": float(residual_norm),
            "residual_norm_by_block": {
                "material_balance": float(mass_residual_norm),
                "charge_balance": float(abs(charge_residual)),
                "reaction_affinity": float(reaction_residual_norm),
            },
            "history": [float(value) for value in history],
            "iterations": int(iteration),
            "state_failure_count": int(state_failure_count),
            "activity_evaluation_count": int(activity_evaluation_count),
            "density_solve_count": int(density_solve_count),
            "activity_coefficients_evaluated": bool(activity_payload["activity_coefficients_map"]),
            "mass_tolerance": float(mass_tolerance),
            "charge_tolerance": float(charge_tolerance),
            "reaction_tolerance": float(reaction_tolerance),
            "mass_residual_norm": float(mass_residual_norm),
            "charge_residual_abs": float(abs(charge_residual)),
            "reaction_residual_norm": float(reaction_residual_norm),
            "named_reaction_residuals": dict(named_reaction_residuals),
            "best_x": dict(x_map),
            "best_activity_coefficients": dict(activity_payload["activity_coefficients_map"]),
            "initial_x_source": str(initial_x_source),
            "continuation_used": str(initial_x_source) != "initial_x",
        }
    )
    diagnostics["derivative_backend_by_block"].update(
        {
            "reaction_residual_jacobian": "analytic",
            "activity_or_fugacity_state": "analytic",
            "activity_fixed_point_outer_iteration": "analytic",
        }
    )
    return diagnostics


def _solve_reactive_speciation_activity_fixed_point(
    *,
    species: list[str],
    mixture_factory: Any,
    T: float,
    P: float,
    balance_matrix: np.ndarray,
    total_vector: np.ndarray,
    balance_names: list[str],
    reactions: list[ReactionDefinition],
    initial_x: np.ndarray,
    options: ReactiveSpeciationOptions,
    initial_x_source: str = "initial_x",
) -> ReactiveSpeciationResult:
    if options.jacobian_backend in {"autodiff", "cppad"}:
        raise InputError(
            "not_available: explicit autodiff/CppAD chemical-equilibrium residual Jacobian is not "
            "implemented for activity-fixed-point speciation."
        )
    if options.solver_backend == "ipopt":
        raise InputError(
            "not_available: solver_backend='ipopt' is not implemented for activity-fixed-point speciation."
        )

    scalar_result = _try_scalar_binary_activity_solve(
        species=species,
        mixture_factory=mixture_factory,
        T=T,
        P=P,
        balance_matrix=balance_matrix,
        total_vector=total_vector,
        balance_names=balance_names,
        reactions=reactions,
        initial_x=initial_x,
        options=options,
        initial_x_source=initial_x_source,
    )
    if scalar_result is not None:
        return scalar_result

    current = np.asarray(initial_x, dtype=float)
    best_result: ReactiveSpeciationResult | None = None
    best_residual_norm = float("inf")
    history: list[float] = []
    activity_evaluation_count = 0
    density_solve_count = 0
    state_failure_count = 0
    last_activity_payload: dict[str, Any] = {}
    fixed_point_success = False

    for iteration in range(options.max_iterations + 1):
        try:
            mixture = mixture_factory(current, T, P)
            _validate_mixture_for_species(mixture, species)
            activity_payload = _activity_fixed_point_payload(
                species=species,
                mixture=mixture,
                T=T,
                P=P,
                x=current,
                reactions=reactions,
                options=options,
            )
        except _SPECIATION_EVALUATION_ERRORS:
            state_failure_count += 1
            raise
        activity_evaluation_count += int(activity_payload.get("activity_evaluation_count", 0))
        density_solve_count += int(activity_payload.get("density_solve_count", 0))
        last_activity_payload = activity_payload
        ideal_reactions = _idealized_reactions_for_activity_fixed_point(
            species,
            reactions,
            activity_payload["log_factors"],
        )
        inner_options = ReactiveSpeciationOptions(
            max_iterations=options.max_iterations,
            tolerance=options.tolerance,
            damping=options.damping,
            min_mole_fraction=options.min_mole_fraction,
            jacobian_backend="auto",
            solver_backend="auto",
            hessian_strategy=options.hessian_strategy,
            phase=options.phase,
            return_best_effort=True,
            error_mode="result",
            activity_output="never",
            mass_tolerance=options.mass_tolerance,
            charge_tolerance=options.charge_tolerance,
            reaction_tolerance=options.reaction_tolerance,
        )
        inner = _solve_reactive_speciation_native(
            species=species,
            mixture_factory=lambda x, t, p, _mixture=mixture: _mixture,
            T=T,
            P=P,
            balance_matrix=balance_matrix,
            total_vector=total_vector,
            balance_names=balance_names,
            reactions=ideal_reactions,
            initial_x=current,
            initial_x_source="activity_fixed_point",
            options=inner_options,
        )
        candidate = _normalize_composition(list(inner.x.values()), len(species), options.min_mole_fraction)
        full_payload = _activity_fixed_point_payload(
            species=species,
            mixture=mixture_factory(candidate, T, P),
            T=T,
            P=P,
            x=candidate,
            reactions=reactions,
            options=options,
        )
        activity_evaluation_count += int(full_payload.get("activity_evaluation_count", 0))
        density_solve_count += int(full_payload.get("density_solve_count", 0))
        last_activity_payload = full_payload
        reaction_residuals = _reaction_residuals_with_log_factors(
            species=species,
            x=candidate,
            reactions=reactions,
            log_factors=full_payload["log_factors"],
            min_mole_fraction=options.min_mole_fraction,
        )
        x_map = {label: float(value) for label, value in zip(species, candidate)}
        mass_balance_residuals = {
            name: float(value)
            for name, value in zip(balance_names, np.asarray(balance_matrix, dtype=float) @ candidate - total_vector)
        }
        charge_residual = _charge_residual_from_mixture(mixture_factory(candidate, T, P), candidate)
        named_reaction_residuals = _named_reaction_residuals(reactions, reaction_residuals)
        mass_tolerance = options.mass_tolerance if options.mass_tolerance is not None else options.tolerance
        charge_tolerance = options.charge_tolerance if options.charge_tolerance is not None else options.tolerance
        reaction_tolerance = options.reaction_tolerance if options.reaction_tolerance is not None else options.tolerance
        mass_residual_norm = float(max((abs(value) for value in mass_balance_residuals.values()), default=0.0))
        reaction_residual_norm = float(max((abs(value) for value in reaction_residuals), default=0.0))
        residual_norm = max(mass_residual_norm, abs(charge_residual), reaction_residual_norm)
        history.append(residual_norm)
        residual_family_success = (
            mass_residual_norm <= mass_tolerance
            and abs(charge_residual) <= charge_tolerance
            and reaction_residual_norm <= reaction_tolerance
        )
        diagnostics = _activity_fixed_point_diagnostics(
            reactions=reactions,
            options=options,
            inner=inner,
            activity_payload=full_payload,
            iteration=iteration,
            history=history,
            initial_x_source=initial_x_source,
            mass_tolerance=mass_tolerance,
            charge_tolerance=charge_tolerance,
            reaction_tolerance=reaction_tolerance,
            mass_residual_norm=mass_residual_norm,
            charge_residual=charge_residual,
            reaction_residual_norm=reaction_residual_norm,
            residual_norm=residual_norm,
            named_reaction_residuals=named_reaction_residuals,
            x_map=x_map,
            activity_evaluation_count=activity_evaluation_count,
            density_solve_count=density_solve_count,
            state_failure_count=state_failure_count,
        )
        best_result = ReactiveSpeciationResult(
            success=bool(inner.success and residual_family_success),
            message="converged" if inner.success and residual_family_success else "activity fixed-point did not converge",
            x=x_map,
            activity_coefficients=full_payload["activity_coefficients_map"],
            mass_balance_residuals=mass_balance_residuals,
            charge_residual=charge_residual,
            reaction_residuals=reaction_residuals,
            named_reaction_residuals=named_reaction_residuals,
            state_failure_count=state_failure_count,
            diagnostics=diagnostics,
            continuation_state=_continuation_state(x=x_map, T=T, P=P, diagnostics=diagnostics),
        )
        if residual_norm < best_residual_norm:
            best_residual_norm = residual_norm
        if best_result.success:
            fixed_point_success = True
            break
        if iteration == options.max_iterations:
            break
        current = _damped_composition_update(current, candidate, options.damping, options.min_mole_fraction)

    if best_result is None:
        raise SolutionError("activity fixed-point speciation did not produce a result", _json_like(last_activity_payload))
    if not fixed_point_success and not options.return_best_effort:
        raise SolutionError(best_result.message, _json_like(best_result.diagnostics))
    return best_result


def _try_scalar_binary_activity_solve(
    *,
    species: list[str],
    mixture_factory: Any,
    T: float,
    P: float,
    balance_matrix: np.ndarray,
    total_vector: np.ndarray,
    balance_names: list[str],
    reactions: list[ReactionDefinition],
    initial_x: np.ndarray,
    options: ReactiveSpeciationOptions,
    initial_x_source: str,
) -> ReactiveSpeciationResult | None:
    if len(species) != 2 or len(reactions) != 1 or balance_matrix.shape != (1, 2):
        return None
    if not np.allclose(balance_matrix[0], np.ones(2), rtol=0.0, atol=1.0e-14):
        return None
    total = float(total_vector[0])
    if not np.isfinite(total) or abs(total - 1.0) > 1.0e-12:
        return None
    reaction = reactions[0]
    if any(label not in species for label in reaction.stoichiometry):
        return None

    def evaluate(x0: float) -> tuple[float, dict[str, Any]]:
        x = _normalize_composition([x0, 1.0 - x0], 2, options.min_mole_fraction)
        mixture = mixture_factory(x, T, P)
        payload = _activity_fixed_point_payload(
            species=species,
            mixture=mixture,
            T=T,
            P=P,
            x=x,
            reactions=reactions,
            options=options,
        )
        residual = _reaction_residuals_with_log_factors(
            species=species,
            x=x,
            reactions=reactions,
            log_factors=payload["log_factors"],
            min_mole_fraction=options.min_mole_fraction,
        )[0]
        return float(residual), payload

    lower = max(options.min_mole_fraction, 1.0e-12)
    upper = 1.0 - lower
    grid = np.linspace(lower, upper, 101)
    values: list[tuple[float, float, dict[str, Any]]] = []
    for point in grid:
        try:
            residual, payload = evaluate(float(point))
        except _SPECIATION_EVALUATION_ERRORS:
            continue
        if np.isfinite(residual):
            values.append((float(point), residual, payload))
    if not values:
        return None
    brackets: list[tuple[float, float, float, float]] = []
    for (x_left, f_left, _), (x_right, f_right, _) in zip(values, values[1:]):
        if f_left == 0.0:
            brackets.append((x_left, x_left, f_left, f_left))
        elif f_left * f_right <= 0.0:
            brackets.append((x_left, x_right, f_left, f_right))
    if not brackets:
        return None
    seed = float(initial_x[0])
    left, right, f_left, f_right = min(brackets, key=lambda item: abs(0.5 * (item[0] + item[1]) - seed))
    payload = values[0][2]
    history: list[float] = []
    root = left
    residual = f_left
    for iteration in range(max(1, options.max_iterations * 4)):
        root = 0.5 * (left + right)
        residual, payload = evaluate(root)
        history.append(abs(residual))
        if abs(residual) <= options.tolerance or abs(right - left) <= options.min_mole_fraction:
            break
        if f_left * residual <= 0.0:
            right = root
            f_right = residual
        else:
            left = root
            f_left = residual
    del f_right
    x = _normalize_composition([root, 1.0 - root], 2, options.min_mole_fraction)
    mixture = mixture_factory(x, T, P)
    mass_balance_residuals = {balance_names[0]: float(np.sum(x) - total)}
    charge_residual = _charge_residual_from_mixture(mixture, x)
    named_reaction_residuals = _named_reaction_residuals(reactions, [residual])
    x_map = {label: float(value) for label, value in zip(species, x)}
    mass_tolerance = options.mass_tolerance if options.mass_tolerance is not None else options.tolerance
    charge_tolerance = options.charge_tolerance if options.charge_tolerance is not None else options.tolerance
    reaction_tolerance = options.reaction_tolerance if options.reaction_tolerance is not None else options.tolerance
    mass_residual_norm = abs(mass_balance_residuals[balance_names[0]])
    reaction_residual_norm = abs(float(residual))
    residual_norm = max(mass_residual_norm, abs(charge_residual), reaction_residual_norm)
    inner = ReactiveSpeciationResult(
        success=True,
        message="scalar binary activity solve",
        x=x_map,
        activity_coefficients={},
        mass_balance_residuals=mass_balance_residuals,
        charge_residual=charge_residual,
        reaction_residuals=[residual],
        named_reaction_residuals=named_reaction_residuals,
        state_failure_count=0,
        diagnostics={"derivative_backend": "analytic", "residual_norm": residual_norm},
    )
    diagnostics = _activity_fixed_point_diagnostics(
        reactions=reactions,
        options=options,
        inner=inner,
        activity_payload=payload,
        iteration=len(history),
        history=history,
        initial_x_source=initial_x_source,
        mass_tolerance=mass_tolerance,
        charge_tolerance=charge_tolerance,
        reaction_tolerance=reaction_tolerance,
        mass_residual_norm=mass_residual_norm,
        charge_residual=charge_residual,
        reaction_residual_norm=reaction_residual_norm,
        residual_norm=residual_norm,
        named_reaction_residuals=named_reaction_residuals,
        x_map=x_map,
        activity_evaluation_count=len(history),
        density_solve_count=len(history) * (1 + len(species)),
        state_failure_count=0,
    )
    diagnostics["selected_solver_backend"] = "scalar_binary_activity_bracket"
    diagnostics["solver_selection_reason"] = "binary_reaction_scalar_activity_residual"
    success = residual_norm <= max(mass_tolerance, charge_tolerance, reaction_tolerance)
    result = ReactiveSpeciationResult(
        success=success,
        message="converged" if success else "scalar binary activity solve did not converge",
        x=x_map,
        activity_coefficients=payload["activity_coefficients_map"],
        mass_balance_residuals=mass_balance_residuals,
        charge_residual=charge_residual,
        reaction_residuals=[residual],
        named_reaction_residuals=named_reaction_residuals,
        state_failure_count=0,
        diagnostics=diagnostics,
        continuation_state=_continuation_state(x=x_map, T=T, P=P, diagnostics=diagnostics),
    )
    if not success and not options.return_best_effort:
        raise SolutionError(result.message, _json_like(result.diagnostics))
    return result


def solve_reactive_speciation_sweep(
    *,
    species: Any,
    mixture_factory: Any,
    points: Any,
    balances: Mapping[str, Mapping[str, float]],
    reactions: Any,
    options: ReactiveSpeciationOptions | None = None,
    continuation: str = "auto",
) -> list[ReactiveSpeciationResult]:
    """Solve an ordered reactive-speciation sweep with fixed-shape results."""

    opts = _normalize_options(options)
    mode = str(continuation).strip().lower()
    if mode not in {"auto", "none"}:
        raise InputError("continuation must be 'auto' or 'none'.")
    labels = [str(label) for label in species]
    if not labels:
        raise InputError("species must include at least one label.")
    reaction_defs = _normalize_reactions(labels, reactions)
    results: list[ReactiveSpeciationResult] = []
    warm_start: dict[str, Any] | None = None
    for index, point in enumerate(points):
        try:
            if "T" not in point or "P" not in point or "totals" not in point:
                raise InputError("Each reactive speciation sweep point requires T, P, and totals.")
            point_warm_start = warm_start if mode == "auto" else None
            result = solve_reactive_speciation(
                species=labels,
                mixture_factory=mixture_factory,
                T=float(point["T"]),
                P=float(point["P"]),
                balances=balances,
                totals=point["totals"],
                reactions=reaction_defs,
                initial_x=point.get("initial_x"),
                options=opts,
                warm_start=point_warm_start,
            )
        except Exception as exc:
            if opts.error_mode != "result":
                raise
            result = _structured_failure_result(
                species=labels,
                point=point,
                index=index,
                exc=exc,
                options=opts,
            )
        results.append(result)
        if result.success:
            warm_start = dict(result.continuation_state)
    return results


def _normalize_options(options: ReactiveSpeciationOptions | None) -> ReactiveSpeciationOptions:
    if options is None:
        return ReactiveSpeciationOptions()
    if not isinstance(options, ReactiveSpeciationOptions):
        raise InputError("options must be a ReactiveSpeciationOptions instance.")
    if options.max_iterations < 0:
        raise InputError("ReactiveSpeciationOptions.max_iterations must be non-negative.")
    if options.tolerance <= 0.0 or options.min_mole_fraction <= 0.0:
        raise InputError("ReactiveSpeciationOptions tolerances must be positive.")
    if not (0.0 < options.damping <= 1.0):
        raise InputError("ReactiveSpeciationOptions.damping must be in (0, 1].")
    if not isinstance(options.return_best_effort, bool):
        raise InputError("ReactiveSpeciationOptions.return_best_effort must be a bool.")
    error_mode = str(options.error_mode).strip().lower()
    if error_mode not in {"raise", "result"}:
        raise InputError("ReactiveSpeciationOptions.error_mode must be 'raise' or 'result'.")
    activity_output = str(options.activity_output).strip().lower()
    if activity_output not in {"auto", "always", "never"}:
        raise InputError("ReactiveSpeciationOptions.activity_output must be 'auto', 'always', or 'never'.")
    return_best_effort = bool(options.return_best_effort or error_mode == "result")
    jacobian_backend = str(options.jacobian_backend).strip().lower()
    if jacobian_backend == "analytic":
        jacobian_backend = "auto"
    if jacobian_backend not in {"auto", "autodiff", "cppad"}:
        raise InputError(
            "ReactiveSpeciationOptions.jacobian_backend must be 'auto', 'autodiff', 'analytic', or 'cppad'."
        )
    solver_backend = str(options.solver_backend).strip().lower()
    if solver_backend not in {"auto", "newton", "ipopt"}:
        raise InputError("ReactiveSpeciationOptions.solver_backend must be 'auto', 'newton', or 'ipopt'.")
    hessian_strategy = str(options.hessian_strategy).strip().lower()
    hessian_aliases = {"gn": "gauss_newton", "gauss-newton": "gauss_newton", "bfgs": "lbfgs"}
    hessian_strategy = hessian_aliases.get(hessian_strategy, hessian_strategy)
    if hessian_strategy not in {"gauss_newton", "lbfgs"}:
        raise InputError("ReactiveSpeciationOptions.hessian_strategy must be 'gauss_newton' or 'lbfgs'.")
    for name in ("mass_tolerance", "charge_tolerance", "reaction_tolerance"):
        value = getattr(options, name)
        if value is not None and value <= 0.0:
            raise InputError(f"ReactiveSpeciationOptions.{name} must be positive when provided.")
    if (
        jacobian_backend == options.jacobian_backend
        and solver_backend == options.solver_backend
        and hessian_strategy == options.hessian_strategy
        and error_mode == options.error_mode
        and activity_output == options.activity_output
        and return_best_effort == options.return_best_effort
    ):
        return options
    return ReactiveSpeciationOptions(
        max_iterations=options.max_iterations,
        tolerance=options.tolerance,
        damping=options.damping,
        min_mole_fraction=options.min_mole_fraction,
        jacobian_backend=jacobian_backend,
        solver_backend=solver_backend,
        hessian_strategy=hessian_strategy,
        phase=options.phase,
        return_best_effort=return_best_effort,
        error_mode=error_mode,
        activity_output=activity_output,
        mass_tolerance=options.mass_tolerance,
        charge_tolerance=options.charge_tolerance,
        reaction_tolerance=options.reaction_tolerance,
    )


def _solve_reactive_speciation_native(
    *,
    species: list[str],
    mixture_factory: Any,
    T: float,
    P: float,
    balance_matrix: np.ndarray,
    total_vector: np.ndarray,
    balance_names: list[str],
    reactions: list[ReactionDefinition],
    initial_x: np.ndarray,
    options: ReactiveSpeciationOptions,
    initial_x_source: str = "initial_x",
) -> ReactiveSpeciationResult:
    from . import _core

    mixture = mixture_factory(initial_x, T, P)
    native = getattr(mixture, "_native", None)
    if native is None:
        raise InputError("native reactive speciation backend requires mixture_factory to return an ePCSAFTMixture.")
    if list(getattr(mixture, "species", species)) != species:
        raise InputError("native reactive speciation backend requires mixture species to match the species argument.")
    reaction_matrix = np.asarray(
        [[float(reaction.stoichiometry.get(label, 0.0)) for label in species] for reaction in reactions],
        dtype=float,
    )
    request = {
        "T": float(T),
        "P": float(P),
        "initial_x": np.asarray(initial_x, dtype=float).tolist(),
        "balance_matrix": np.asarray(balance_matrix, dtype=float).reshape(-1).tolist(),
        "balance_rows": int(balance_matrix.shape[0]),
        "total_vector": np.asarray(total_vector, dtype=float).tolist(),
        "reaction_stoichiometry": reaction_matrix.reshape(-1).tolist(),
        "reaction_rows": int(reaction_matrix.shape[0]),
        "log_equilibrium_constants": [float(reaction.log_equilibrium_constant) for reaction in reactions],
        "reaction_standard_states": [reaction.convention.native_standard_state_code for reaction in reactions],
        "options": {
            "max_iterations": int(options.max_iterations),
            "tolerance": float(options.tolerance),
            "damping": float(options.damping),
            "min_mole_fraction": float(options.min_mole_fraction),
            "jacobian_backend": str(options.jacobian_backend),
            "solver_backend": str(options.solver_backend),
            "hessian_strategy": str(options.hessian_strategy),
            "phase": str(options.phase),
            "activity_output": str(options.activity_output),
        },
    }
    try:
        payload = _core._solve_chemical_equilibrium_native(native, request)
    except _core.NativeValueError as exc:
        raise InputError(str(exc)) from exc
    x = {label: float(value) for label, value in zip(species, payload["composition"])}
    activity_coefficients = {label: float(value) for label, value in zip(species, payload["activity_coefficients"])}
    mass_balance_residuals = {
        name: float(value) for name, value in zip(balance_names, payload["mass_balance_residuals"])
    }
    reaction_residuals = [float(value) for value in payload["reaction_residuals"]]
    named_reaction_residuals = _named_reaction_residuals(reactions, reaction_residuals)
    mass_tolerance = options.mass_tolerance if options.mass_tolerance is not None else options.tolerance
    charge_tolerance = options.charge_tolerance if options.charge_tolerance is not None else options.tolerance
    reaction_tolerance = options.reaction_tolerance if options.reaction_tolerance is not None else options.tolerance
    mass_residual_norm = float(max((abs(value) for value in mass_balance_residuals.values()), default=0.0))
    reaction_residual_norm = float(max((abs(value) for value in reaction_residuals), default=0.0))
    charge_residual = float(payload["charge_residual"])
    residual_family_success = (
        mass_residual_norm <= mass_tolerance
        and abs(charge_residual) <= charge_tolerance
        and reaction_residual_norm <= reaction_tolerance
    )
    diagnostics = dict(payload["diagnostics"])
    _normalize_reactive_derivative_diagnostics(diagnostics)
    activity_basis = _reaction_standard_state_summary(reactions)
    handoff = dict(diagnostics.get("phase_equilibrium_handoff", {}))
    handoff.setdefault("composition", [float(value) for value in payload["composition"]])
    handoff.setdefault("activity_coefficients", [float(value) for value in payload["activity_coefficients"]])
    handoff["composition_map"] = dict(x)
    handoff["activity_coefficients_map"] = dict(activity_coefficients)
    handoff["activity_basis"] = activity_basis
    diagnostics["phase_equilibrium_handoff"] = handoff
    diagnostics["reaction_standard_states"] = [reaction.standard_state for reaction in reactions]
    diagnostics["reaction_constant_conventions"] = _reaction_constant_conventions(reactions)
    diagnostics["reaction_constant_sources"] = _reaction_constant_sources(reactions)
    diagnostics["reaction_constant_policy"] = "fixed_literature_constants_first"
    diagnostics.update(
        {
            "activity_basis": activity_basis,
            "success": bool(payload["success"] and residual_family_success),
            "native_success": bool(payload["success"]),
            "residual_family_success": bool(residual_family_success),
            "message": str(payload["message"]),
            "backend": "native",
            "activity_output": str(options.activity_output),
            "initial_x_source": str(initial_x_source),
            "continuation_used": str(initial_x_source) != "initial_x",
            "requested_solver_backend": str(options.solver_backend),
            "requested_hessian_strategy": str(options.hessian_strategy),
            "selected_solver_backend": "native",
            "solver_selection_reason": "default_native" if options.solver_backend == "auto" else "explicit_request",
            "default_auto_uses_ipopt": False,
            "mass_tolerance": float(mass_tolerance),
            "charge_tolerance": float(charge_tolerance),
            "reaction_tolerance": float(reaction_tolerance),
            "mass_residual_norm": mass_residual_norm,
            "charge_residual_abs": abs(charge_residual),
            "reaction_residual_norm": reaction_residual_norm,
            "named_reaction_residuals": dict(named_reaction_residuals),
            "best_x": dict(x),
            "best_activity_coefficients": dict(activity_coefficients),
        }
    )
    success = bool(payload["success"] and residual_family_success)
    message = str(payload["message"])
    if bool(payload["success"]) and not residual_family_success:
        message = "reactive speciation residual family tolerances were not met"
        diagnostics["message"] = message
    result = ReactiveSpeciationResult(
        success=success,
        message=message,
        x=x,
        activity_coefficients=activity_coefficients,
        mass_balance_residuals=mass_balance_residuals,
        charge_residual=charge_residual,
        reaction_residuals=reaction_residuals,
        named_reaction_residuals=named_reaction_residuals,
        state_failure_count=int(diagnostics.get("state_failure_count", 0)),
        diagnostics=diagnostics,
        continuation_state=_continuation_state(
            x=x,
            T=T,
            P=P,
            diagnostics=diagnostics,
        ),
    )
    if not success and not options.return_best_effort:
        raise SolutionError(message, _json_like(diagnostics))
    return result


def _normalize_reactive_derivative_diagnostics(diagnostics: dict[str, Any]) -> None:
    derivative_backend = str(diagnostics.get("derivative_backend", "not_available"))
    solved_state_backend = implicit_backend_for_residual_backend(derivative_backend)
    diagnostics.setdefault("thermodynamic_backend", "epcsaft_state_activity_chemical_potential_api")
    diagnostics.setdefault("solver_backend", diagnostics.get("nonlinear_solver", "native_newton"))
    diagnostics.setdefault("derivative_backend", derivative_backend)
    diagnostics.setdefault("derivative_status", derivative_backend)
    diagnostics.setdefault("jacobian_fallback_used", False)
    diagnostics.setdefault("hessian_fallback_used", False)
    if derivative_backend == "not_available":
        diagnostics.setdefault(
            "not_available_reason",
            "not_available: reactive speciation sensitivities are not implemented for this route.",
        )
        diagnostics.setdefault("derivative_available", False)
        diagnostics.setdefault("jacobian_available", False)
    if "residual_norm" in diagnostics:
        diagnostics.setdefault("residual_norm_by_block", {"reactive_speciation": float(diagnostics["residual_norm"])})
    else:
        diagnostics.setdefault("residual_norm_by_block", {})
    diagnostics.setdefault("solved_internal_states", ["reactive_speciation_log_amounts", "density_roots"])
    diagnostics.setdefault(
        "derivative_backend_by_block",
        {
            "reaction_residual_jacobian": derivative_backend,
            "density_root": "not_available",
            "activity_or_fugacity_state": "analytic" if derivative_backend == "analytic" else derivative_backend,
        },
    )
    diagnostics.setdefault("implicit_sensitivity_blocks", [])
    diagnostics.setdefault("best_state_available", True)
    diagnostics.setdefault("best_state", {"source": "native_chemical_equilibrium_result"})
    diagnostics.setdefault("row_failure_count", int(diagnostics.get("state_failure_count", 0)))
    diagnostics.setdefault("association_solver_status", "not_available_if_active")
    diagnostics.setdefault(
        "derivative_policy",
        {
            "numerical_derivative_backend_available": False,
            "accepted_derivative_backends": [
                "auto",
                "analytic",
                "cppad",
                "analytic_implicit",
                "cppad_implicit",
                "not_available",
            ],
            "unsupported_derivative_status": "not_available",
        },
    )
    diagnostics.setdefault(
        "solved_state_derivative_blocks",
        [
            "association_site_fractions",
            "reactive_speciation_variables",
            "density_roots",
            "bubble_pressure_roots",
            "phase_equilibrium_variables",
        ],
    )
    diagnostics["derivative_backend_by_block"].setdefault(
        "reactive_speciation_variables",
        solved_state_backend,
    )
    diagnostics["derivative_backend_by_block"].setdefault("association_site_fractions", "not_available")
    if solved_state_backend == "not_available":
        reactive_implicit_result = not_available_implicit_result(
            reason="reactive speciation implicit sensitivities are unavailable for this residual backend.",
            diagnostics={"residual_backend": derivative_backend},
        )
    else:
        reactive_implicit_result = ImplicitSolveResult(
            state=(),
            residual=(),
            jacobians={},
            sensitivity=(),
            backend=solved_state_backend,
            status="residual_jacobian_available",
            diagnostics={
                "residual_backend": derivative_backend,
                "sensitivity_scope": "generic implicit solved-state contract",
            },
        )
    diagnostics.setdefault(
        "implicit_solve_results",
        {
            "reactive_speciation_variables": reactive_implicit_result.to_dict(),
            "association_site_fractions": not_available_implicit_result(
                reason="association site-fraction implicit sensitivities are unavailable for this reactive speciation route.",
                diagnostics={"residual_backend": "not_available"},
            ).to_dict(),
        },
    )


def _named_reaction_residuals(reactions: list[ReactionDefinition], reaction_residuals: list[float]) -> dict[str, float]:
    names: list[str] = []
    counts: dict[str, int] = {}
    for index, reaction in enumerate(reactions):
        base = reaction.name.strip() if reaction.name.strip() else f"reaction_{index}"
        count = counts.get(base, 0)
        counts[base] = count + 1
        names.append(base if count == 0 else f"{base}_{count}")
    return {name: float(value) for name, value in zip(names, reaction_residuals)}


def _reaction_constant_sources(reactions: list[ReactionDefinition]) -> dict[str, str]:
    sources: dict[str, str] = {}
    counts: dict[str, int] = {}
    for index, reaction in enumerate(reactions):
        base = reaction.name.strip() if reaction.name.strip() else f"reaction_{index}"
        count = counts.get(base, 0)
        counts[base] = count + 1
        name = base if count == 0 else f"{base}_{count}"
        sources[name] = str(reaction.metadata.get("constant_source", "fixed_input"))
    return sources


def _reaction_constant_conventions(reactions: list[ReactionDefinition]) -> dict[str, dict[str, Any]]:
    conventions: dict[str, dict[str, Any]] = {}
    counts: dict[str, int] = {}
    for index, reaction in enumerate(reactions):
        base = reaction.name.strip() if reaction.name.strip() else f"reaction_{index}"
        count = counts.get(base, 0)
        counts[base] = count + 1
        name = base if count == 0 else f"{base}_{count}"
        conventions[name] = reaction.convention.to_dict()
    return conventions


def _reaction_standard_state_summary(reactions: list[ReactionDefinition]) -> str:
    standard_states = [reaction.standard_state for reaction in reactions]
    if not standard_states:
        return "mole_fraction"
    first = standard_states[0]
    if any(value != first for value in standard_states[1:]):
        return "mixed_standard_state"
    if first == "mole_fraction_activity":
        return "mole_fraction"
    return first


def _initial_composition_from_inputs(
    *,
    initial_x: Any,
    warm_start: Any,
    ncomp: int,
    min_value: float,
) -> tuple[np.ndarray, str]:
    if warm_start is not None:
        if isinstance(warm_start, Mapping):
            for key in ("composition", "x"):
                if key in warm_start:
                    value = warm_start[key]
                    if isinstance(value, Mapping):
                        value = list(value.values())
                    return _normalize_composition(value, ncomp, min_value), "previous_successful_result"
        return _normalize_composition(warm_start, ncomp, min_value), "warm_start"
    if initial_x is None:
        raise InputError("initial_x is required when no warm_start composition is supplied.")
    return _normalize_composition(initial_x, ncomp, min_value), "initial_x"


def _continuation_state(
    *, x: Mapping[str, float], T: float, P: float, diagnostics: Mapping[str, Any]
) -> dict[str, Any]:
    composition = {str(label): float(value) for label, value in x.items()}
    return {
        "composition": composition,
        "T": float(T),
        "P": float(P),
        "density_solve_count": int(diagnostics.get("density_solve_count", 0)),
        "activity_evaluation_count": int(diagnostics.get("activity_evaluation_count", 0)),
    }


def _structured_failure_result(
    *,
    species: list[str],
    point: Mapping[str, Any],
    index: int,
    exc: Exception,
    options: ReactiveSpeciationOptions,
) -> ReactiveSpeciationResult:
    try:
        x_array = _normalize_composition(point.get("initial_x"), len(species), options.min_mole_fraction)
    except _COMPOSITION_NORMALIZATION_ERRORS:
        x_array = np.full(len(species), 1.0 / max(len(species), 1), dtype=float)
    x = {label: float(value) for label, value in zip(species, x_array)}
    diagnostics = {
        "success": False,
        "structured_failure": True,
        "sweep_index": int(index),
        "message": str(exc),
        "exception_type": type(exc).__name__,
        "backend": "python_sweep",
        "selected_solver_backend": "not_run",
        "initial_x_source": "failure_payload",
        "continuation_used": False,
    }
    return ReactiveSpeciationResult(
        success=False,
        message=str(exc),
        x=x,
        activity_coefficients={},
        mass_balance_residuals={},
        charge_residual=0.0,
        reaction_residuals=[],
        named_reaction_residuals={},
        state_failure_count=0,
        diagnostics=diagnostics,
        continuation_state={},
    )


def _normalize_composition(value: Any, ncomp: int, min_value: float) -> np.ndarray:
    x = np.asarray(value, dtype=float).flatten()
    if x.size != int(ncomp):
        raise InputError("initial_x length must match species length.")
    if np.any(~np.isfinite(x)) or np.any(x < -min_value):
        raise InputError("initial_x values must be finite and non-negative.")
    x = np.clip(x, 0.0, None)
    total = float(np.sum(x))
    if total <= 0.0:
        raise InputError("initial_x must have a positive sum.")
    return x / total


def _normalize_balances(
    species: list[str],
    balances: Mapping[str, Mapping[str, float]],
    totals: Mapping[str, float],
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    rows: list[list[float]] = []
    names: list[str] = []
    total_values: list[float] = []
    for name, coeffs in balances.items():
        row = [0.0 for _ in species]
        for label, coeff in coeffs.items():
            if str(label) not in species:
                raise InputError(f"Unknown species '{label}' in balance '{name}'.")
            row[species.index(str(label))] = float(coeff)
        if name not in totals:
            raise InputError(f"Missing total for balance '{name}'.")
        rows.append(row)
        names.append(str(name))
        total_values.append(float(totals[name]))
    if not rows:
        raise InputError("At least one material balance is required.")
    return np.asarray(rows, dtype=float), np.asarray(total_values, dtype=float), names


def _normalize_reactions(species: list[str], reactions: Any) -> list[ReactionDefinition]:
    out: list[ReactionDefinition] = []
    for reaction in reactions:
        if not isinstance(reaction, ReactionDefinition):
            raise InputError("reactions must contain ReactionDefinition instances.")
        for label in reaction.stoichiometry:
            if label not in species:
                raise InputError(f"Unknown species '{label}' in reaction stoichiometry.")
        out.append(reaction)
    return out


def _json_like(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(key): _json_like(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_like(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value
