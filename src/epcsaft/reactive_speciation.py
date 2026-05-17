"""Homogeneous reactive speciation helpers using ePC-SAFT activity states."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ._types import InputError, SolutionError
from .implicit_sensitivity import (
    ImplicitSolveResult,
    implicit_backend_for_residual_backend,
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


def _raise_native_ipopt_reactive_speciation_required() -> None:
    raise InputError(
        "reactive_speciation requires a native Ipopt homogeneous reactive-speciation NLP route. "
        "No package-owned alternate homogeneous reactive-speciation solver is available for this public route."
    )


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
                "unsupported reaction constant convention "
                f"'{self.standard_state}' is defined but is not supported by the native speciation backend."
            )
        return int(code)

    @property
    def requires_activity_coefficients(self) -> bool:
        """Whether the convention needs activity coefficients."""
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
    phase_stoichiometry: Mapping[str, Mapping[str, float]] | None = None

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
        phase_stoichiometry = None
        if self.phase_stoichiometry is not None:
            phase_stoichiometry = {
                str(phase): {str(label): float(coeff) for label, coeff in dict(coeffs).items()}
                for phase, coeffs in dict(self.phase_stoichiometry).items()
            }
            metadata.setdefault("reaction_phase_scope", "phase_tagged_cross_phase")
        object.__setattr__(self, "standard_state", standard_state)
        object.__setattr__(self, "metadata", metadata)
        object.__setattr__(self, "convention", convention)
        object.__setattr__(self, "phase_stoichiometry", phase_stoichiometry)

    @classmethod
    def from_literature_constant(
        cls,
        stoichiometry: Mapping[str, float],
        *,
        log_equilibrium_constant: float,
        name: str = "",
        standard_state: str = "mole_fraction_activity",
        convention: ReactionConstantConvention | Mapping[str, Any] | None = None,
        phase_stoichiometry: Mapping[str, Mapping[str, float]] | None = None,
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
            phase_stoichiometry=phase_stoichiometry,
        )

    @classmethod
    def from_fitted_constant(
        cls,
        stoichiometry: Mapping[str, float],
        *,
        log_equilibrium_constant: float,
        name: str = "",
        convention: ReactionConstantConvention | Mapping[str, Any] | None = None,
        phase_stoichiometry: Mapping[str, Mapping[str, float]] | None = None,
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
            phase_stoichiometry=phase_stoichiometry,
        )


@dataclass(frozen=True, slots=True)
class ReactiveSpeciationOptions:
    """Numerical controls for homogeneous reactive speciation."""

    max_iterations: int = 50
    tolerance: float = 1.0e-8
    min_mole_fraction: float = 1.0e-14
    jacobian_backend: str = "auto"
    solver_backend: str = "auto"
    phase: str = "liq"
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
) -> ReactiveSpeciationResult:
    """Solve homogeneous reactive speciation through the accepted native route."""
    opts = _normalize_options(options)
    requested_solver_backend = opts.solver_backend
    labels = [str(label) for label in species]
    if not labels:
        raise InputError("species must include at least one label.")
    initial = _normalize_required_initial_composition(initial_x, len(labels), opts.min_mole_fraction)
    balance_matrix, total_vector, balance_names = _normalize_balances(labels, balances, totals)
    reaction_defs = _normalize_reactions(labels, reactions)
    if opts.jacobian_backend == "cppad":
        raise InputError("CppAD reactive-speciation Jacobians require an implemented CppAD NLP derivative route.")
    if opts.solver_backend == "auto":
        opts = _options_with_solver_backend(opts, "ipopt")
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
        initial_x_source="initial_x",
        options=opts,
        requested_solver_backend=requested_solver_backend,
    )


def solve_reactive_speciation_sweep(
    *,
    species: Any,
    mixture_factory: Any,
    points: Any,
    balances: Mapping[str, Mapping[str, float]],
    reactions: Any,
    options: ReactiveSpeciationOptions | None = None,
) -> list[ReactiveSpeciationResult]:
    """Solve an ordered reactive-speciation sweep with fixed-shape results."""

    opts = _normalize_options(options)
    labels = [str(label) for label in species]
    if not labels:
        raise InputError("species must include at least one label.")
    reaction_defs = _normalize_reactions(labels, reactions)
    results: list[ReactiveSpeciationResult] = []
    for index, point in enumerate(points):
        try:
            if "T" not in point or "P" not in point or "totals" not in point:
                raise InputError("Each reactive speciation sweep point requires T, P, and totals.")
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
    error_mode = str(options.error_mode).strip().lower()
    if error_mode not in {"raise", "result"}:
        raise InputError("ReactiveSpeciationOptions.error_mode must be 'raise' or 'result'.")
    activity_output = str(options.activity_output).strip().lower()
    if activity_output not in {"auto", "always", "never"}:
        raise InputError("ReactiveSpeciationOptions.activity_output must be 'auto', 'always', or 'never'.")
    jacobian_backend = str(options.jacobian_backend).strip().lower()
    if jacobian_backend == "analytic":
        jacobian_backend = "auto"
    if jacobian_backend not in {"auto", "cppad"}:
        raise InputError(
            "ReactiveSpeciationOptions.jacobian_backend must be 'auto', 'analytic', or 'cppad'."
        )
    solver_backend = str(options.solver_backend).strip().lower()
    if solver_backend not in {"auto", "ipopt"}:
        raise InputError("ReactiveSpeciationOptions.solver_backend must be 'auto' or 'ipopt'.")
    for name in ("mass_tolerance", "charge_tolerance", "reaction_tolerance"):
        value = getattr(options, name)
        if value is not None and value <= 0.0:
            raise InputError(f"ReactiveSpeciationOptions.{name} must be positive when provided.")
    if (
        jacobian_backend == options.jacobian_backend
        and solver_backend == options.solver_backend
        and error_mode == options.error_mode
        and activity_output == options.activity_output
    ):
        return options
    return ReactiveSpeciationOptions(
        max_iterations=options.max_iterations,
        tolerance=options.tolerance,
        min_mole_fraction=options.min_mole_fraction,
        jacobian_backend=jacobian_backend,
        solver_backend=solver_backend,
        phase=options.phase,
        error_mode=error_mode,
        activity_output=activity_output,
        mass_tolerance=options.mass_tolerance,
        charge_tolerance=options.charge_tolerance,
        reaction_tolerance=options.reaction_tolerance,
    )

def _options_with_solver_backend(
    options: ReactiveSpeciationOptions,
    solver_backend: str,
) -> ReactiveSpeciationOptions:
    return ReactiveSpeciationOptions(
        max_iterations=options.max_iterations,
        tolerance=options.tolerance,
        min_mole_fraction=options.min_mole_fraction,
        jacobian_backend=options.jacobian_backend,
        solver_backend=solver_backend,
        phase=options.phase,
        error_mode=options.error_mode,
        activity_output=options.activity_output,
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
    requested_solver_backend: str | None = None,
) -> ReactiveSpeciationResult:
    if options.solver_backend != "ipopt":
        _raise_native_ipopt_reactive_speciation_required()
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
            "min_mole_fraction": float(options.min_mole_fraction),
            "jacobian_backend": str(options.jacobian_backend),
            "solver_backend": str(options.solver_backend),
            "phase": str(options.phase),
            "activity_output": str(options.activity_output),
        },
    }
    try:
        payload = _core._solve_chemical_equilibrium_native(native, request)
    except _core.NativeValueError as exc:
        raise InputError(str(exc)) from exc
    except _core.NativeSolutionError as exc:
        raise SolutionError(str(exc)) from exc
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
            "requested_solver_backend": str(requested_solver_backend or options.solver_backend),
            "selected_solver_backend": str(diagnostics.get("selected_solver_backend", "native_ipopt")),
            "solver_selection_reason": str(
                "auto_selected_native_ipopt"
                if requested_solver_backend == "auto"
                else diagnostics.get("solver_selection_reason", "explicit_request")
            ),
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
    )
    if not success and options.error_mode != "result":
        raise SolutionError(message, _json_like(diagnostics))
    return result


def _normalize_reactive_derivative_diagnostics(diagnostics: dict[str, Any]) -> None:
    derivative_backend = str(diagnostics.get("derivative_backend", "")).strip().lower()
    if not derivative_backend:
        raise InputError("native reactive speciation diagnostics must include a derivative backend.")
    solved_state_backend = implicit_backend_for_residual_backend(derivative_backend)
    diagnostics.setdefault("thermodynamic_backend", "epcsaft_state_activity_chemical_potential_api")
    diagnostics.setdefault(
        "solver_backend",
        diagnostics.get("selected_solver_backend", diagnostics.get("nonlinear_solver", "native_ipopt")),
    )
    diagnostics.setdefault("derivative_backend", derivative_backend)
    if "residual_norm" in diagnostics:
        diagnostics.setdefault("residual_norm_by_block", {"reactive_speciation": float(diagnostics["residual_norm"])})
    else:
        diagnostics.setdefault("residual_norm_by_block", {})
    diagnostics.setdefault(
        "derivative_backend_by_block",
        {
            "reaction_residual_jacobian": derivative_backend,
            "activity_or_fugacity_state": "analytic" if derivative_backend == "analytic" else derivative_backend,
        },
    )
    diagnostics.setdefault("implicit_sensitivity_blocks", [])
    diagnostics.setdefault(
        "derivative_policy",
        {
            "accepted_derivative_backends": [
                "auto",
                "analytic",
                "cppad",
                "analytic_implicit",
                "cppad_implicit",
            ],
            "unsupported_derivative_behavior": "raise",
        },
    )
    diagnostics["derivative_backend_by_block"].setdefault(
        "reactive_speciation_variables",
        solved_state_backend,
    )
    if _has_native_reactive_implicit_payload(diagnostics):
        reactive_implicit_result = _native_reactive_implicit_result(
            diagnostics,
            solved_state_backend=solved_state_backend,
            derivative_backend=derivative_backend,
        )
        blocks = list(diagnostics.get("implicit_sensitivity_blocks", []))
        if "reactive_speciation_variables" not in blocks:
            blocks.append("reactive_speciation_variables")
        diagnostics["implicit_sensitivity_blocks"] = blocks
        diagnostics.setdefault(
            "implicit_solve_results",
            {
                "reactive_speciation_variables": reactive_implicit_result.to_dict(),
            },
        )


def _has_native_reactive_implicit_payload(diagnostics: dict[str, Any]) -> bool:
    rows = int(diagnostics.get("reactive_speciation_residual_rows", 0))
    state_size = int(diagnostics.get("reactive_speciation_state_size", 0))
    parameter_size = int(diagnostics.get("reactive_speciation_parameter_size", 0))
    return rows > 0 and state_size > 0 and parameter_size > 0


def _native_reactive_implicit_result(
    diagnostics: dict[str, Any],
    *,
    solved_state_backend: str,
    derivative_backend: str,
) -> ImplicitSolveResult:
    rows = int(diagnostics.get("reactive_speciation_residual_rows", 0))
    state_size = int(diagnostics.get("reactive_speciation_state_size", 0))
    parameter_size = int(diagnostics.get("reactive_speciation_parameter_size", 0))
    if rows <= 0 or state_size <= 0 or parameter_size <= 0:
        raise InputError("native reactive speciation did not return implicit sensitivity matrices.")
    state = np.asarray(diagnostics.get("reactive_speciation_state", []), dtype=float)
    residual = np.asarray(diagnostics.get("reactive_speciation_residual", []), dtype=float)
    residual_state = np.asarray(
        diagnostics.get("reactive_speciation_residual_state_jacobian_row_major", []),
        dtype=float,
    ).reshape(rows, state_size)
    residual_parameter = np.asarray(
        diagnostics.get("reactive_speciation_residual_parameter_jacobian_row_major", []),
        dtype=float,
    ).reshape(rows, parameter_size)
    sensitivity = np.asarray(
        diagnostics.get("reactive_speciation_log_amount_sensitivity_to_log_k_row_major", []),
        dtype=float,
    ).reshape(state_size, parameter_size)
    return ImplicitSolveResult(
        state=state,
        residual=residual,
        jacobians={
            "residual_state": residual_state,
            "residual_parameter": residual_parameter,
        },
        sensitivity=sensitivity,
        backend=solved_state_backend,
        diagnostics={
            "residual_backend": derivative_backend,
            "sensitivity_scope": "log_amount_response_to_reaction_constants",
            "parameter": str(diagnostics.get("reactive_speciation_sensitivity_parameter", "log_equilibrium_constants")),
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


def _normalize_required_initial_composition(initial_x: Any, ncomp: int, min_value: float) -> np.ndarray:
    if initial_x is None:
        raise InputError("initial_x is required.")
    return _normalize_composition(initial_x, ncomp, min_value)


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
        if reaction.phase_stoichiometry is not None:
            for phase, coeffs in reaction.phase_stoichiometry.items():
                if not coeffs:
                    raise InputError(f"phase_stoichiometry for phase '{phase}' must include at least one species.")
                for label in coeffs:
                    if label not in species:
                        raise InputError(f"Unknown species '{label}' in reaction phase_stoichiometry.")
        _ = reaction.convention.native_standard_state_code
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
