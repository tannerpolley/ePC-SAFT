from copy import deepcopy

import numpy as np

from . import _core
from ._types import ActivityCoefficientResult, InputError, SolutionError, phase_to_int, vector_to_array

_NATIVE_CALL_ERRORS = (
    InputError,
    SolutionError,
    ValueError,
    RuntimeError,
    ArithmeticError,
    getattr(_core, "NativeValueError", ValueError),
    getattr(_core, "NativeSolutionError", SolutionError),
)
_DERIVATIVE_VALUE_ERRORS = (
    InputError,
    SolutionError,
    KeyError,
    ValueError,
    RuntimeError,
    ArithmeticError,
)

STATE_METHOD_ALIAS_MAP = {
    "pressure": "p",
    "density": "rho",
    "molar_density": "rho_molar",
    "mass_density": "rho_mass",
    "compressibility_factor": "z",
    "residual_helmholtz": "ares",
    "temperature_derivative_residual_helmholtz": "dadt",
    "composition_derivative_residual_helmholtz": "dadx",
    "residual_enthalpy": "hres",
    "residual_entropy": "sres",
    "residual_gibbs": "gres",
    "residual_chemical_potential": "mures",
    "activity_coefficient": "gamma",
    "fugacity_coefficient": "fugcoef",
    "relative_permittivity": "epsr",
    "osmotic_coefficient": "osmotic_coef",
    "state_diagnostics": "diag",
    "solvation_free_energy": "gsolv",
}
_STATE_METHOD_ALIAS_LOOKUP = {alias: name for name, alias in STATE_METHOD_ALIAS_MAP.items()}
_CONTRIBUTION_NAMES = ("hc", "disp", "assoc", "ion", "born")
_CONTRIBUTION_PUBLIC_NAMES = {
    "hc": "hard_chain",
    "disp": "dispersion",
    "assoc": "association",
    "ion": "ionic",
    "born": "born",
}
_GAS_CONSTANT = 8.31446261815324


def _result_with_route_diagnostics(result, route):
    diagnostics = dict(result.diagnostics)
    diagnostics["equilibrium_route"] = route["route"]
    diagnostics["route_reason"] = route["reason"]
    return type(result)(
        backend=result.backend,
        problem_kind=result.problem_kind,
        phases=result.phases,
        stable=result.stable,
        split_detected=result.split_detected,
        diagnostics=diagnostics,
    )


def _state_construction_error_message(T, x, phase, ncomp, mode_name, variable_name, variable_value, exc):
    return (
        f"{mode_name}-based state solve failed for "
        f"T={float(T)}, {variable_name}={float(variable_value)}, phase={phase}, ncomp={ncomp}, x={x.tolist()}: {exc}"
    )


def _sum_vector_terms(terms):
    total = None
    for name in _CONTRIBUTION_NAMES:
        value = np.asarray(terms[name], dtype=float)
        total = value.copy() if total is None else total + value
    return total


def _public_contribution_terms(terms):
    out = {}
    for internal, public in _CONTRIBUTION_PUBLIC_NAMES.items():
        out[public] = terms[internal]
    if "ideal" in terms:
        out["ideal"] = terms["ideal"]
    return out


def _scalar_terms_dict(terms):
    return {
        "hc": float(terms.hc),
        "disp": float(terms.disp),
        "assoc": float(terms.assoc),
        "ion": float(terms.ion),
        "born": float(terms.born),
        "total": float(terms.total),
    }


def _vector_terms_dict(terms, expected_size, label):
    out = {}
    blocks = {
        "hc": vector_to_array(terms.hc),
        "disp": vector_to_array(terms.disp),
        "assoc": vector_to_array(terms.assoc),
        "ion": vector_to_array(terms.ion),
        "born": vector_to_array(terms.born),
        "total": vector_to_array(terms.total),
    }
    for name, arr in blocks.items():
        arr = np.asarray(arr, dtype=float)
        if arr.size != expected_size:
            raise SolutionError(
                f"Unexpected {label} payload size for {name}: expected {int(expected_size)}, got {int(arr.size)}."
            )
        out[name] = arr
    return out


def _derivative_result_payload(
    *,
    supported,
    backend,
    message,
    value,
    jacobian,
    shape=None,
    outputs=None,
    variables=None,
    **extra,
):
    value_arr = np.asarray(value, dtype=float)
    jac_arr = np.asarray(jacobian, dtype=float)
    if shape is None:
        if jac_arr.ndim == 2:
            shape = [int(jac_arr.shape[0]), int(jac_arr.shape[1])]
        elif value_arr.ndim == 0:
            shape = [1, int(jac_arr.size)]
        else:
            shape = [int(value_arr.size), int(jac_arr.size // max(int(value_arr.size), 1))]
    if outputs is None:
        outputs = extra.get("output_order", extra.get("component_order", ()))
    if variables is None:
        variables = extra.get("variable_order", extra.get("parameter_order", ()))
    payload = {
        "supported": bool(supported),
        "backend": str(backend),
        "derivative_backend": str(backend),
        "message": str(message),
        "value": value_arr,
        "jacobian": jac_arr.reshape(tuple(shape)) if bool(supported) and jac_arr.size else jac_arr,
        "outputs": [str(item) for item in (outputs or [])],
        "variables": [str(item) for item in (variables or [])],
        "shape": [int(shape[0]), int(shape[1])],
    }
    payload.update(extra)
    return payload


def _backend_from_contribution_details(details, *, available=True, reason=""):
    if not available:
        _unsupported_derivative(reason or "composition derivative contribution backend is unsupported")
    labels = {str(v) for v in dict(details).values()}
    if "unsupported" in labels:
        _unsupported_derivative(reason or "composition derivative contribution backend is unsupported")
    if "cppad" in labels:
        return "cppad", "mixed CppAD/analytic derivative contributions"
    if "cppad_implicit" in labels:
        return "cppad_implicit", "mixed implicit CppAD/analytic derivative contributions"
    if "analytic_implicit" in labels:
        return "analytic_implicit", "mixed implicit analytic/analytic derivative contributions"
    return "analytic", "analytic derivative contributions"


def _unsupported_derivative(message):
    text = str(message)
    if not text.startswith("unsupported:"):
        text = f"unsupported: {text}"
    raise InputError(text)


def _as_rule_number(value, aliases, default):
    if value is None:
        return int(default)
    if isinstance(value, (int, np.integer)):
        return int(value)
    token = str(value).strip().lower()
    if token in aliases:
        return int(aliases[token])
    if token.isdigit() or (token.startswith("-") and token[1:].isdigit()):
        return int(token)
    return int(default)


def _state_rel_perm_rule_and_mode(params):
    z = np.asarray(params.get("z", []), dtype=float).flatten()
    default_rule = 1 if z.size else 0
    elec_model = params.get("elec_model") if isinstance(params.get("elec_model"), dict) else {}
    rel_perm = elec_model.get("rel_perm", {}) if isinstance(elec_model.get("rel_perm", {}), dict) else {}
    rule_alias = {
        "constant": 0,
        "rule0": 0,
        "linear": 1,
        "linear-molefraction": 1,
        "linear-mixing-mole": 1,
        "rule1": 1,
        "linear-massfraction": 2,
        "linear-mixing-weight": 2,
        "rule2": 2,
        "combined": 3,
        "rule3": 3,
        "empirical": 4,
        "rule4": 4,
        "rule5": 5,
        "rule6": 6,
        "aqueous-organic": 8,
        "aqueous_organic": 8,
        "mixed-aqueous-organic": 8,
        "mixed_aqueous_organic": 8,
        "rule8": 8,
        "salt-free-massfraction": 9,
        "salt_free_massfraction": 9,
        "salt-free-solvent-massfraction": 9,
        "salt_free_solvent_weight": 9,
        "rule9": 9,
    }
    diff_alias = {
        "auto": 3,
        "automatic": 3,
        "analytic": 0,
        "analytical": 0,
        "cppad": 2,
    }
    rule = _as_rule_number(rel_perm.get("rule"), rule_alias, default_rule)
    mode = _as_rule_number(rel_perm.get("differential_mode"), diff_alias, 3)
    return rule, mode


def _relative_permittivity_backend(params):
    rule, mode = _state_rel_perm_rule_and_mode(params)
    if mode == 2 or (mode == 3 and rule == 8):
        return "cppad"
    return "analytic"


class ePCSAFTMixture:
    """Native-backed ePC-SAFT parameter model and state factory."""

    def __init__(self, params=None, species=None):
        """Create a mixture from a resolved parameter payload."""
        self._native = None
        self._params = None
        self._species = None
        if params is not None:
            self._init_from_params(params, species)

    @classmethod
    def from_params(cls, params, species=None):
        """Construct a mixture from an already-resolved parameter dict."""
        if hasattr(params, "to_legacy_dict"):
            if species is None and hasattr(params, "components"):
                species = list(params.components)
            params = params.to_legacy_dict()
        return cls(params=params, species=species)

    @classmethod
    def from_dataset(cls, dataset_name, species, x, T, user_options=None):
        """Construct a mixture by resolving packaged dataset parameters."""
        from .parameters import get_prop_dict

        params = get_prop_dict(dataset_name, species, x, T, user_options=user_options)
        return cls(params=params, species=species)

    def _init_from_params(self, params, species=None):
        """Initialize the native mixture from a normalized parameter payload."""
        if hasattr(params, "to_legacy_dict"):
            if species is None and hasattr(params, "components"):
                species = list(params.components)
            params = params.to_legacy_dict()
        params = check_association(params)
        cppargs = create_struct(params)
        self._native = _core.NativeMixture(cppargs)
        self._params = params
        if species is None:
            ncomp = int(np.asarray(params["m"], dtype=float).size)
            self._species = [str(i) for i in range(ncomp)]
        else:
            self._species = [str(s) for s in species]

    @property
    def species(self):
        """Return the species labels in the mixture order."""
        return list(self._species)

    @property
    def parameters(self):
        """Return a deep copy of the resolved parameter payload."""
        return deepcopy(self._params)

    @property
    def ncomp(self):
        """Return the number of components in the mixture."""
        return int(self._native.ncomp())

    def clear_runtime_caches(self):
        """Clear internal runtime caches used for repeated state/reference evaluations."""
        self._native.clear_runtime_caches()

    def reset_runtime_cache_stats(self):
        """Reset runtime cache hit/rejection counters without clearing cached values."""
        self._native.reset_runtime_cache_stats()

    def runtime_cache_stats(self):
        """Return native runtime cache counters for profiling and validation."""
        return {
            "reference_state_cache_hits": int(self._native.reference_state_cache_hits()),
            "reference_state_cache_misses": int(self._native.reference_state_cache_misses()),
            "density_warm_start_hits": int(self._native.density_warm_start_hits()),
            "density_warm_start_rejections": int(self._native.density_warm_start_rejections()),
        }

    def state(self, T, x, P=None, rho=None, phase="liq", rho_guess=None, rho_seed=None):
        """Create an immutable thermodynamic state for the mixture.

        States built from pressure resolve and cache density during construction.
        ``rho_guess``/``rho_seed`` may seed that pressure-density solve but does not replace
        pressure closure.
        """
        if (P is None) == (rho is None):
            raise InputError("Provide exactly one of P or rho when constructing a state.")
        if rho_guess is not None and rho_seed is not None:
            raise InputError("Provide only one of rho_guess or rho_seed.")
        if rho_seed is not None:
            rho_guess = rho_seed
        return ePCSAFTState(self, T, x, P=P, rho=rho, phase=phase, rho_guess=rho_guess)

    def check_density(self, T, x, P, rho, *, phase="liq", rtol=1.0e-6, atol=1.0e-3):
        """Return pressure-consistency diagnostics for an externally supplied density."""
        if float(P) <= 0.0:
            raise InputError("P must be positive when checking density consistency.")
        if not np.isfinite(float(rtol)) or float(rtol) < 0.0:
            raise InputError("rtol must be finite and non-negative.")
        if not np.isfinite(float(atol)) or float(atol) < 0.0:
            raise InputError("atol must be finite and non-negative.")
        state = self.state(T=T, x=x, rho=rho, phase=phase)
        pressure_target = float(P)
        pressure_from_density = float(state.pressure())
        pressure_residual = pressure_from_density - pressure_target
        scale = max(abs(pressure_target), float(atol))
        relative_pressure_residual = abs(pressure_residual) / scale
        within_tolerance = bool(abs(pressure_residual) <= float(atol) + float(rtol) * abs(pressure_target))
        return {
            "density": float(rho),
            "pressure_target": pressure_target,
            "pressure_from_density": pressure_from_density,
            "pressure_residual": float(pressure_residual),
            "relative_pressure_residual": float(relative_pressure_residual),
            "within_tolerance": within_tolerance,
            "state": state,
        }

    def solve_equilibrium(self, problem):
        """Solve a typed equilibrium problem object against this mixture."""
        solve = getattr(problem, "solve", None)
        if solve is None or not callable(solve):
            raise InputError("solve_equilibrium requires an EquilibriumProblem object with a solve(mixture) method.")
        return solve(self)

    def flash_tp(self, T, P, z, *, options=None):
        """Validate a neutral TP flash request with explicit thermodynamic-method naming."""
        from .equilibrium import tp_flash

        return tp_flash(self, T=T, P=P, z=z, options=options)

    def bubble_p(self, T, x, *, options=None):
        """Solve a neutral bubble pressure from liquid composition at fixed temperature."""
        from .equilibrium import bubble_p

        return bubble_p(self, T=T, x=x, options=options)

    def bubble_t(self, P, x, *, options=None):
        """Solve a neutral bubble temperature from liquid composition at fixed pressure."""
        from .equilibrium import bubble_t

        return bubble_t(self, P=P, x=x, options=options)

    def dew_p(self, T, y, *, options=None):
        """Solve a neutral dew pressure from vapor composition at fixed temperature."""
        from .equilibrium import dew_p

        return dew_p(self, T=T, y=y, options=options)

    def dew_t(self, P, y, *, options=None):
        """Solve a neutral dew temperature from vapor composition at fixed pressure."""
        from .equilibrium import dew_t

        return dew_t(self, P=P, y=y, options=options)

    def lle_tp(self, T, P, z, *, options=None, initial_phases=None):
        """Validate a neutral liquid-liquid TP flash request."""
        from .equilibrium import lle_flash

        return lle_flash(self, T=T, P=P, z=z, options=options, initial_phases=initial_phases)

    def stability_tp(self, T, P, z, *, options=None, parent_phase=None, trial_phases=None):
        """Validate a neutral TP stability request."""
        from .equilibrium import neutral_stability

        return neutral_stability(
            self,
            T=T,
            P=P,
            z=z,
            options=options,
            parent_phase=parent_phase,
            trial_phases=trial_phases,
        )

    def electrolyte_stability_tp(self, T, P, *, z=None, solvent_feed=None, salt_molality=None, options=None):
        """Validate a transformed-basis electrolyte TP stability request."""
        from .equilibrium import electrolyte_stability

        return electrolyte_stability(
            self,
            T=T,
            P=P,
            z=z,
            solvent_feed=solvent_feed,
            salt_molality=salt_molality,
            options=options,
        )

    def electrolyte_lle_tp(
        self,
        T,
        P,
        *,
        z=None,
        solvent_feed=None,
        salt_molality=None,
        initial_phases=None,
        options=None,
    ):
        """Validate a charge-constrained electrolyte LLE request."""
        from .equilibrium import electrolyte_lle_flash_native

        return electrolyte_lle_flash_native(
            self,
            T=T,
            P=P,
            z=z,
            solvent_feed=solvent_feed,
            salt_molality=salt_molality,
            initial_phases=initial_phases,
            options=options,
        )

    def electrolyte_bubble_p(
        self,
        T,
        *,
        x_liq=None,
        z=None,
        vapor_species=None,
        volatile_species=None,
        nonvolatile_species=None,
        options=None,
    ):
        """Solve fixed-liquid electrolyte bubble pressure."""
        from .electrolyte_bubble import electrolyte_bubble_pressure

        return electrolyte_bubble_pressure(
            self,
            T=T,
            x_liq=x_liq,
            z=z,
            vapor_species=vapor_species,
            volatile_species=volatile_species,
            nonvolatile_species=nonvolatile_species,
            options=options,
        )

    def chemical_equilibrium(
        self,
        T,
        P,
        *,
        balances,
        totals,
        reactions,
        initial_x=None,
        z=None,
        options=None,
    ):
        """Validate homogeneous reactive speciation and use the accepted Ipopt route."""
        from .reactive_speciation import ReactiveSpeciationOptions, solve_reactive_speciation

        if initial_x is None:
            if z is None:
                raise InputError("chemical_equilibrium requires z or initial_x as the initial composition.")
            initial_x = z
        reactive_options = options if options is not None else ReactiveSpeciationOptions()
        if not isinstance(reactive_options, ReactiveSpeciationOptions):
            raise InputError("chemical_equilibrium options must be a ReactiveSpeciationOptions instance.")
        return solve_reactive_speciation(
            species=self.species,
            mixture_factory=lambda x, T, P: self,
            T=T,
            P=P,
            balances=balances,
            totals=totals,
            reactions=reactions,
            initial_x=initial_x,
            options=reactive_options,
        )

    def reactive_electrolyte_bubble_p(
        self,
        T,
        *,
        P_seed,
        balances,
        totals,
        reactions,
        initial_x=None,
        z=None,
        vapor_species=None,
        volatile_species=None,
        nonvolatile_species=None,
        options=None,
    ):
        """Run staged reactive speciation followed by fixed-liquid electrolyte bubble pressure."""
        from .reactive_electrolyte import ReactiveElectrolyteBubbleOptions, solve_reactive_electrolyte_bubble

        if initial_x is None:
            if z is None:
                raise InputError(
                    "reactive_electrolyte_bubble_pressure requires z or initial_x as the initial composition."
                )
            initial_x = z
        reactive_bubble_options = options if options is not None else ReactiveElectrolyteBubbleOptions()
        if not isinstance(reactive_bubble_options, ReactiveElectrolyteBubbleOptions):
            raise InputError(
                "reactive_electrolyte_bubble_pressure options must be a ReactiveElectrolyteBubbleOptions instance."
            )
        return solve_reactive_electrolyte_bubble(
            species=self.species,
            mixture_factory=lambda x, T, P: self,
            T=T,
            P_seed=P_seed,
            balances=balances,
            totals=totals,
            reactions=reactions,
            initial_x=initial_x,
            vapor_species=vapor_species if vapor_species is not None else volatile_species,
            volatile_species=volatile_species,
            nonvolatile_species=nonvolatile_species,
            options=reactive_bubble_options,
        )

    def reactive_staged_equilibrium(
        self,
        T,
        P,
        *,
        balances,
        totals,
        reactions,
        initial_x=None,
        z=None,
        phase_kind="auto",
        options=None,
        phase_options=None,
        phase_kwargs=None,
    ):
        """Run staged chemical equilibrium followed by an explicit phase route."""
        from .reactive_speciation import ReactiveSpeciationOptions
        from .reactive_staged import solve_reactive_staged_equilibrium

        if initial_x is None:
            if z is None:
                raise InputError("reactive_staged_equilibrium requires z or initial_x as the initial composition.")
            initial_x = z
        reactive_options = options if options is not None else ReactiveSpeciationOptions()
        if not isinstance(reactive_options, ReactiveSpeciationOptions):
            raise InputError("reactive_staged_equilibrium options must be a ReactiveSpeciationOptions instance.")
        return solve_reactive_staged_equilibrium(
            species=self.species,
            mixture_factory=lambda x, T, P: self,
            T=T,
            P=P,
            balances=balances,
            totals=totals,
            reactions=reactions,
            initial_x=initial_x,
            phase_kind=phase_kind,
            speciation_options=reactive_options,
            phase_options=phase_options,
            phase_kwargs=phase_kwargs,
        )

    def equilibrium(
        self,
        kind="tp_flash",
        T=None,
        P=None,
        z=None,
        x_liq=None,
        volatile_species=None,
        vapor_species=None,
        nonvolatile_species=None,
        solvent_feed=None,
        salt_molality=None,
        balances=None,
        totals=None,
        reactions=None,
        initial_x=None,
        options=None,
        phase_kind=None,
        phase_options=None,
        backend=None,
        initial_phases=None,
        parent_phase=None,
        trial_phases=None,
    ):
        """Run a native-backed equilibrium calculation for this mixture."""
        from .equilibrium_core.classify import classify_equilibrium_route
        from .reactive_speciation import ReactiveSpeciationOptions

        if kind in {"reactive_staged_equilibrium", "reactive_staged"}:
            if balances is None or totals is None or reactions is None:
                raise InputError(f"{kind} requires balances, totals, and reactions.")
            phase_kwargs = {
                "vapor_species": vapor_species,
                "volatile_species": volatile_species,
                "nonvolatile_species": nonvolatile_species,
                "solvent_feed": solvent_feed,
                "salt_molality": salt_molality,
                "initial_phases": initial_phases,
                "parent_phase": parent_phase,
                "trial_phases": trial_phases,
            }
            return self.reactive_staged_equilibrium(
                T=T,
                P=P,
                balances=balances,
                totals=totals,
                reactions=reactions,
                initial_x=initial_x,
                z=z,
                phase_kind=phase_kind if phase_kind is not None else ("auto" if backend is None else str(backend)),
                options=options,
                phase_options=phase_options,
                phase_kwargs=phase_kwargs,
            )

        if kind in {
            "reactive_lle",
            "reactive_lle_flash",
            "reactive_electrolyte_lle",
            "reactive_electrolyte_lle_flash",
        }:
            from .equilibrium import reactive_phase_equilibrium

            if balances is None or totals is None or reactions is None:
                raise InputError(f"{kind} requires balances, totals, and reactions.")
            if kind in {"reactive_lle", "reactive_lle_flash"}:
                if phase_kind is not None and phase_kind not in {"lle_flash", "lle_tp"}:
                    raise InputError("reactive_lle phase_kind must be 'lle_flash' when provided.")
                phase_kind = "lle_flash"
            if kind in {"reactive_electrolyte_lle", "reactive_electrolyte_lle_flash"}:
                if phase_kind is not None and phase_kind not in {"electrolyte_lle", "electrolyte_lle_tp"}:
                    raise InputError("reactive_electrolyte_lle phase_kind must be 'electrolyte_lle' when provided.")
                phase_kind = "electrolyte_lle"
            phase_kwargs = {
                "vapor_species": vapor_species,
                "volatile_species": volatile_species,
                "nonvolatile_species": nonvolatile_species,
                "solvent_feed": solvent_feed,
                "salt_molality": salt_molality,
                "initial_phases": initial_phases,
                "parent_phase": parent_phase,
                "trial_phases": trial_phases,
            }
            return reactive_phase_equilibrium(
                self,
                T=T,
                P=P,
                balances=balances,
                totals=totals,
                reactions=reactions,
                initial_x=initial_x,
                z=z,
                phase_kind=phase_kind,
                options=options,
                phase_options=phase_options,
                phase_kwargs=phase_kwargs,
            )

        if kind in {"reactive_electrolyte_bubble_pressure", "reactive_electrolyte_bubble"}:
            if balances is None or totals is None or reactions is None:
                raise InputError("reactive_electrolyte_bubble_pressure requires balances, totals, and reactions.")
            if initial_x is None:
                if z is None:
                    raise InputError(
                        "reactive_electrolyte_bubble_pressure requires z or initial_x as the initial composition."
                    )
                initial_x = z
            if P is None:
                raise InputError("reactive_electrolyte_bubble_pressure requires P as the pressure seed.")
            return self.reactive_electrolyte_bubble_p(
                T=T,
                P_seed=P,
                balances=balances,
                totals=totals,
                reactions=reactions,
                initial_x=initial_x,
                z=z,
                vapor_species=vapor_species if vapor_species is not None else volatile_species,
                volatile_species=volatile_species,
                nonvolatile_species=nonvolatile_species,
                options=options,
            )

        if kind in {"chemical_equilibrium", "reactive_speciation"}:
            if balances is None or totals is None or reactions is None:
                raise InputError("chemical_equilibrium requires balances, totals, and reactions.")
            if initial_x is None:
                if z is None:
                    raise InputError("chemical_equilibrium requires z or initial_x as the initial composition.")
                initial_x = z
            if backend not in (None, "native"):
                raise InputError("chemical_equilibrium backend must be None or 'native'.")
            return self.chemical_equilibrium(
                T=T,
                P=P,
                balances=balances,
                totals=totals,
                reactions=reactions,
                initial_x=initial_x,
                z=z,
                options=options,
            )

        if kind in {"reactive_stability", "chemical_stability"}:
            if balances is None or totals is None or reactions is None:
                raise InputError("reactive_stability requires balances, totals, and reactions.")
            if initial_x is None:
                if z is None:
                    raise InputError("reactive_stability requires z or initial_x as the initial composition.")
                initial_x = z
            if backend not in (None, "native"):
                raise InputError("reactive_stability backend must be None or 'native'.")
            if initial_phases is not None:
                raise InputError("initial_phases is not supported for kind='reactive_stability'.")
            if (
                x_liq is not None
                or volatile_species is not None
                or vapor_species is not None
                or nonvolatile_species is not None
            ):
                raise InputError(
                    "x_liq and vapor species controls are only supported for kind='electrolyte_bubble_pressure'."
                )
            if solvent_feed is not None or salt_molality is not None:
                raise InputError("solvent_feed and salt_molality are not supported for kind='reactive_stability'.")
            reactive_options = options if options is not None else ReactiveSpeciationOptions()
            if not isinstance(reactive_options, ReactiveSpeciationOptions):
                raise InputError("reactive_stability options must be a ReactiveSpeciationOptions instance.")
            chemical = self.chemical_equilibrium(
                T=T,
                P=P,
                balances=balances,
                totals=totals,
                reactions=reactions,
                initial_x=initial_x,
                z=z,
                options=reactive_options,
            )
            equilibrated_z = [chemical.x[label] for label in self.species]
            charges = np.asarray(self.parameters.get("z", []), dtype=float).flatten()
            if charges.size and np.any(np.abs(charges) > 1.0e-12):
                if parent_phase is not None or trial_phases is not None:
                    raise InputError("parent_phase and trial_phases are not supported for ionic reactive_stability.")
                result = self.electrolyte_stability_tp(T=T, P=P, z=equilibrated_z, options=phase_options)
            else:
                result = self.stability_tp(
                    T=T,
                    P=P,
                    z=equilibrated_z,
                    options=phase_options,
                    parent_phase=parent_phase,
                    trial_phases=trial_phases,
                )
            diagnostics = dict(result.diagnostics)
            diagnostics["reactive_phase_method"] = "chemical_equilibrium_then_native_stability"
            diagnostics["reactive_feed_composition"] = equilibrated_z
            diagnostics["reactive_chemical_equilibrium"] = chemical.to_dict()
            diagnostics["phase_equilibrium_handoff"] = chemical.diagnostics.get("phase_equilibrium_handoff")
            return type(result)(
                backend=result.backend,
                problem_kind=result.problem_kind,
                stable=result.stable,
                min_tpd=result.min_tpd,
                parent_phase=result.parent_phase,
                trial_phase=result.trial_phase,
                trial_composition=result.trial_composition,
                trials=result.trials,
                diagnostics=diagnostics,
            )

        if kind == "auto":
            route = classify_equilibrium_route(self, kind, backend)
            if route["route"] == "electrolyte_bubble":
                return self.electrolyte_bubble_p(
                    T=T,
                    x_liq=x_liq,
                    z=z,
                    volatile_species=volatile_species,
                    vapor_species=vapor_species,
                    nonvolatile_species=nonvolatile_species,
                    options=options,
                )
            if route["route"] == "electrolyte_lle":
                try:
                    result = self.electrolyte_lle_tp(
                        T=T,
                        P=P,
                        z=z,
                        solvent_feed=solvent_feed,
                        salt_molality=salt_molality,
                        initial_phases=initial_phases,
                        options=options,
                    )
                except SolutionError as exc:
                    diagnostics = dict(getattr(exc, "diagnostics", None) or {})
                    diagnostics["equilibrium_route"] = route["route"]
                    diagnostics["route_reason"] = route["reason"]
                    raise SolutionError(exc.message, diagnostics) from exc
                diagnostics = dict(result.diagnostics)
                diagnostics["equilibrium_route"] = route["route"]
                diagnostics["route_reason"] = route["reason"]
                return type(result)(
                    backend=result.backend,
                    problem_kind=result.problem_kind,
                    phases=result.phases,
                    stable=result.stable,
                    split_detected=result.split_detected,
                    diagnostics=diagnostics,
                )
            if route["route"] == "neutral_lle":
                return self.lle_tp(T=T, P=P, z=z, options=options, initial_phases=initial_phases)
            return _result_with_route_diagnostics(self.flash_tp(T=T, P=P, z=z, options=options), route)
        if kind in {
            "bubble_p",
            "neutral_bubble_p",
            "bubble_t",
            "neutral_bubble_t",
            "dew_p",
            "neutral_dew_p",
            "dew_t",
            "neutral_dew_t",
        }:
            route = classify_equilibrium_route(self, kind, backend)
            if route["route"] != "neutral_vle":
                raise InputError(f"kind='{kind}' is not a neutral vapor-liquid route.")
            if solvent_feed is not None or salt_molality is not None:
                raise InputError("solvent_feed and salt_molality are only supported for kind='electrolyte_lle'.")
            if volatile_species is not None or vapor_species is not None or nonvolatile_species is not None:
                raise InputError("vapor species controls are only supported for kind='electrolyte_bubble_pressure'.")
            if initial_phases is not None:
                raise InputError("initial_phases is only supported for kind='lle_flash'.")
            if parent_phase is not None or trial_phases is not None:
                raise InputError("parent_phase and trial_phases are only supported for kind='stability'.")
            if backend not in (None, "native", "neutral_vle"):
                raise InputError("Neutral bubble/dew backend must be None, 'native', or 'neutral_vle'.")
            if kind in {"bubble_p", "neutral_bubble_p"}:
                if P is not None:
                    raise InputError("P is solved by kind='bubble_p'.")
                composition = x_liq if x_liq is not None else z
                if composition is None:
                    raise InputError("kind='bubble_p' requires x_liq or z as the liquid composition.")
                return _result_with_route_diagnostics(self.bubble_p(T=T, x=composition, options=options), route)
            if kind in {"bubble_t", "neutral_bubble_t"}:
                if T is not None:
                    raise InputError("T is solved by kind='bubble_t'.")
                composition = x_liq if x_liq is not None else z
                if composition is None:
                    raise InputError("kind='bubble_t' requires x_liq or z as the liquid composition.")
                return _result_with_route_diagnostics(self.bubble_t(P=P, x=composition, options=options), route)
            if x_liq is not None:
                raise InputError("x_liq is only supported for bubble-point routes.")
            if z is None:
                raise InputError(f"kind='{kind}' requires z as the vapor composition.")
            if kind in {"dew_p", "neutral_dew_p"}:
                if P is not None:
                    raise InputError("P is solved by kind='dew_p'.")
                return _result_with_route_diagnostics(self.dew_p(T=T, y=z, options=options), route)
            if T is not None:
                raise InputError("T is solved by kind='dew_t'.")
            return _result_with_route_diagnostics(self.dew_t(P=P, y=z, options=options), route)
        if kind in {"electrolyte_bubble_pressure", "electrolyte_bubble", "bubble_pressure"}:
            if P is not None:
                raise InputError("P is solved by kind='electrolyte_bubble_pressure'.")
            if solvent_feed is not None or salt_molality is not None:
                raise InputError(
                    "solvent_feed and salt_molality are not supported for kind='electrolyte_bubble_pressure'."
                )
            if initial_phases is not None:
                raise InputError("initial_phases is not supported for kind='electrolyte_bubble_pressure'.")
            if parent_phase is not None or trial_phases is not None:
                raise InputError("parent_phase and trial_phases are only supported for kind='stability'.")
            if backend not in (None, "native"):
                raise InputError("Electrolyte bubble backend must be None or 'native'.")
            return self.electrolyte_bubble_p(
                T=T,
                x_liq=x_liq,
                z=z,
                volatile_species=volatile_species,
                vapor_species=vapor_species,
                nonvolatile_species=nonvolatile_species,
                options=options,
            )
        if kind == "tp_flash":
            if solvent_feed is not None or salt_molality is not None:
                raise InputError("solvent_feed and salt_molality are only supported for kind='electrolyte_lle'.")
            if (
                x_liq is not None
                or volatile_species is not None
                or vapor_species is not None
                or nonvolatile_species is not None
            ):
                raise InputError(
                    "x_liq and vapor species controls are only supported for kind='electrolyte_bubble_pressure'."
                )
            if initial_phases is not None:
                raise InputError("initial_phases is only supported for kind='lle_flash'.")
            if parent_phase is not None or trial_phases is not None:
                raise InputError("parent_phase and trial_phases are only supported for kind='stability'.")
            if backend not in (None, "native", "neutral_vle"):
                raise InputError("TP flash backend must be None, 'native', or 'neutral_vle'.")
            route = classify_equilibrium_route(self, kind, backend)
            return _result_with_route_diagnostics(self.flash_tp(T=T, P=P, z=z, options=options), route)
        if kind == "lle_flash":
            if solvent_feed is not None or salt_molality is not None:
                raise InputError("solvent_feed and salt_molality are only supported for kind='electrolyte_lle'.")
            if (
                x_liq is not None
                or volatile_species is not None
                or vapor_species is not None
                or nonvolatile_species is not None
            ):
                raise InputError(
                    "x_liq and vapor species controls are only supported for kind='electrolyte_bubble_pressure'."
                )
            if parent_phase is not None or trial_phases is not None:
                raise InputError("parent_phase and trial_phases are only supported for kind='stability'.")
            if backend not in (None, "native", "neutral_lle"):
                raise InputError("LLE flash backend must be None, 'native', or 'neutral_lle'.")
            return self.lle_tp(T=T, P=P, z=z, options=options, initial_phases=initial_phases)
        if kind in {"electrolyte_lle", "electrolyte_lle_flash"}:
            if (
                x_liq is not None
                or volatile_species is not None
                or vapor_species is not None
                or nonvolatile_species is not None
            ):
                raise InputError(
                    "x_liq and vapor species controls are only supported for kind='electrolyte_bubble_pressure'."
                )
            if parent_phase is not None or trial_phases is not None:
                raise InputError("parent_phase and trial_phases are only supported for kind='stability'.")
            if backend not in (None, "native", "electrolyte_lle"):
                raise InputError("Electrolyte LLE backend must be None, 'native', or 'electrolyte_lle'.")
            return self.electrolyte_lle_tp(
                T=T,
                P=P,
                z=z,
                solvent_feed=solvent_feed,
                salt_molality=salt_molality,
                initial_phases=initial_phases,
                options=options,
            )
        if kind == "electrolyte_stability":
            if (
                x_liq is not None
                or volatile_species is not None
                or vapor_species is not None
                or nonvolatile_species is not None
            ):
                raise InputError(
                    "x_liq and vapor species controls are only supported for kind='electrolyte_bubble_pressure'."
                )
            if initial_phases is not None:
                raise InputError("initial_phases is not supported for kind='electrolyte_stability'.")
            if parent_phase is not None or trial_phases is not None:
                raise InputError("parent_phase and trial_phases are not supported for kind='electrolyte_stability'.")
            if backend not in (None, "native"):
                raise InputError("Electrolyte stability backend must be None or 'native'.")
            return self.electrolyte_stability_tp(
                T=T,
                P=P,
                z=z,
                solvent_feed=solvent_feed,
                salt_molality=salt_molality,
                options=options,
            )
        if kind == "stability":
            if solvent_feed is not None or salt_molality is not None:
                raise InputError("solvent_feed and salt_molality are only supported for kind='electrolyte_lle'.")
            if (
                x_liq is not None
                or volatile_species is not None
                or vapor_species is not None
                or nonvolatile_species is not None
            ):
                raise InputError(
                    "x_liq and vapor species controls are only supported for kind='electrolyte_bubble_pressure'."
                )
            if initial_phases is not None:
                raise InputError("initial_phases is only supported for kind='lle_flash'.")
            if backend not in (None, "native"):
                raise InputError("Stability backend must be None or 'native'.")
            return self.stability_tp(
                T=T,
                P=P,
                z=z,
                options=options,
                parent_phase=parent_phase,
                trial_phases=trial_phases,
            )
        raise InputError(
            "Only kind='tp_flash', kind='auto', kind='bubble_p', kind='bubble_t', kind='dew_p', kind='dew_t', kind='lle_flash', kind='electrolyte_lle', kind='electrolyte_bubble_pressure', kind='electrolyte_stability', kind='stability', kind='chemical_equilibrium', kind='reactive_staged_equilibrium', kind='reactive_lle', kind='reactive_electrolyte_lle', kind='reactive_stability', or kind='reactive_electrolyte_bubble_pressure' is supported by equilibrium."
        )

    def equilibrium_curve(self, points, *, kind="electrolyte_lle", T=None, P=None, options=None, initial_phases=None):
        """Solve an ordered equilibrium curve with route-owned canonical initial points."""

        if initial_phases is not None:
            raise InputError("equilibrium_curve uses route-owned canonical initial points; initial_phases is not accepted.")
        results = []
        for point in points:
            payload = dict(point)
            if "z" not in payload and "solvent_feed" not in payload:
                raise InputError("equilibrium_curve points must provide z or solvent_feed.")
            result = self.equilibrium(
                kind=kind,
                T=payload.pop("T", T),
                P=payload.pop("P", P),
                options=payload.pop("options", options),
                **payload,
            )
            results.append(result)
        return results

    def equilibrium_sweep(
        self,
        points,
        *,
        kind="reactive_electrolyte_bubble_pressure",
        balances=None,
        reactions=None,
        vapor_species=None,
        volatile_species=None,
        nonvolatile_species=None,
        options=None,
    ):
        """Solve an ordered reactive-electrolyte bubble sweep."""
        if kind not in {"reactive_electrolyte_bubble_pressure", "reactive_electrolyte_bubble"}:
            raise InputError("equilibrium_sweep currently supports kind='reactive_electrolyte_bubble_pressure'.")
        if balances is None or reactions is None:
            raise InputError("reactive electrolyte equilibrium_sweep requires balances and reactions.")
        if vapor_species is None and volatile_species is None:
            raise InputError("reactive electrolyte equilibrium_sweep requires vapor_species or volatile_species.")
        from .reactive_electrolyte import solve_reactive_electrolyte_bubble_sweep

        return solve_reactive_electrolyte_bubble_sweep(
            species=self.species,
            mixture_factory=lambda x, T, P: self,
            points=points,
            balances=balances,
            reactions=reactions,
            vapor_species=vapor_species if vapor_species is not None else volatile_species,
            volatile_species=volatile_species,
            nonvolatile_species=nonvolatile_species,
            options=options,
        )

    def __repr__(self):
        """Return a short debugging representation of the mixture."""
        return f"ePCSAFTMixture(ncomp={self.ncomp}, species={self._species})"


class ePCSAFTState:
    """Immutable thermodynamic state bound to one mixture."""

    def __init__(self, mixture, T, x, P=None, rho=None, phase="liq", rho_guess=None):
        """Create a state with exactly one intensive variable fixed.

        Pressure-based states solve the internal T, P, x -> rho closure eagerly.
        """
        if not isinstance(mixture, ePCSAFTMixture):
            raise InputError("mixture must be a ePCSAFTMixture instance.")
        mix = mixture
        x, _params = ensure_numpy_input(x, mix._params)
        x = np.asarray(x, dtype=float).flatten()
        ncomp = int(mix.ncomp)
        if x.size != ncomp:
            raise InputError(f"State composition length ({int(x.size)}) must match mixture component count ({ncomp}).")
        # ensure_numpy_input may normalize a scalar mixture parameter path, but the
        # state should retain the original mixture data unchanged.
        phase_num = phase_to_int(phase)
        has_p = P is not None
        has_rho = rho is not None
        if has_p == has_rho:
            raise InputError("Provide exactly one of P or rho when constructing a state.")
        has_rho_guess = rho_guess is not None
        if has_rho_guess and not has_p:
            raise InputError("rho_guess is only supported for pressure-based states constructed with P.")
        if has_rho_guess and (not np.isfinite(float(rho_guess)) or float(rho_guess) <= 0.0):
            raise InputError("rho_guess must be finite and positive.")
        if has_p:
            check_input(x, {"temperature": T, "pressure": P})
        else:
            check_input(x, {"temperature": T, "density": rho})
        cpp_x = np_to_vector_double(x)
        try:
            self._native = _core.NativeState(
                mix._native,
                float(T),
                cpp_x,
                phase_num,
                has_p,
                float(P) if P is not None else 0.0,
                has_rho,
                float(rho) if rho is not None else 0.0,
                has_rho_guess,
                float(rho_guess) if rho_guess is not None else 0.0,
            )
        except _NATIVE_CALL_ERRORS as exc:
            variable_name = "P" if has_p else "rho"
            mode_name = "pressure" if has_p else "density"
            variable_value = P if has_p else rho
            message = _state_construction_error_message(
                T, x, phase, ncomp, mode_name, variable_name, variable_value, exc
            )
            diagnostics = None
            if has_p and hasattr(mix._native, "last_density_diagnostics"):
                diagnostics = dict(mix._native.last_density_diagnostics())
                for context in diagnostics.get("density_failure_contexts", []):
                    context["phase_label"] = "state"
                    context["phase_kind"] = "vap" if phase_num == 1 else "liq"
                if diagnostics.get("density_scan_summary"):
                    diagnostics["density_scan_summary"]["phase_label"] = "state"
                    diagnostics["density_scan_summary"]["phase_kind"] = "vap" if phase_num == 1 else "liq"
            raise SolutionError(message, diagnostics) from exc
        self._mixture = mixture
        self._x = np.asarray(x, dtype=float)
        self._T = float(T)
        self._P = None if P is None else float(P)
        self._rho = None if rho is None else float(rho)
        self._phase = phase_num

    @property
    def mixture(self):
        """Return the parent mixture."""
        return self._mixture

    @property
    def T(self):
        """Return the state temperature in kelvin."""
        return self._T

    @property
    def x(self):
        """Return the state composition as a NumPy array."""
        return np.asarray(self._x, dtype=float)

    @property
    def phase(self):
        """Return the native liquid/vapor phase flag."""
        return self._phase

    def method_aliases(self):
        """Return the canonical full-name -> abbreviation map for state methods."""
        return dict(STATE_METHOD_ALIAS_MAP)

    def pressure(self):
        """Return the pressure of the bound state."""
        return float(self._native.pressure())

    def density(self, units="molar"):
        """Return the state density in the requested units."""
        token = str(units or "molar").strip().lower()
        if token in {"molar", "mol", "molar_density", "mol/m^3", "mol/m3"}:
            return float(self._native.density())
        if token in {"mass", "mass_density", "kg/m^3", "kg/m3"}:
            return self.mass_density()
        raise InputError("density units must be 'molar'/'mol/m^3' or 'mass'/'kg/m^3'.")

    def molar_density(self):
        """Return the molar density of the bound state in mol/m^3."""
        return float(self._native.density())

    def mass_density(self):
        """Return the mass density of the bound state in kg/m^3."""
        mix = self._mixture
        mw = np.asarray(mix._params.get("MW", []), dtype=float).flatten()
        if mw.size == 0:
            raise InputError("Mass density requires component molecular weights in kg/mol.")
        if mw.size != self._x.size:
            raise InputError("Mass density requires one molecular-weight value per component.")
        return float(self.molar_density() * float(np.dot(np.asarray(self._x, dtype=float), mw)))

    def _native_args_copy(self):
        return create_struct(self._mixture._params)

    def _neutral_binary_kij_property_derivatives(self):
        if int(self._x.size) != 2:
            return None
        try:
            raw = _core._native_cppad_neutral_binary_pair_properties(
                self._T,
                self.molar_density(),
                np_to_vector_double(self._x),
                self._native_args_copy(),
            )
        except _NATIVE_CALL_ERRORS:
            return None
        return raw if bool(raw.get("supported", raw.get("cppad_used", False))) else None

    def _neutral_binary_pair_parameter_order(self, raw):
        species = tuple(self._mixture.species)
        order = []
        for name in tuple(raw.get("parameter_names", ("k_ij",))):
            order.append(f"{name}:{species[0]}:{species[1]}")
        return tuple(order)

    @staticmethod
    def _neutral_binary_pair_jacobian(raw, quantity, nrows):
        columns = []
        for name in tuple(raw.get("parameter_names", ("k_ij",))):
            columns.append(np.asarray(raw[f"{name}_{quantity}_derivative"], dtype=float).reshape((nrows, 1)))
        return np.concatenate(columns, axis=1) if columns else np.empty((nrows, 0), dtype=float)

    def _pure_neutral_parameter_derivatives(self):
        if int(self._x.size) != 1:
            return None
        try:
            raw = _core._native_cppad_pure_neutral_parameters(
                self._T,
                self.molar_density(),
                self._native_args_copy(),
            )
        except _NATIVE_CALL_ERRORS:
            return None
        return raw if bool(raw.get("supported", raw.get("cppad_used", False))) else None

    def _pure_neutral_parameter_order(self):
        species = str(self._mixture.species[0])
        return (f"m:{species}", f"sigma:{species}", f"epsilon:{species}")

    def pressure_density_derivative_result(self):
        """Return the pressure derivative with respect to density in result-contract form."""
        ncomp = int(self._x.size)
        try:
            raw = _core._native_cppad_pressure_density(
                self._T,
                self.molar_density(),
                np_to_vector_double(self._x),
                self._native_args_copy(),
            )
        except _NATIVE_CALL_ERRORS as exc:
            _unsupported_derivative(f"pressure-density derivative backend failed: {exc}")
        if not bool(raw.get("cppad_compiled", False)):
            _unsupported_derivative("CppAD derivative backend did not report required compiled support.")
        return _derivative_result_payload(
            supported=bool(raw.get("cppad_used", False)),
            backend=str(raw.get("derivative_backend", "cppad")),
            message=str(raw.get("message", "CppAD pressure-density derivative available")),
            value=[self.pressure()],
            jacobian=np.asarray(raw.get("jacobian_row_major", []), dtype=float),
            shape=[1, 1],
            output_order=("pressure",),
            variable_order=("density",),
            component_order=tuple(self._mixture.species[:ncomp]),
        )

    def pressure_composition_derivative_result(self):
        """Return pressure composition derivatives when production support exists."""
        _unsupported_derivative("pressure composition derivatives require EOS pressure sensitivity routing.")

    def pressure_parameter_derivative_result(self):
        """Return pressure parameter-derivative support status."""
        raw = self._pure_neutral_parameter_derivatives()
        if raw is not None:
            jacobian = np.asarray(raw["jacobian_row_major"], dtype=float).reshape((3, 3))
            return _derivative_result_payload(
                supported=True,
                backend=str(raw.get("derivative_backend", "cppad")),
                message=str(raw.get("message", "CppAD pure-neutral m/sigma/epsilon pressure derivatives available")),
                value=[float(np.asarray(raw["value"], dtype=float)[0])],
                jacobian=jacobian[0:1, :],
                shape=[1, 3],
                output_order=("pressure",),
                parameter_order=self._pure_neutral_parameter_order(),
                source_equation_ids=("pressure_from_z", "ares_hc", "ares_disp"),
            )
        raw = self._neutral_binary_kij_property_derivatives()
        if raw is not None:
            parameter_order = self._neutral_binary_pair_parameter_order(raw)
            return _derivative_result_payload(
                supported=True,
                backend=str(raw.get("backend", "cppad")),
                message=str(raw.get("message", "CppAD neutral binary pair-parameter pressure derivative available")),
                value=[float(raw.get("k_ij_pressure", self.pressure()))],
                jacobian=self._neutral_binary_pair_jacobian(raw, "pressure", 1),
                shape=[1, len(parameter_order)],
                output_order=("pressure",),
                parameter_order=parameter_order,
                source_equation_ids=("pressure_from_z", "ares_disp"),
            )
        _unsupported_derivative("pressure parameter derivatives require an analytic or CppAD route for this state.")

    def density_pressure_derivative_result(self):
        """Return density pressure derivatives when production support exists."""
        _unsupported_derivative("density-root implicit pressure derivatives require an analytic or CppAD implicit route.")

    def ares_composition_derivative_result(self):
        """Return residual-Helmholtz composition derivatives in the public result shape."""
        result = self._composition_derivative_residual_helmholtz_result()
        backend, message = _backend_from_contribution_details(
            result["derivative_backend"],
            available=bool(result["derivative_available"]),
        )
        return _derivative_result_payload(
            supported=backend != "unsupported",
            backend=backend,
            message=message,
            value=[self.residual_helmholtz()],
            jacobian=np.asarray(result["total"], dtype=float).reshape((1, int(self._x.size))),
            shape=[1, int(self._x.size)],
            output_order=("ares_residual",),
            variable_order=tuple(self._mixture.species),
            backend_details=dict(result["derivative_backend"]),
            source_equation_ids=("ares_total",),
        )

    def chemical_potential_composition_derivative_result(self):
        """Return chemical-potential composition derivatives when production support exists."""
        _unsupported_derivative("chemical-potential composition Jacobians require second composition derivatives.")

    def chemical_potential_parameter_derivative_result(self):
        """Return residual-chemical-potential parameter derivatives where production support exists."""
        ncomp = int(self._x.size)
        value = self.residual_chemical_potential()
        raw = self._pure_neutral_parameter_derivatives()
        if raw is not None:
            jacobian = np.asarray(raw["jacobian_row_major"], dtype=float).reshape((3, 3))
            return _derivative_result_payload(
                supported=True,
                backend=str(raw.get("derivative_backend", "cppad")),
                message=str(
                    raw.get(
                        "message",
                        "CppAD pure-neutral m/sigma/epsilon residual-chemical-potential derivatives available",
                    )
                ),
                value=[float(np.asarray(raw["value"], dtype=float)[1])],
                jacobian=jacobian[1:2, :],
                shape=[1, 3],
                component_order=tuple(self._mixture.species),
                parameter_order=self._pure_neutral_parameter_order(),
                value_basis="residual_chemical_potential",
                source_equation_ids=("mu_res", "ares_hc", "ares_disp"),
            )
        raw = self._neutral_binary_kij_property_derivatives()
        if raw is not None:
            species = tuple(self._mixture.species)
            parameter_order = self._neutral_binary_pair_parameter_order(raw)
            return _derivative_result_payload(
                supported=True,
                backend=str(raw.get("backend", "cppad")),
                message=str(
                    raw.get(
                        "message",
                        "CppAD neutral binary pair-parameter residual-chemical-potential derivative available",
                    )
                ),
                value=np.asarray(raw.get("k_ij_residual_chemical_potential", value), dtype=float),
                jacobian=self._neutral_binary_pair_jacobian(raw, "residual_chemical_potential", ncomp),
                shape=[ncomp, len(parameter_order)],
                component_order=species,
                parameter_order=parameter_order,
                value_basis="residual_chemical_potential",
                source_equation_ids=("mu_res", "ares_disp"),
            )
        _unsupported_derivative(
            "chemical-potential parameter derivatives require an analytic or CppAD route for this state."
        )

    def ln_fugacity_composition_derivative_result(self):
        """Return ln-fugacity composition derivatives for pressure-based native states."""
        ncomp = int(self._x.size)
        if self._P is None:
            _unsupported_derivative(
                "ln-fugacity composition Jacobians require second composition derivatives through a pressure-based phase state."
            )
        try:
            raw = _core._native_phase_state_ln_fugacity_composition_sensitivity(
                self._T,
                self._P,
                np_to_vector_double(self._x),
                self._phase,
                self._native_args_copy(),
            )
        except _NATIVE_CALL_ERRORS as exc:
            _unsupported_derivative(f"phase-state ln-fugacity composition sensitivity backend failed: {exc}")
        if not bool(raw.get("supported", False)):
            _unsupported_derivative(
                str(
                    raw.get(
                        "message",
                        "phase-state ln-fugacity composition sensitivities require a CppAD implicit route.",
                    )
                )
            )
        shape = [int(raw["shape"][0]), int(raw["shape"][1])]
        if shape != [ncomp, ncomp]:
            _unsupported_derivative("invalid phase-state sensitivity payload.")
        return _derivative_result_payload(
            supported=True,
            backend=str(raw.get("backend", "cppad_implicit")),
            message=str(raw.get("message", "CppAD phase-state ln-fugacity composition sensitivities available.")),
            value=np.asarray(raw.get("ln_fugacity", self.fugacity_coefficient(natural_log=True)), dtype=float),
            jacobian=np.asarray(raw.get("jacobian_row_major", []), dtype=float),
            shape=shape,
            component_order=tuple(self._mixture.species),
            variable_order=tuple(self._mixture.species),
            value_basis="natural_log_fugacity_coefficient_at_fixed_pressure",
            density=float(raw.get("density", self.molar_density())),
            density_backend=str(raw.get("density_backend", "implicit_density_root")),
            density_composition_derivative=np.asarray(raw.get("density_composition_derivative", []), dtype=float),
            pressure_density_derivative=float(raw.get("pressure_density_derivative", 0.0)),
            pressure_composition_fixed_density_derivative=np.asarray(
                raw.get("pressure_composition_fixed_density_derivative", []),
                dtype=float,
            ),
            ln_fugacity_density_derivative=np.asarray(raw.get("ln_fugacity_density_derivative", []), dtype=float),
            fixed_density_jacobian=np.asarray(raw.get("fixed_density_jacobian_row_major", []), dtype=float).reshape(
                (ncomp, ncomp)
            ),
            source_equation_ids=("lnphi_total", "pressure_from_z", "density_root_implicit"),
        )

    def ln_fugacity_parameter_derivative_result(self):
        """Return ln-fugacity parameter derivatives where production support exists."""
        ncomp = int(self._x.size)
        value = self.fugacity_coefficient(natural_log=True)
        raw = self._pure_neutral_parameter_derivatives()
        if raw is not None:
            jacobian = np.asarray(raw["jacobian_row_major"], dtype=float).reshape((3, 3))
            return _derivative_result_payload(
                supported=True,
                backend=str(raw.get("derivative_backend", "cppad")),
                message=str(raw.get("message", "CppAD pure-neutral m/sigma/epsilon ln-fugacity derivatives available")),
                value=[float(np.asarray(raw["value"], dtype=float)[2])],
                jacobian=jacobian[2:3, :],
                shape=[1, 3],
                component_order=tuple(self._mixture.species),
                parameter_order=self._pure_neutral_parameter_order(),
                value_basis="natural_log_fugacity_coefficient",
                source_equation_ids=("lnphi_total", "ares_hc", "ares_disp"),
            )
        raw = self._neutral_binary_kij_property_derivatives()
        if raw is not None:
            species = tuple(self._mixture.species)
            parameter_order = self._neutral_binary_pair_parameter_order(raw)
            return _derivative_result_payload(
                supported=True,
                backend=str(raw.get("backend", "cppad")),
                message=str(raw.get("message", "CppAD neutral binary pair-parameter ln-fugacity derivative available")),
                value=np.asarray(raw.get("k_ij_ln_fugacity", value), dtype=float),
                jacobian=self._neutral_binary_pair_jacobian(raw, "ln_fugacity", ncomp),
                shape=[ncomp, len(parameter_order)],
                component_order=species,
                parameter_order=parameter_order,
                value_basis="natural_log_fugacity_coefficient",
                source_equation_ids=("lnphi_total", "ares_disp"),
            )
        born = self.born_ssmds_liquid_derivatives()
        if not bool(born.get("supported", False)):
            _unsupported_derivative(
                str(born.get("message", "ln-fugacity parameter derivatives require an analytic or CppAD route."))
            )
        jacobian = np.concatenate(
            [
                np.asarray(born["lnfug_d_d_born"], dtype=float),
                np.asarray(born["lnfug_d_f_solv"], dtype=float),
            ],
            axis=1,
        )
        parameter_order = tuple(f"d_born:{name}" for name in self._mixture.species) + tuple(
            f"f_solv:{name}" for name in self._mixture.species
        )
        return _derivative_result_payload(
            supported=True,
            backend=str(born.get("backend", "analytic")),
            message=str(born.get("message", "analytic liquid-electrolyte SSM+DS Born parameter derivatives available")),
            value=value,
            jacobian=jacobian,
            shape=[ncomp, 2 * ncomp],
            component_order=tuple(self._mixture.species),
            parameter_order=parameter_order,
            phase_scope="liquid_electrolyte_only",
        )

    def activity_composition_derivative_result(self, species=None):
        """Return activity-coefficient composition derivatives when production support exists."""
        _ = self._mixture.species if species is None else [str(s) for s in species]
        _unsupported_derivative("activity-coefficient composition Jacobians require an analytic or CppAD route.")

    def activity_parameter_derivative_result(self, species=None):
        """Return activity-coefficient parameter derivatives where production support exists."""
        species = self._mixture.species if species is None else [str(s) for s in species]
        ncomp = int(self._x.size)
        born = self.born_ssmds_liquid_derivatives()
        try:
            gamma = self.activity_coefficient(species=species)
            value = np.log(np.asarray([gamma[label] for label in species], dtype=float))
        except _DERIVATIVE_VALUE_ERRORS:
            value = np.asarray([], dtype=float)
        if not bool(born.get("supported", False)):
            _unsupported_derivative(
                str(born.get("message", "activity parameter derivatives require an analytic or CppAD route."))
            )
        jacobian = np.concatenate(
            [
                np.asarray(born["lngamma_d_d_born"], dtype=float),
                np.asarray(born["lngamma_d_f_solv"], dtype=float),
            ],
            axis=1,
        )
        parameter_order = tuple(f"d_born:{name}" for name in self._mixture.species) + tuple(
            f"f_solv:{name}" for name in self._mixture.species
        )
        return _derivative_result_payload(
            supported=True,
            backend=str(born.get("backend", "analytic")),
            message=str(born.get("message", "analytic liquid-electrolyte SSM+DS Born activity derivatives available")),
            value=value,
            jacobian=jacobian,
            shape=[ncomp, 2 * ncomp],
            component_order=tuple(species),
            parameter_order=parameter_order,
            value_basis="natural_log_activity_coefficient",
            phase_scope="liquid_electrolyte_only",
        )

    def relative_permittivity_composition_derivative_result(self):
        """Return relative-permittivity composition derivatives in result-contract form."""
        try:
            epsilon, deps_dx = self.relative_permittivity()
        except _DERIVATIVE_VALUE_ERRORS as exc:
            _unsupported_derivative(f"relative-permittivity derivative backend failed: {exc}")
        return _derivative_result_payload(
            supported=True,
            backend=_relative_permittivity_backend(self._mixture._params),
            message="relative-permittivity composition derivative available",
            value=[float(epsilon)],
            jacobian=np.asarray(deps_dx, dtype=float).reshape((1, int(self._x.size))),
            shape=[1, int(self._x.size)],
            output_order=("relative_permittivity",),
            variable_order=tuple(self._mixture.species),
            source_equation_ids=("relative_permittivity",),
        )

    def relative_permittivity_parameter_derivative_result(self):
        """Return relative-permittivity parameter derivatives for supported dielectric rules."""
        params = self._mixture._params
        try:
            epsilon, _ = self.relative_permittivity()
        except _DERIVATIVE_VALUE_ERRORS as exc:
            _unsupported_derivative(f"relative-permittivity derivative backend failed: {exc}")
        rule, _mode = _state_rel_perm_rule_and_mode(params)
        dielc = np.asarray(params.get("dielc", []), dtype=float).flatten()
        if rule != 1 or dielc.size != self._x.size:
            _unsupported_derivative(
                "relative-permittivity parameter derivatives are implemented for linear mole-fraction mixing only."
            )
        return _derivative_result_payload(
            supported=True,
            backend="analytic",
            message="analytic linear-mixing relative-permittivity parameter derivatives available",
            value=[float(epsilon)],
            jacobian=np.asarray(self._x, dtype=float).reshape((1, int(self._x.size))),
            shape=[1, int(self._x.size)],
            output_order=("relative_permittivity",),
            parameter_order=tuple(f"relative_permittivity:{name}" for name in self._mixture.species),
            source_equation_ids=("relative_permittivity",),
        )

    def derivative_coverage_matrix(self):
        """Return structured implemented derivative backend coverage for EOS property workflows."""
        ncomp = int(self._x.size)
        rows = []

        def classify(*, supported, out_of_scope, override=None):
            if override:
                return str(override)
            if out_of_scope:
                return "out_of_scope"
            if supported:
                return "production_supported"
            raise SolutionError("Derivative coverage rows must be supported or explicitly out of scope.")

        def add(quantity, derivative, result, *, out_of_scope=False, source_equation_ids=()):
            supported = bool(result.get("supported", False))
            row_out_of_scope = bool(out_of_scope)
            if not supported and not row_out_of_scope:
                return
            rows.append(
                {
                    "quantity": quantity,
                    "derivative": derivative,
                    "backend": "out_of_scope" if row_out_of_scope else str(result.get("backend", "unspecified")),
                    "supported": supported,
                    "classification": classify(
                        supported=supported,
                        out_of_scope=row_out_of_scope,
                        override=result.get("classification"),
                    ),
                    "source_equation_ids": list(source_equation_ids),
                    "parameter_family": str(result.get("parameter_family", "")),
                }
            )

        def add_result_factory(quantity, derivative, factory, *, source_equation_ids=()):
            try:
                result = factory()
            except _NATIVE_CALL_ERRORS:
                return
            add(quantity, derivative, result, source_equation_ids=source_equation_ids)

        composition = self._composition_derivative_residual_helmholtz_result()
        details = dict(composition["derivative_backend"])
        assoc_active = bool(np.asarray(self._mixture._params.get("assoc_num", []), dtype=int).size)
        for key, quantity, eqid in (
            ("hc", "hard_chain", "ares_hc"),
            ("disp", "dispersion", "ares_disp"),
            ("assoc", "association", "ares_assoc"),
            ("ion", "debye_huckel / ion", "ares_dh"),
            ("born", "born_direct", "ares_born"),
        ):
            backend = details.get(key, "unsupported")
            if key == "assoc" and assoc_active and backend == "analytic":
                backend = "analytic_implicit"
            inactive_term = key in {"assoc", "ion", "born"} and np.allclose(composition["terms"][key], 0.0)
            add(
                quantity,
                "composition",
                {"supported": backend != "unsupported", "backend": backend, "message": ""},
                out_of_scope=inactive_term,
                source_equation_ids=(eqid,),
            )

        add_result_factory(
            "born_ssmds_liquid",
            "d_born/f_solv",
            self.born_ssmds_liquid_derivatives,
            source_equation_ids=("ares_born",),
        )
        add_result_factory("relative_permittivity", "composition", self.relative_permittivity_composition_derivative_result)
        add_result_factory("relative_permittivity", "parameter", self.relative_permittivity_parameter_derivative_result)
        add_result_factory(
            "pressure",
            "density",
            self.pressure_density_derivative_result,
            source_equation_ids=("pressure_from_z",),
        )
        add_result_factory(
            "pressure",
            "parameter",
            self.pressure_parameter_derivative_result,
            source_equation_ids=("pressure_from_z",),
        )
        add_result_factory(
            "fugacity",
            "composition",
            self.ln_fugacity_composition_derivative_result,
            source_equation_ids=("lnphi_total",),
        )
        add_result_factory(
            "fugacity",
            "parameter",
            self.ln_fugacity_parameter_derivative_result,
            source_equation_ids=("lnphi_total",),
        )
        add_result_factory(
            "activity",
            "composition",
            self.activity_composition_derivative_result,
            source_equation_ids=("lngamma_sym",),
        )
        add_result_factory(
            "activity",
            "parameter",
            self.activity_parameter_derivative_result,
            source_equation_ids=("lngamma_sym",),
        )
        add_result_factory(
            "chemical_potential",
            "composition",
            self.chemical_potential_composition_derivative_result,
            source_equation_ids=("mu_res",),
        )
        add_result_factory(
            "chemical_potential",
            "parameter",
            self.chemical_potential_parameter_derivative_result,
            source_equation_ids=("mu_res",),
        )
        add_result_factory(
            "density_root",
            "pressure",
            self.density_pressure_derivative_result,
            source_equation_ids=("density_root",),
        )
        if ncomp == 0:
            raise SolutionError("Derivative coverage matrix requires at least one component.")
        return rows

    def _temperature_derivative_residual_helmholtz_term_result(self):
        result = self._native.temperature_derivative_residual_helmholtz_result()
        return _scalar_terms_dict(result)

    def _composition_derivative_residual_helmholtz_result(self):
        ncomp = int(self._x.size)
        result = self._native.composition_derivative_residual_helmholtz_result()
        dadx = _vector_terms_dict(result.dadx, ncomp, "dadx")
        ares = _scalar_terms_dict(result.ares)
        sum_x = _scalar_terms_dict(result.sum_x_dadx)
        z_raw = _scalar_terms_dict(result.z_raw)
        z_terms = _scalar_terms_dict(result.z)
        terms = {name: np.asarray(dadx[name], dtype=float) for name in _CONTRIBUTION_NAMES}
        return {
            "total": _sum_vector_terms(terms),
            "terms": terms,
            "ares_terms": {name: float(ares[name]) for name in _CONTRIBUTION_NAMES},
            "sum_x_terms": {name: float(sum_x[name]) for name in _CONTRIBUTION_NAMES},
            "z_raw_terms": {name: float(z_raw[name]) for name in _CONTRIBUTION_NAMES},
            "z_terms": {name: float(z_terms[name]) for name in _CONTRIBUTION_NAMES},
            "z_total": float(z_terms["total"]),
            "derivative_backend": {str(k): str(v) for k, v in dict(result.derivative_backend).items()},
            "derivative_available": bool(result.derivative_available),
        }

    def _fugacity_coefficient_term_result(self):
        ncomp = int(self._x.size)
        result = self._native.fugacity_coefficient_result()
        mu = _vector_terms_dict(result.mu, ncomp, "mu")
        lnfug = _vector_terms_dict(result.lnfugcoef, ncomp, "lnfugcoef")
        dadx = _vector_terms_dict(result.composition.dadx, ncomp, "dadx")
        ares = _scalar_terms_dict(result.composition.ares)
        sum_x = _scalar_terms_dict(result.composition.sum_x_dadx)
        z_raw = _scalar_terms_dict(result.composition.z_raw)
        z_terms = _scalar_terms_dict(result.composition.z)
        return {
            "mu_hc": mu["hc"],
            "mu_disp": mu["disp"],
            "mu_assoc": mu["assoc"],
            "mu_ion": mu["ion"],
            "mu_born": mu["born"],
            "mu_total": mu["total"],
            "lnfugcoef_total": lnfug["total"],
            "lnfugcoef_hc": lnfug["hc"],
            "lnfugcoef_disp": lnfug["disp"],
            "lnfugcoef_assoc": lnfug["assoc"],
            "lnfugcoef_ion": lnfug["ion"],
            "lnfugcoef_born": lnfug["born"],
            "dadx_hc": dadx["hc"],
            "dadx_disp": dadx["disp"],
            "dadx_assoc": dadx["assoc"],
            "dadx_ion": dadx["ion"],
            "dadx_born": dadx["born"],
            "a_hc": ares["hc"],
            "a_disp": ares["disp"],
            "a_assoc": ares["assoc"],
            "a_ion": ares["ion"],
            "a_born": ares["born"],
            "sum_x_dadx_hc": sum_x["hc"],
            "sum_x_dadx_disp": sum_x["disp"],
            "sum_x_dadx_assoc": sum_x["assoc"],
            "sum_x_dadx_ion": sum_x["ion"],
            "sum_x_dadx_born": sum_x["born"],
            "z_raw_hc": z_raw["hc"],
            "z_raw_disp": z_raw["disp"],
            "z_raw_assoc": z_raw["assoc"],
            "z_raw_ion": z_raw["ion"],
            "z_raw_born": z_raw["born"],
            "z_hc": z_terms["hc"],
            "z_disp": z_terms["disp"],
            "z_assoc": z_terms["assoc"],
            "z_ion": z_terms["ion"],
            "z_born": z_terms["born"],
            "z_total": z_terms["total"],
        }

    def compressibility_factor(self, return_contribution_terms=False):
        """Return the compressibility factor."""
        if not return_contribution_terms:
            return float(self._native.compressibility_factor())
        result = self._native.compressibility_factor_result()
        payload_dict = _scalar_terms_dict(result.terms)
        terms = {name: float(payload_dict[name]) for name in _CONTRIBUTION_NAMES}
        terms["ideal"] = 1.0
        return {
            "total": float(payload_dict["total"]),
            "terms": terms,
        }

    def pressure_contributions(self):
        """Return pressure contributions derived from native compressibility-factor terms."""
        z_payload = self.compressibility_factor(return_contribution_terms=True)
        pressure_scale = float(self.molar_density() * _GAS_CONSTANT * self._T)
        return {
            "total": float(self.pressure()),
            "terms": {
                name: float(value) * pressure_scale
                for name, value in _public_contribution_terms(z_payload["terms"]).items()
            },
            "term_basis": "pressure_from_compressibility_factor",
            "compressibility_factor": z_payload,
        }

    def residual_helmholtz(self, return_contribution_terms=False):
        """Return the residual Helmholtz energy."""
        if not return_contribution_terms:
            return float(self._native.residual_helmholtz())
        result = self._native.residual_helmholtz_result()
        payload_dict = _scalar_terms_dict(result)
        terms = {name: float(payload_dict[name]) for name in _CONTRIBUTION_NAMES}
        return {
            "total": float(payload_dict["total"]),
            "terms": terms,
        }

    def helmholtz_contributions(self):
        """Return the native residual Helmholtz contribution map."""
        payload = self.residual_helmholtz(return_contribution_terms=True)
        return {
            **payload,
            "terms": _public_contribution_terms(payload["terms"]),
            "term_basis": "dimensionless_residual_helmholtz",
            "ideal": {"available": False, "reason": "ideal Helmholtz decomposition is not exposed by the native API"},
        }

    def residual_helmholtz_contributions(self):
        """Return the native residual Helmholtz contribution map."""
        payload = self.residual_helmholtz(return_contribution_terms=True)
        return {
            **payload,
            "terms": _public_contribution_terms(payload["terms"]),
            "term_basis": "dimensionless_residual_helmholtz",
        }

    def temperature_derivative_residual_helmholtz(self, return_contribution_terms=False):
        """Return the temperature derivative of the residual Helmholtz energy."""
        if not return_contribution_terms:
            return float(self._native.temperature_derivative_residual_helmholtz())
        result = self._temperature_derivative_residual_helmholtz_term_result()
        terms = {name: float(result[name]) for name in _CONTRIBUTION_NAMES}
        return {
            "total": float(result["total"]),
            "terms": terms,
        }

    def composition_derivative_residual_helmholtz(self):
        """Return the composition-derivative contribution breakdown."""
        return self._composition_derivative_residual_helmholtz_result()

    def residual_enthalpy(self):
        """Return the residual enthalpy."""
        return float(self._native.residual_enthalpy())

    def residual_entropy(self):
        """Return the residual entropy."""
        return float(self._native.residual_entropy())

    def residual_gibbs(self):
        """Return the residual Gibbs energy."""
        return float(self._native.residual_gibbs())

    def residual_chemical_potential(self, return_contribution_terms=False):
        """Return the residual chemical potentials."""
        if not return_contribution_terms:
            return vector_to_array(self._native.residual_chemical_potential())
        result = self._native.residual_chemical_potential_result()
        payload_dict = _vector_terms_dict(result.mu, int(self._x.size), "mu")
        terms = {name: np.asarray(payload_dict[name], dtype=float) for name in _CONTRIBUTION_NAMES}
        return {
            "total": np.asarray(payload_dict["total"], dtype=float),
            "terms": terms,
        }

    def chemical_potential_contributions(self):
        """Return the native residual chemical-potential contribution map."""
        payload = self.residual_chemical_potential(return_contribution_terms=True)
        return {
            **payload,
            "terms": _public_contribution_terms(payload["terms"]),
            "term_basis": "residual_chemical_potential",
        }

    def _activity_coefficient_bundle(self, species=None, solvent=None, include_aux=False):
        """Build the full activity-coefficient payload for internal reuse."""
        species = self._mixture.species if species is None else [str(s) for s in species]
        if len(species) != self._x.size:
            raise InputError(f"species length ({len(species)}) must match composition length ({self._x.size}).")
        has_solvent_override, solvent_index = _resolve_solvent_override(self._mixture, species, solvent)
        include_aux_c = bool(include_aux)
        has_solvent_override_c = bool(has_solvent_override)
        solvent_index_c = int(solvent_index)
        out = self._native.activity_coefficient_native(include_aux_c, has_solvent_override_c, solvent_index_c)
        pair_cat = np.asarray(out.pair_cation_indices, dtype=int)
        pair_an = np.asarray(out.pair_anion_indices, dtype=int)
        pair_labels = tuple(species[int(ic)] + species[int(ia)] for ic, ia in zip(pair_cat.tolist(), pair_an.tolist()))
        ion_idx = np.sort(
            np.unique(
                np.concatenate([np.asarray(out.cation_indices, dtype=int), np.asarray(out.anion_indices, dtype=int)])
            )
        )
        ion_labels = tuple(species[int(i)] for i in ion_idx.tolist())
        return ActivityCoefficientResult(
            species=tuple(species),
            component_activity_coefficients=np.asarray(out.component_activity_coefficients, dtype=float),
            solvation_free_energy_values=np.asarray(out.solvation_free_energy, dtype=float),
            mean_ionic_activity_coefficients_mole_fraction_values=np.asarray(
                out.mean_ionic_activity_coefficients_mole_fraction, dtype=float
            ),
            mean_ionic_activity_coefficients_molality_values=np.asarray(
                out.mean_ionic_activity_coefficients_molality, dtype=float
            ),
            pair_labels=pair_labels,
            ion_labels=ion_labels,
            ion_indices=ion_idx,
            cation_indices=np.asarray(out.cation_indices, dtype=int),
            anion_indices=np.asarray(out.anion_indices, dtype=int),
            solvent_indices=np.asarray(out.solvent_indices, dtype=int),
            pair_cation_indices=pair_cat,
            pair_anion_indices=pair_an,
            pair_nu_cation=np.asarray(out.pair_nu_cation, dtype=int),
            pair_nu_anion=np.asarray(out.pair_nu_anion, dtype=int),
            pair_molality=np.asarray(out.pair_molality, dtype=float),
            pair_conversion_factor=np.asarray(out.pair_conversion_factor, dtype=float),
            solvent_index=int(out.solvent_index),
            osmotic_coefficient=float(out.osmotic_coefficient),
        )

    def activity_coefficient(self, species=None, solvent=None, mean_ionic_form=False, basis="mole"):
        """Return activity coefficients in the requested form."""
        species = self._mixture.species if species is None else [str(s) for s in species]
        if len(species) != self._x.size:
            raise InputError(f"species length ({len(species)}) must match composition length ({self._x.size}).")
        has_solvent_override, solvent_index = _resolve_solvent_override(self._mixture, species, solvent)
        include_aux_c = False
        has_solvent_override_c = bool(has_solvent_override)
        solvent_index_c = int(solvent_index)
        out = self._native.activity_coefficient_native(include_aux_c, has_solvent_override_c, solvent_index_c)
        if mean_ionic_form:
            token = str(basis).strip().lower()
            if token in {"mole", "mole_fraction", "molefraction", "x"}:
                values = np.asarray(out.mean_ionic_activity_coefficients_mole_fraction, dtype=float)
            elif token in {"molality", "m"}:
                values = np.asarray(out.mean_ionic_activity_coefficients_molality, dtype=float)
            else:
                raise InputError("basis must be one of: 'mole', 'mole_fraction', 'x', 'molality', 'm'.")
            pair_cat = np.asarray(out.pair_cation_indices, dtype=int)
            pair_an = np.asarray(out.pair_anion_indices, dtype=int)
            return {
                species[int(ic)] + species[int(ia)]: float(value)
                for ic, ia, value in zip(pair_cat.tolist(), pair_an.tolist(), values.tolist())
            }
        values = np.asarray(out.component_activity_coefficients, dtype=float)
        return {label: float(value) for label, value in zip(species, values.tolist())}

    def activity_coefficient_contributions(self, *args, **kwargs):
        """Reject unsupported additive activity-coefficient decomposition explicitly."""
        raise InputError(
            "Activity coefficients are available through activity_coefficient(), but additive "
            "hard_chain/dispersion/association/ionic/born activity-coefficient contributions are not exposed by "
            "the native activity API."
        )

    def fugacity_coefficient(self, natural_log=True, return_contribution_terms=False):
        """Return fugacity coefficients, defaulting to natural-log form."""
        ln_total = vector_to_array(self._native.ln_fugacity_coefficient())
        if not return_contribution_terms:
            if natural_log:
                return ln_total
            return np.exp(ln_total)
        result = self._fugacity_coefficient_term_result()
        ln_term_total = np.asarray(result["lnfugcoef_total"], dtype=float)
        terms = {name: np.asarray(result["lnfugcoef_" + name], dtype=float) for name in _CONTRIBUTION_NAMES}
        out = {
            "total": ln_total if natural_log else np.exp(ln_total),
            "terms": terms,
            "term_basis": "natural_log",
            "terms_total_natural_log": ln_term_total,
        }
        return out

    def ln_fugacity_coefficient_contributions(self):
        """Return native natural-log fugacity-coefficient contribution terms."""
        payload = self.fugacity_coefficient(natural_log=True, return_contribution_terms=True)
        return {
            **payload,
            "terms": _public_contribution_terms(payload["terms"]),
        }

    def born_ssmds_liquid_derivatives(self):
        """Return analytic liquid-electrolyte SSM+DS Born parameter derivatives."""
        if self._phase != 0:
            raise InputError("unsupported: SSM+DS Born derivatives are liquid-electrolyte only.")
        try:
            payload = dict(self._native.born_ssmds_liquid_derivatives())
        except _NATIVE_CALL_ERRORS as exc:
            _unsupported_derivative(str(exc))
        ncomp = int(payload.get("ncomp", self._x.size))
        for key in (
            "a_born_d_d_born",
            "a_born_d_f_solv",
            "mu_res_d_d_born_row_major",
            "mu_res_d_f_solv_row_major",
            "lnfug_d_d_born_row_major",
            "lnfug_d_f_solv_row_major",
            "lngamma_d_d_born_row_major",
            "lngamma_d_f_solv_row_major",
        ):
            payload[key] = np.asarray(payload[key], dtype=float)
        for key in (
            "mu_res_d_d_born",
            "mu_res_d_f_solv",
            "lnfug_d_d_born",
            "lnfug_d_f_solv",
            "lngamma_d_d_born",
            "lngamma_d_f_solv",
        ):
            row_major_key = key + "_row_major"
            payload[key] = np.asarray(payload[row_major_key], dtype=float).reshape((ncomp, ncomp))
        payload["parameter_order"] = tuple(self._mixture.species)
        payload["component_order"] = tuple(self._mixture.species)
        payload["phase_scope"] = "liquid_electrolyte_only"
        payload["parameters"] = ("d_born", "f_solv")
        payload["vapor_support"] = False
        return payload

    def relative_permittivity(self):
        """Return the dielectric model evaluation for the current state."""
        flat = vector_to_array(self._native.relative_permittivity())
        return flat[0], flat[1:]

    def osmotic_coefficient(self):
        """Return the osmotic coefficient."""
        return np.asarray([self._native.osmotic_coefficient()], dtype=float)

    def state_diagnostics(self, species=None):
        """Return a diagnostic dictionary of the main state properties."""
        mix = self._mixture
        species = self._mixture.species if species is None else species
        z = np.asarray(mix._params.get("z", []), dtype=float).flatten()
        has_ions = bool(np.any(np.abs(z) > 1e-12))
        terms = self._fugacity_coefficient_term_result()
        fugacity_coefficient = self.fugacity_coefficient(natural_log=False)
        if has_ions:
            activity_coefficient = self.activity_coefficient(species=species, mean_ionic_form=False)
            relative_permittivity = self.relative_permittivity()
            osmotic_coefficient = self.osmotic_coefficient()
            mean_ionic_activity_coefficient_molality = self.activity_coefficient(
                species=species,
                mean_ionic_form=True,
                basis="molality",
            )
            mean_ionic_activity_coefficient_mole = self.activity_coefficient(
                species=species,
                mean_ionic_form=True,
                basis="mole",
            )
            solvation_free_energy = self.solvation_free_energy(species=species)
        else:
            activity_coefficient = {}
            relative_permittivity = None
            osmotic_coefficient = None
            mean_ionic_activity_coefficient_molality = {}
            mean_ionic_activity_coefficient_mole = {}
            solvation_free_energy = {}
        return {
            "T": self._T,
            "phase": self._phase,
            "x": np.asarray(self._x, dtype=float),
            "pressure": self.pressure(),
            "density": self.molar_density(),
            "density_molar": self.molar_density(),
            "mass_density": self.mass_density() if np.asarray(mix._params.get("MW", []), dtype=float).size else None,
            "compressibility_factor": self.compressibility_factor(),
            "residual_helmholtz": self.residual_helmholtz(),
            "residual_enthalpy": self.residual_enthalpy(),
            "residual_entropy": self.residual_entropy(),
            "residual_gibbs": self.residual_gibbs(),
            "residual_chemical_potential": self.residual_chemical_potential(),
            "fugacity_coefficient": fugacity_coefficient,
            "activity_coefficient": activity_coefficient,
            "fugacity_coefficient_terms": terms,
            "relative_permittivity": relative_permittivity,
            "osmotic_coefficient": osmotic_coefficient,
            "mean_ionic_activity_coefficient_molality": mean_ionic_activity_coefficient_molality,
            "mean_ionic_activity_coefficient_mole": mean_ionic_activity_coefficient_mole,
            "solvation_free_energy": solvation_free_energy,
        }

    def solvation_free_energy(self, species=None):
        """Return ion solvation free-energy values keyed by species."""
        return self._activity_coefficient_bundle(species=species, include_aux=True).solvation_free_energy()

    def __getattr__(self, name):
        """Resolve supported short scientific aliases for state methods."""
        target = _STATE_METHOD_ALIAS_LOOKUP.get(str(name))
        if target is not None:
            return getattr(self, target)
        raise AttributeError(f"{type(self).__name__} has no attribute '{name}'")

    def __repr__(self):
        """Return a short debugging representation of the state."""
        return f"ePCSAFTState(T={self._T}, phase={self._phase}, x={self._x})"


def check_input(x, vars):
    if abs(np.sum(x) - 1) > 1e-7:
        raise InputError(f"The mole fractions do not sum to 1. x = {x}")
    if "temperature" in vars:
        if vars["temperature"] <= 0:
            raise InputError(
                "The {} must be a positive number. {} = {}".format("temperature", "temperature", vars["temperature"])
            )
    if "density" in vars:
        if vars["density"] <= 0:
            raise InputError("The {} must be a positive number. {} = {}".format("density", "density", vars["density"]))
    if "pressure" in vars:
        if vars["pressure"] <= 0:
            raise InputError(
                "The {} must be a positive number. {} = {}".format("pressure", "pressure", vars["pressure"])
            )
    if "Q" in vars:
        if (vars["Q"] < 0) or (vars["Q"] > 1):
            raise InputError("{} must be <= 1 and >= 0. {} = {}".format("Q", "Q", vars["Q"]))


def check_association(params):
    params = deepcopy(params)
    if ("e_assoc" in params) and ("vol_a" not in params):
        raise InputError("e_assoc was given, but not vol_a.")
    elif ("vol_a" in params) and ("e_assoc" not in params):
        raise InputError("vol_a was given, but not e_assoc.")

    if ("e_assoc" in params) and ("assoc_scheme" not in params):
        params["assoc_scheme"] = []
        for a in params["vol_a"]:
            if a != 0:
                params["assoc_scheme"].append("2b")
            else:
                params["assoc_scheme"].append(None)

    if "e_assoc" in params:
        params = create_assoc_matrix(params)

    return params


def create_assoc_matrix(params):
    charge = (
        []
    )  # whether the association site has a partial positive charge (i.e. hydrogen), negative charge, or elements of both (e.g. for acids modelled as type 1)

    scheme_charges = {
        "1": [0],
        "2a": [0, 0],
        "2b": [-1, 1],
        "3a": [0, 0, 0],
        "3b": [-1, -1, 1],
        "4a": [0, 0, 0, 0],
        "4b": [1, 1, 1, -1],
        "4c": [-1, -1, 1, 1],
    }

    assoc_num = []
    for comp in params["assoc_scheme"]:
        if comp is None:
            assoc_num.append(0)
        elif type(comp) is list:
            num = 0
            for site in comp:
                if site.lower() not in scheme_charges:
                    raise InputError(f"{site} is not a valid association type.")
                charge.extend(scheme_charges[site.lower()])
                num += len(scheme_charges[site.lower()])
            assoc_num.append(num)
        else:
            if comp.lower() not in scheme_charges:
                raise InputError(f"{comp} is not a valid association type.")
            charge.extend(scheme_charges[comp.lower()])
            assoc_num.append(len(scheme_charges[comp.lower()]))
    params["assoc_num"] = np.asarray(assoc_num)

    params["assoc_matrix"] = np.zeros(len(charge) * len(charge))
    ctr = 0
    for c1 in charge:
        for c2 in charge:
            if c1 == 0 or c2 == 0:
                params["assoc_matrix"][ctr] = 1
            elif c1 == 1 and c2 == -1:
                params["assoc_matrix"][ctr] = 1
            elif c1 == -1 and c2 == 1:
                params["assoc_matrix"][ctr] = 1
            else:
                params["assoc_matrix"][ctr] = 0
            ctr += 1

    return params


def ensure_numpy_input(x, params):
    if np.isscalar(x):
        x = np.asarray([x], dtype=float)
    if np.isscalar(params["m"]):
        params["m"] = np.asarray([params["m"]], dtype=float)
    if np.isscalar(params["s"]):
        params["s"] = np.asarray([params["s"]], dtype=float)
    if np.isscalar(params["e"]):
        params["e"] = np.asarray([params["e"]], dtype=float)
    return x, params


def _safe_unit_log(values, floor=1e-300):
    vals = np.asarray(values, dtype=float)
    return np.log(np.maximum(vals, floor))


def _resolve_solvent_override(mixture, species, solvent):
    z = np.asarray(mixture._params.get("z", []), dtype=float).flatten()
    if z.size == 0:
        if solvent is not None:
            raise InputError("solvent override requires ionic parameters with params['z'].")
        return False, -1
    neutral_idx = np.where(np.abs(z) <= 1e-12)[0]
    if neutral_idx.size == 0:
        if solvent is not None:
            raise InputError("solvent override requires at least one neutral solvent species.")
        return False, -1
    if solvent is None:
        return False, -1
    if isinstance(solvent, (int, np.integer)):
        idx = int(solvent)
    else:
        token = str(solvent)
        if token not in species:
            raise InputError(f"Unknown solvent label '{token}'. Available species={list(species)}")
        idx = species.index(token)
    if idx < 0 or idx >= z.size:
        raise InputError(f"solvent index out of bounds: {idx} (ncomp={int(z.size)})")
    if abs(float(z[idx])) > 1e-12:
        raise InputError("solvent override must reference a neutral species (z=0).")
    return True, idx


def np_to_vector_double(np_array):
    """Take a numpy array and return a C++ vector."""
    return np.asarray(np_array, dtype=float).ravel().tolist()


def np_to_vector_int(np_array):
    """Take a numpy array and return a C++ vector."""
    return np.asarray(np_array, dtype=int).ravel().tolist()


def create_struct(params):
    """Convert ePC-SAFT parameters to a C++ struct."""

    cppargs = _core.NativeArgs()

    for removed_key in ("dipm", "dip_num"):
        if removed_key in params:
            raise ValueError(
                f'Removed polar parameter "{removed_key}" is not supported by the active ePC-SAFT package.'
            )

    cppargs.mixed_rel_perm_water_index = -1
    cppargs.m = np_to_vector_double(params["m"])
    ncomp = len(np.asarray(params["m"]).flatten())
    cppargs.s = np_to_vector_double(params["s"])
    cppargs.e = np_to_vector_double(params["e"])
    if "k_ij" in params:
        cppargs.k_ij = np_to_vector_double(params["k_ij"])
    if ("e_assoc" in params) and np.any(params["e_assoc"]):
        cppargs.e_assoc = np_to_vector_double(params["e_assoc"])
    if ("vol_a" in params) and np.any(params["vol_a"]):
        cppargs.vol_a = np_to_vector_double(params["vol_a"])
    z_arr = None
    if "z" in params:
        z_arr = np.asarray(params["z"], dtype=float).flatten()
        if z_arr.size not in (0, ncomp):
            raise ValueError(f'params["z"] must have length {ncomp} (or be empty), got {z_arr.size}.')
        if z_arr.size == ncomp:
            cppargs.z = np_to_vector_double(z_arr)
    if "dielc" in params:
        dielc_arr = np.asarray(params["dielc"], dtype=float).flatten()
        if dielc_arr.size != ncomp:
            raise ValueError(f'params["dielc"] must have length {ncomp}, got {dielc_arr.size}.')
        cppargs.dielc = np_to_vector_double(dielc_arr)
    if "MW" in params:
        mw_arr = np.asarray(params["MW"], dtype=float).flatten()
        if mw_arr.size != ncomp:
            raise ValueError(f'params["MW"] must have length {ncomp}, got {mw_arr.size}.')
        cppargs.mw = np_to_vector_double(mw_arr)
    if "mixed_rel_perm_a" in params:
        mixed_a_arr = np.asarray(params["mixed_rel_perm_a"], dtype=float).flatten()
        if mixed_a_arr.size != ncomp:
            raise ValueError(f'params["mixed_rel_perm_a"] must have length {ncomp}, got {mixed_a_arr.size}.')
        cppargs.mixed_rel_perm_a = np_to_vector_double(mixed_a_arr)
    if "mixed_rel_perm_b" in params:
        mixed_b_arr = np.asarray(params["mixed_rel_perm_b"], dtype=float).flatten()
        if mixed_b_arr.size != ncomp:
            raise ValueError(f'params["mixed_rel_perm_b"] must have length {ncomp}, got {mixed_b_arr.size}.')
        cppargs.mixed_rel_perm_b = np_to_vector_double(mixed_b_arr)
    if "mixed_rel_perm_c" in params:
        mixed_c_arr = np.asarray(params["mixed_rel_perm_c"], dtype=float).flatten()
        if mixed_c_arr.size != ncomp:
            raise ValueError(f'params["mixed_rel_perm_c"] must have length {ncomp}, got {mixed_c_arr.size}.')
        cppargs.mixed_rel_perm_c = np_to_vector_double(mixed_c_arr)
    if "mixed_rel_perm_mask" in params:
        mixed_mask_arr = np.asarray(params["mixed_rel_perm_mask"], dtype=int).flatten()
        if mixed_mask_arr.size != ncomp:
            raise ValueError(f'params["mixed_rel_perm_mask"] must have length {ncomp}, got {mixed_mask_arr.size}.')
        cppargs.mixed_rel_perm_mask = np_to_vector_int(mixed_mask_arr)
    if "mixed_rel_perm_water_index" in params:
        cppargs.mixed_rel_perm_water_index = int(params["mixed_rel_perm_water_index"])
    if len(cppargs.z) > 0 and len(cppargs.dielc) == 0:
        raise ValueError('Electrolyte parameters require params["dielc"] as a per-species array.')
    d_born_arr = None
    if "d_born" in params:
        d_born_arr = np.asarray(params["d_born"], dtype=float).flatten()
        if d_born_arr.size != ncomp:
            raise ValueError(f'params["d_born"] must have length {ncomp}, got {d_born_arr.size}.')
        cppargs.d_born = np_to_vector_double(d_born_arr)
    if "f_solv" in params:
        cppargs.f_solv = np_to_vector_double(np.asarray(params["f_solv"], dtype=float))

    unsupported_flat_elec_keys = {
        "dielc_rule",
        "dielc_diff_mode",
        "born_model",
        "born_radius_model",
        "born_diff_mode",
        "born_eps_mode",
        "DH_model",
        "bjeruum_treatment",
        "d_ion_mode",
        "include_born_model",
        "d_Born_mode",
        "born_solvation_shell_model",
        "born_dielectric_saturation",
        "born_bulk_mode",
        "mu_DH_diff_mode",
        "mu_DH_comp_dep_rel_perm",
        "mu_DH_include_sum_term",
        "mu_born_diff_mode",
        "mu_born_comp_dep_rel_perm",
        "mu_born_include_sum_term",
        "mu_born_comp_dep_delta_d",
    }
    if any((k in params) for k in unsupported_flat_elec_keys):
        raise ValueError(
            'Flat electrostatic params are no longer supported; provide nested params["elec_model"] schema.'
        )

    elec_model = params.get("elec_model", None)
    if elec_model is None:
        if len(cppargs.z) > 0:
            # Apply the canonical electrolyte defaults when ionic parameters are present.
            elec_model = {
                "rel_perm": {"rule": 1, "differential_mode": "auto"},
                "hc_model": {"dadx_differential_mode": "auto"},
                "disp_model": {"dadx_differential_mode": "auto"},
                "assoc_model": {"dadx_differential_mode": "auto"},
                "DH_model": {
                    "d_ion_mode": 1,
                    "bjeruum_treatment": False,
                    "mu_DH_model": {
                        "differential_mode": "auto",
                        "comp_dep_rel_perm": True,
                        "include_sum_term": True,
                    },
                },
                "include_born_model": True,
                "born_model": {
                    "d_Born_mode": 0,
                    "solvation_shell_model": False,
                    "dielectric_saturation": False,
                    "bulk_mode": "mix",
                    "mu_born_model": {
                        "differential_mode": "auto",
                        "comp_dep_rel_perm": True,
                        "include_sum_term": True,
                        "comp_dep_delta_d": False,
                    },
                },
            }
        else:
            elec_model = {}
    if not isinstance(elec_model, dict):
        raise ValueError('params["elec_model"] must be a dict when provided.')

    def _reject_unknown_keys(mapping, allowed, label):
        unknown = sorted(set(mapping) - set(allowed))
        if unknown:
            raise ValueError(f"{label} contains unsupported key(s): {unknown}.")

    _reject_unknown_keys(
        elec_model,
        {"rel_perm", "hc_model", "disp_model", "assoc_model", "DH_model", "include_born_model", "born_model"},
        'params["elec_model"]',
    )

    def _as_bool(v):
        if isinstance(v, bool):
            return v
        if isinstance(v, (int, np.integer)):
            return bool(v)
        if isinstance(v, str):
            s = v.strip().lower()
            if s in {"1", "true", "yes", "y", "on"}:
                return True
            if s in {"0", "false", "no", "n", "off"}:
                return False
        raise ValueError(f"Could not coerce value to bool: {v}")

    def _as_int_alias(v, aliases):
        if isinstance(v, (int, np.integer)):
            return int(v)
        if isinstance(v, str):
            s = v.strip().lower()
            if s in aliases:
                return int(aliases[s])
            if s.isdigit() or (s.startswith("-") and s[1:].isdigit()):
                return int(s)
        raise ValueError(f"Unknown option value: {v}")

    rule_alias = {
        "constant": 0,
        "rule0": 0,
        "linear": 1,
        "linear-molefraction": 1,
        "linear-mixing-mole": 1,
        "rule1": 1,
        "linear-massfraction": 2,
        "linear-mixing-weight": 2,
        "rule2": 2,
        "combined": 3,
        "rule3": 3,
        "empirical": 4,
        "rule4": 4,
        "rule5": 5,
        "rule6": 6,
        "aqueous-organic": 8,
        "aqueous_organic": 8,
        "mixed-aqueous-organic": 8,
        "mixed_aqueous_organic": 8,
        "rule8": 8,
        "salt-free-massfraction": 9,
        "salt_free_massfraction": 9,
        "salt-free-solvent-massfraction": 9,
        "salt_free_solvent_massfraction": 9,
        "salt-free-solvent-weight": 9,
        "salt_free_solvent_weight": 9,
        "rule9": 9,
    }
    diff_alias = {
        "auto": 3,
        "automatic": 3,
        "analytic": 0,
        "analytical": 0,
        "cppad": 2,
    }
    d_ion_alias = {"t_indep": 0, "t_dep_1": 1, "t_dep_2": 2}
    d_born_alias = {"t_indep": 0, "t_dep_1": 1, "t_dep_2": 2, "fitted_param": 3}
    bulk_alias = {"mix": 0, "bulk": 0, "solvent": 1}

    rel_perm = elec_model.get("rel_perm", {})
    if not isinstance(rel_perm, dict):
        raise ValueError('params["elec_model"]["rel_perm"] must be a dict.')
    _reject_unknown_keys(rel_perm, {"rule", "differential_mode"}, 'params["elec_model"]["rel_perm"]')
    hc_model_dict = elec_model.get("hc_model", {})
    if not isinstance(hc_model_dict, dict):
        raise ValueError('params["elec_model"]["hc_model"] must be a dict.')
    _reject_unknown_keys(hc_model_dict, {"dadx_differential_mode"}, 'params["elec_model"]["hc_model"]')
    disp_model_dict = elec_model.get("disp_model", {})
    if not isinstance(disp_model_dict, dict):
        raise ValueError('params["elec_model"]["disp_model"] must be a dict.')
    _reject_unknown_keys(disp_model_dict, {"dadx_differential_mode"}, 'params["elec_model"]["disp_model"]')
    assoc_model_dict = elec_model.get("assoc_model", {})
    if not isinstance(assoc_model_dict, dict):
        raise ValueError('params["elec_model"]["assoc_model"] must be a dict.')
    _reject_unknown_keys(assoc_model_dict, {"dadx_differential_mode"}, 'params["elec_model"]["assoc_model"]')
    dh_model_dict = elec_model.get("DH_model", {})
    if not isinstance(dh_model_dict, dict):
        raise ValueError('params["elec_model"]["DH_model"] must be a dict.')
    _reject_unknown_keys(
        dh_model_dict, {"d_ion_mode", "bjeruum_treatment", "mu_DH_model"}, 'params["elec_model"]["DH_model"]'
    )
    born_model_dict = elec_model.get("born_model", {})
    if not isinstance(born_model_dict, dict):
        raise ValueError('params["elec_model"]["born_model"] must be a dict.')
    _reject_unknown_keys(
        born_model_dict,
        {"d_Born_mode", "solvation_shell_model", "dielectric_saturation", "bulk_mode", "mu_born_model"},
        'params["elec_model"]["born_model"]',
    )
    mu_dh = dh_model_dict.get("mu_DH_model", {})
    if not isinstance(mu_dh, dict):
        raise ValueError('params["elec_model"]["DH_model"]["mu_DH_model"] must be a dict.')
    _reject_unknown_keys(
        mu_dh,
        {"differential_mode", "comp_dep_rel_perm", "include_sum_term"},
        'params["elec_model"]["DH_model"]["mu_DH_model"]',
    )
    mu_born = born_model_dict.get("mu_born_model", {})
    if not isinstance(mu_born, dict):
        raise ValueError('params["elec_model"]["born_model"]["mu_born_model"] must be a dict.')
    _reject_unknown_keys(
        mu_born,
        {"differential_mode", "comp_dep_rel_perm", "include_sum_term", "comp_dep_delta_d"},
        'params["elec_model"]["born_model"]["mu_born_model"]',
    )

    cppargs.dielc_rule = _as_int_alias(rel_perm.get("rule", 1), rule_alias)
    cppargs.dielc_diff_mode = _as_int_alias(rel_perm.get("differential_mode", "auto"), diff_alias)
    if cppargs.dielc_diff_mode not in (0, 2, 3):
        raise ValueError("Unknown rel_perm differential_mode. Supported values are analytic/cppad/auto (0/2/3).")
    cppargs.hc_dadx_diff_mode = _as_int_alias(hc_model_dict.get("dadx_differential_mode", "auto"), diff_alias)
    if cppargs.hc_dadx_diff_mode not in (0, 2, 3):
        raise ValueError(
            "Unknown hc_model dadx_differential_mode. Supported values are analytic/cppad/auto (0/2/3)."
        )
    cppargs.disp_dadx_diff_mode = _as_int_alias(disp_model_dict.get("dadx_differential_mode", "auto"), diff_alias)
    if cppargs.disp_dadx_diff_mode not in (0, 2, 3):
        raise ValueError(
            "Unknown disp_model dadx_differential_mode. Supported values are analytic/cppad/auto (0/2/3)."
        )
    cppargs.assoc_dadx_diff_mode = _as_int_alias(assoc_model_dict.get("dadx_differential_mode", "auto"), diff_alias)
    if cppargs.assoc_dadx_diff_mode not in (0, 2, 3):
        raise ValueError(
            "Unknown assoc_model dadx_differential_mode. Supported values are analytic/cppad/auto (0/2/3)."
        )
    if cppargs.dielc_rule < 0 or cppargs.dielc_rule > 9:
        raise ValueError("Unknown rel_perm rule. Supported values are 0..9.")

    cppargs.d_ion_mode = _as_int_alias(dh_model_dict.get("d_ion_mode", 1), d_ion_alias)
    if cppargs.d_ion_mode not in (0, 1, 2):
        raise ValueError("Unknown d_ion_mode. Supported values are 0,1,2.")
    bjeruum = _as_bool(dh_model_dict.get("bjeruum_treatment", False))
    cppargs.mu_DH_diff_mode = _as_int_alias(mu_dh.get("differential_mode", "auto"), diff_alias)
    if cppargs.mu_DH_diff_mode not in (0, 2, 3):
        raise ValueError("Unknown mu_DH differential_mode. Supported values are analytic/cppad/auto (0/2/3).")
    cppargs.mu_DH_comp_dep_rel_perm = int(_as_bool(mu_dh.get("comp_dep_rel_perm", True)))
    cppargs.mu_DH_include_sum_term = int(_as_bool(mu_dh.get("include_sum_term", True)))

    cppargs.include_born_model = int(_as_bool(elec_model.get("include_born_model", True)))
    cppargs.d_born_mode = _as_int_alias(born_model_dict.get("d_Born_mode", 0), d_born_alias)
    if cppargs.d_born_mode not in (0, 1, 2, 3):
        raise ValueError("Unknown d_Born_mode. Supported values are 0,1,2,3.")
    cppargs.born_solvation_shell_model = int(_as_bool(born_model_dict.get("solvation_shell_model", False)))
    cppargs.born_dielectric_saturation = int(_as_bool(born_model_dict.get("dielectric_saturation", False)))
    cppargs.born_bulk_mode = _as_int_alias(born_model_dict.get("bulk_mode", "mix"), bulk_alias)
    cppargs.mu_born_diff_mode = _as_int_alias(mu_born.get("differential_mode", "auto"), diff_alias)
    if cppargs.mu_born_diff_mode not in (0, 2, 3):
        raise ValueError("Unknown mu_born differential_mode. Supported values are analytic/cppad/auto (0/2/3).")
    cppargs.mu_born_comp_dep_rel_perm = int(_as_bool(mu_born.get("comp_dep_rel_perm", True)))
    cppargs.mu_born_include_sum_term = int(_as_bool(mu_born.get("include_sum_term", True)))
    cppargs.mu_born_comp_dep_delta_d = int(_as_bool(mu_born.get("comp_dep_delta_d", False)))

    if cppargs.include_born_model == 0:
        cppargs.born_model = 0
        cppargs.born_radius_model = 1
        cppargs.born_diff_mode = 0
        cppargs.born_eps_mode = cppargs.born_bulk_mode
    else:
        if cppargs.born_solvation_shell_model or cppargs.born_dielectric_saturation:
            cppargs.born_model = 2
        else:
            cppargs.born_model = 1

        # Project the canonical Born radius mode into the current native runtime encoding.
        if cppargs.d_born_mode == 0:
            cppargs.born_radius_model = 1
        elif cppargs.d_born_mode == 1:
            cppargs.born_radius_model = 2
        elif cppargs.d_born_mode == 2:
            cppargs.born_radius_model = 3
        else:
            cppargs.born_radius_model = 5

        cppargs.born_eps_mode = cppargs.born_bulk_mode

        if cppargs.mu_born_diff_mode == 1:
            cppargs.born_diff_mode = 1
        elif cppargs.mu_born_diff_mode == 2:
            cppargs.born_diff_mode = 4
        elif cppargs.mu_born_diff_mode == 3:
            cppargs.born_diff_mode = 5
        elif cppargs.mu_born_comp_dep_rel_perm == 0:
            cppargs.born_diff_mode = 3
        elif cppargs.mu_born_include_sum_term == 0:
            cppargs.born_diff_mode = 2
        else:
            cppargs.born_diff_mode = 0

    if cppargs.born_model > 0 and cppargs.born_radius_model in (4, 5):
        if z_arr is None:
            raise ValueError("fitted d_Born_mode requires params['z'] as a per-species array.")
        if d_born_arr is None:
            raise ValueError("fitted d_Born_mode requires params['d_born'] as a per-species array.")
        ion_mask = np.abs(z_arr) > 1e-12
        if np.any(d_born_arr[ion_mask] <= 0.0):
            raise ValueError("fitted d_Born_mode requires positive ionic params['d_born'] values.")

    cppargs.DH_model = 2 if bjeruum else 1
    if cppargs.DH_model == 2:
        raise ValueError("Bjerrum treatment is reserved; DH_model=2 has no active public route.")
    if cppargs.DH_model < 0 or cppargs.DH_model > 2:
        raise ValueError("Unknown DH_model. Supported values are 0, 1, and reserved 2.")
    if "assoc_num" in params:
        cppargs.assoc_num = np_to_vector_int(params["assoc_num"])
    if "assoc_matrix" in params:
        cppargs.assoc_matrix = np_to_vector_int(params["assoc_matrix"])
    if "k_hb" in params:
        cppargs.k_hb = np_to_vector_double(params["k_hb"])
    if "l_ij" in params:
        cppargs.l_ij = np_to_vector_double(params["l_ij"])

    return cppargs


def _fit_pure_neutral_native_debug(
    fixed_payload,
    density_T,
    density_P,
    density_rho_exp,
    density_phase,
    density_scale,
    vle_T,
    vle_P,
    pure_vle_scale,
    x,
):
    params = check_association(dict(fixed_payload))
    cppargs = create_struct(params)
    result = _core._fit_pure_neutral_native_debug(
        cppargs,
        np.asarray(density_T, dtype=float),
        np.asarray(density_P, dtype=float),
        np.asarray(density_rho_exp, dtype=float),
        np.asarray(density_phase, dtype=int),
        float(density_scale),
        np.asarray(vle_T, dtype=float),
        np.asarray(vle_P, dtype=float),
        float(pure_vle_scale),
        np.asarray(x, dtype=float),
    )
    return {
        "objective": float(result["objective"]),
        "gradient": vector_to_array(result["gradient"]),
        "residuals": vector_to_array(result["residuals"]),
        "jacobian_row_major": vector_to_array(result["jacobian_row_major"]),
        "jacobian_shape": (result["jacobian_shape"]),
        "jacobian_available": bool(result.get("jacobian_available", True)),
        "jacobian_backend": str(result.get("jacobian_backend", "unspecified")),
        "density_raw_residuals": vector_to_array(result["density_raw_residuals"]),
        "pure_vle_raw_residuals": vector_to_array(result["pure_vle_raw_residuals"]),
        "residual_evaluations": int(result["residual_evaluations"]),
        "density_solves": int(result["density_solves"]),
        "fused_state_evaluations": int(result["fused_state_evaluations"]),
        "callback_wall_time_s": float(result["callback_wall_time_s"]),
    }


def _native_args_sequence(payloads):
    return [create_struct(check_association(dict(payload))) for payload in payloads]


def _fit_generic_native_ceres(
    fixed_payloads,
    records,
    target_kinds,
    target_indices,
    target_indices_2,
    x0,
    lower,
    upper,
    max_nfev=200,
):
    max_nfev = int(max_nfev)
    if max_nfev < 1:
        raise InputError("Native Ceres generic regression requires max_nfev >= 1.")
    result = _core._fit_generic_native_ceres(
        _native_args_sequence(fixed_payloads),
        list(records),
        np.asarray(target_kinds, dtype=int),
        np.asarray(target_indices, dtype=int),
        np.asarray(target_indices_2, dtype=int),
        np.asarray(x0, dtype=float),
        np.asarray(lower, dtype=float),
        np.asarray(upper, dtype=float),
        max_nfev,
    )
    return {
        "x": vector_to_array(result["x"]),
        "cost": float(result["cost"]),
        "residual_norm": float(result["residual_norm"]),
        "initial_cost": float(result["initial_cost"]),
        "initial_residual_norm": float(result["initial_residual_norm"]),
        "metrics_by_term": {str(k): float(v) for k, v in dict(result["metrics_by_term"]).items()},
        "success": bool(result["success"]),
        "status": int(result["status"]),
        "nfev": int(result["nfev"]),
        "iterations": int(result["iterations"]),
        "starts_tried": int(result["starts_tried"]),
        "message": str(result["message"]),
        "backend": str(result["backend"]),
        "optimizer_backend": str(result["optimizer_backend"]),
        "derivative_backend": str(result["derivative_backend"]),
        "jacobian_available": bool(result.get("jacobian_available", True)),
        "jacobian_backend": str(result.get("jacobian_backend", "unspecified")),
    }


def _evaluate_generic_native_debug(
    fixed_payloads,
    records,
    target_kinds,
    target_indices,
    target_indices_2,
    x,
):
    result = _core._evaluate_generic_native_debug(
        _native_args_sequence(fixed_payloads),
        list(records),
        np.asarray(target_kinds, dtype=int),
        np.asarray(target_indices, dtype=int),
        np.asarray(target_indices_2, dtype=int),
        np.asarray(x, dtype=float),
    )
    return {
        "cost": float(result["cost"]),
        "residual_norm": float(result["residual_norm"]),
        "residuals": vector_to_array(result["residuals"]),
        "metrics_by_term": {str(k): float(v) for k, v in dict(result["metrics_by_term"]).items()},
        "jacobian_row_major": vector_to_array(result.get("jacobian_row_major", [])),
        "jacobian_shape": tuple(result.get("jacobian_shape", (len(result["residuals"]), 0))),
        "jacobian_available": bool(result.get("jacobian_available", True)),
        "jacobian_backend": str(result.get("jacobian_backend", "unspecified")),
    }
