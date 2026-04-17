# -*- coding: utf-8 -*-
# setuptools: language=c++

import numpy as np
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr
from copy import deepcopy
from . cimport epcsaft
from ._types import ActivityCoefficientResult
from ._types import InputError
from ._types import SolutionError
from ._types import phase_to_int
from ._types import vector_to_array

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


def _sum_vector_terms(terms):
    total = None
    for name in _CONTRIBUTION_NAMES:
        value = np.asarray(terms[name], dtype=float)
        total = value.copy() if total is None else total + value
    return total


cdef dict _scalar_terms_dict(ScalarContributionTerms terms):
    return {
        "hc": float(terms.hc),
        "disp": float(terms.disp),
        "assoc": float(terms.assoc),
        "ion": float(terms.ion),
        "born": float(terms.born),
        "total": float(terms.total),
    }


cdef dict _vector_terms_dict(VectorContributionTerms terms, Py_ssize_t expected_size, str label):
    cdef dict out = {}
    cdef dict blocks = {
        "hc": vector_to_array(terms.hc),
        "disp": vector_to_array(terms.disp),
        "assoc": vector_to_array(terms.assoc),
        "ion": vector_to_array(terms.ion),
        "born": vector_to_array(terms.born),
        "total": vector_to_array(terms.total),
    }
    cdef object name
    cdef object arr
    for name, arr in blocks.items():
        arr = np.asarray(arr, dtype=float)
        if arr.size != expected_size:
            raise SolutionError("Unexpected {} payload size for {}: expected {}, got {}.".format(label, name, int(expected_size), int(arr.size)))
        out[name] = arr
    return out


cdef class ePCSAFTMixture:
    """Native-backed ePC-SAFT parameter model and state factory."""
    cdef shared_ptr[ePCSAFTMixtureNative] _native
    cdef object _params
    cdef object _species

    def __init__(self, params=None, species=None):
        """Create a mixture from a resolved parameter payload."""
        self._native = shared_ptr[ePCSAFTMixtureNative]()
        self._params = None
        self._species = None
        if params is not None:
            self._init_from_params(params, species)

    @classmethod
    def from_params(cls, params, species=None):
        """Construct a mixture from an already-resolved parameter dict."""
        return cls(params=params, species=species)

    @classmethod
    def from_dataset(cls, dataset_name, species, x, T, user_options=None):
        """Construct a mixture by resolving packaged dataset parameters."""
        from .parameters import get_prop_dict

        params = get_prop_dict(dataset_name, species, x, T, user_options=user_options)
        return cls(params=params, species=species)

    cdef void _init_from_params(self, params, species=None):
        """Initialize the native mixture from a normalized parameter payload."""
        params = check_association(params)
        cppargs = create_struct(params)
        self._native = shared_ptr[ePCSAFTMixtureNative](new ePCSAFTMixtureNative(cppargs))
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
        return int(self._native.get().ncomp())

    def clear_runtime_caches(self):
        """Clear internal runtime caches used for repeated state/reference evaluations."""
        self._native.get().clear_runtime_caches()

    def reset_runtime_cache_stats(self):
        """Reset runtime cache hit/fallback counters without clearing cached values."""
        self._native.get().reset_runtime_cache_stats()

    def runtime_cache_stats(self):
        """Return native runtime cache counters for profiling and validation."""
        return {
            "reference_state_cache_hits": int(self._native.get().reference_state_cache_hits()),
            "reference_state_cache_misses": int(self._native.get().reference_state_cache_misses()),
            "density_warm_start_hits": int(self._native.get().density_warm_start_hits()),
            "density_warm_start_fallbacks": int(self._native.get().density_warm_start_fallbacks()),
        }

    def state(self, T, x, P=None, rho=None, phase="liq"):
        """Create an immutable thermodynamic state for the mixture.

        States built from pressure resolve and cache density during construction.
        """
        if (P is None) == (rho is None):
            raise InputError("Provide exactly one of P or rho when constructing a state.")
        return ePCSAFTState(self, T, x, P=P, rho=rho, phase=phase)

    def __repr__(self):
        """Return a short debugging representation of the mixture."""
        return f"ePCSAFTMixture(ncomp={self.ncomp}, species={self._species})"


cdef class ePCSAFTState:
    """Immutable thermodynamic state bound to one mixture."""
    cdef shared_ptr[ePCSAFTStateNative] _native
    cdef object _mixture
    cdef object _x
    cdef double _T
    cdef object _P
    cdef object _rho
    cdef int _phase

    def __init__(self, mixture, T, x, P=None, rho=None, phase="liq"):
        """Create a state with exactly one intensive variable fixed.

        Pressure-based states solve the internal T, P, x -> rho closure eagerly.
        """
        if not isinstance(mixture, ePCSAFTMixture):
            raise InputError("mixture must be a ePCSAFTMixture instance.")
        cdef ePCSAFTMixture mix = mixture
        x, params = ensure_numpy_input(x, mix._params)
        # ensure_numpy_input may normalize a scalar mixture parameter path, but the
        # state should retain the original mixture data unchanged.
        phase_num = phase_to_int(phase)
        has_p = P is not None
        has_rho = rho is not None
        if has_p == has_rho:
            raise InputError("Provide exactly one of P or rho when constructing a state.")
        if has_p:
            check_input(x, {"temperature": T, "pressure": P})
        else:
            check_input(x, {"temperature": T, "density": rho})
        cpp_x = np_to_vector_double(x)
        try:
            self._native = shared_ptr[ePCSAFTStateNative](new ePCSAFTStateNative(
                mix._native,
                float(T),
                cpp_x,
                phase_num,
                has_p,
                float(P) if P is not None else 0.0,
                has_rho,
                float(rho) if rho is not None else 0.0,
            ))
        except Exception as exc:
            if has_p:
                raise SolutionError(str(exc))
            raise
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
        return float(self._native.get().pressure())

    def density(self, units="molar"):
        """Return the state density in the requested units."""
        token = str(units or "molar").strip().lower()
        if token in {"molar", "mol", "molar_density", "mol/m^3", "mol/m3"}:
            return float(self._native.get().density())
        if token in {"mass", "mass_density", "kg/m^3", "kg/m3"}:
            return self.mass_density()
        raise InputError("density units must be 'molar'/'mol/m^3' or 'mass'/'kg/m^3'.")

    def molar_density(self):
        """Return the molar density of the bound state in mol/m^3."""
        return float(self._native.get().density())

    def mass_density(self):
        """Return the mass density of the bound state in kg/m^3."""
        cdef ePCSAFTMixture mix = self._mixture
        mw = np.asarray(mix._params.get("MW", []), dtype=float).flatten()
        if mw.size == 0:
            raise InputError("Mass density requires component molecular weights in kg/mol.")
        if mw.size != self._x.size:
            raise InputError("Mass density requires one molecular-weight value per component.")
        return float(self.molar_density() * float(np.dot(np.asarray(self._x, dtype=float), mw)))

    def _temperature_derivative_residual_helmholtz_term_result(self):
        cdef ScalarContributionTerms result = self._native.get().temperature_derivative_residual_helmholtz_result()
        return _scalar_terms_dict(result)

    def _composition_derivative_residual_helmholtz_result(self):
        ncomp = int(self._x.size)
        cdef CompositionContributionResult result = self._native.get().composition_derivative_residual_helmholtz_result()
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
        }

    def _fugacity_coefficient_term_result(self):
        ncomp = int(self._x.size)
        cdef FugacityContributionResult result = self._native.get().fugacity_coefficient_result()
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
            return float(self._native.get().compressibility_factor())
        cdef CompressibilityFactorResult result = self._native.get().compressibility_factor_result()
        payload_dict = _scalar_terms_dict(result.terms)
        terms = {name: float(payload_dict[name]) for name in _CONTRIBUTION_NAMES}
        terms["ideal"] = 1.0
        return {
            "total": float(payload_dict["total"]),
            "terms": terms,
        }

    def residual_helmholtz(self, return_contribution_terms=False):
        """Return the residual Helmholtz energy."""
        if not return_contribution_terms:
            return float(self._native.get().residual_helmholtz())
        cdef ScalarContributionTerms result = self._native.get().residual_helmholtz_result()
        payload_dict = _scalar_terms_dict(result)
        terms = {name: float(payload_dict[name]) for name in _CONTRIBUTION_NAMES}
        return {
            "total": float(payload_dict["total"]),
            "terms": terms,
        }

    def temperature_derivative_residual_helmholtz(self, return_contribution_terms=False):
        """Return the temperature derivative of the residual Helmholtz energy."""
        if not return_contribution_terms:
            return float(self._native.get().temperature_derivative_residual_helmholtz())
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
        return float(self._native.get().residual_enthalpy())

    def residual_entropy(self):
        """Return the residual entropy."""
        return float(self._native.get().residual_entropy())

    def residual_gibbs(self):
        """Return the residual Gibbs energy."""
        return float(self._native.get().residual_gibbs())

    def residual_chemical_potential(self, return_contribution_terms=False):
        """Return the residual chemical potentials."""
        if not return_contribution_terms:
            return vector_to_array(self._native.get().residual_chemical_potential())
        cdef ResidualChemicalPotentialResult result = self._native.get().residual_chemical_potential_result()
        payload_dict = _vector_terms_dict(result.mu, int(self._x.size), "mu")
        terms = {name: np.asarray(payload_dict[name], dtype=float) for name in _CONTRIBUTION_NAMES}
        return {
            "total": np.asarray(payload_dict["total"], dtype=float),
            "terms": terms,
        }

    def _activity_coefficient_bundle(self, species=None, solvent=None, include_aux=False):
        """Build the full activity-coefficient payload for internal reuse."""
        species = self._mixture.species if species is None else [str(s) for s in species]
        if len(species) != self._x.size:
            raise InputError("species length ({}) must match composition length ({}).".format(len(species), self._x.size))
        has_solvent_override, solvent_index = _resolve_solvent_override(self._mixture, species, solvent)
        cdef bint include_aux_c = bool(include_aux)
        cdef bint has_solvent_override_c = bool(has_solvent_override)
        cdef int solvent_index_c = int(solvent_index)
        cdef ActivityCoefficientNative out = self._native.get().activity_coefficient_native(include_aux_c, has_solvent_override_c, solvent_index_c)
        pair_cat = np.asarray(out.pair_cation_indices, dtype=int)
        pair_an = np.asarray(out.pair_anion_indices, dtype=int)
        pair_labels = tuple(species[int(ic)] + species[int(ia)] for ic, ia in zip(pair_cat.tolist(), pair_an.tolist()))
        ion_idx = np.sort(np.unique(np.concatenate([np.asarray(out.cation_indices, dtype=int), np.asarray(out.anion_indices, dtype=int)])))
        ion_labels = tuple(species[int(i)] for i in ion_idx.tolist())
        return ActivityCoefficientResult(
            species=tuple(species),
            component_activity_coefficients=np.asarray(out.component_activity_coefficients, dtype=float),
            solvation_free_energy_values=np.asarray(out.solvation_free_energy, dtype=float),
            mean_ionic_activity_coefficients_mole_fraction_values=np.asarray(out.mean_ionic_activity_coefficients_mole_fraction, dtype=float),
            mean_ionic_activity_coefficients_molality_values=np.asarray(out.mean_ionic_activity_coefficients_molality, dtype=float),
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
            raise InputError("species length ({}) must match composition length ({}).".format(len(species), self._x.size))
        has_solvent_override, solvent_index = _resolve_solvent_override(self._mixture, species, solvent)
        cdef bint include_aux_c = False
        cdef bint has_solvent_override_c = bool(has_solvent_override)
        cdef int solvent_index_c = int(solvent_index)
        cdef ActivityCoefficientNative out = self._native.get().activity_coefficient_native(include_aux_c, has_solvent_override_c, solvent_index_c)
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

    def fugacity_coefficient(self, natural_log=True, return_contribution_terms=False):
        """Return fugacity coefficients, defaulting to natural-log form."""
        ln_total = vector_to_array(self._native.get().ln_fugacity_coefficient())
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

    def relative_permittivity(self):
        """Return the dielectric model evaluation for the current state."""
        flat = vector_to_array(self._native.get().relative_permittivity())
        return flat[0], flat[1:]

    def osmotic_coefficient(self):
        """Return the osmotic coefficient."""
        return np.asarray([self._native.get().osmotic_coefficient()], dtype=float)

    def state_diagnostics(self, species=None):
        """Return a diagnostic dictionary of the main state properties."""
        cdef ePCSAFTMixture mix = self._mixture
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
        raise AttributeError("{} has no attribute '{}'".format(type(self).__name__, name))

    def __repr__(self):
        """Return a short debugging representation of the state."""
        return f"ePCSAFTState(T={self._T}, phase={self._phase}, x={self._x})"


cdef double _epcsaft_density_checked(double t, double p, vector[double] x, int phase, add_args &cppargs) except *:
    try:
        return den_cpp(t, p, x, phase, cppargs)
    except Exception as exc:
        raise SolutionError(str(exc))

def check_input(x, vars):
    if abs(np.sum(x) - 1) > 1e-7:
        raise InputError('The mole fractions do not sum to 1. x = {}'.format(x))
    if 'temperature' in vars:
        if vars['temperature'] <= 0:
            raise InputError('The {} must be a positive number. {} = {}'.format('temperature', 'temperature', vars['temperature']))
    if 'density' in vars:
        if vars['density'] <= 0:
            raise InputError('The {} must be a positive number. {} = {}'.format('density', 'density', vars['density']))
    if 'pressure' in vars:
        if vars['pressure'] <= 0:
            raise InputError('The {} must be a positive number. {} = {}'.format('pressure', 'pressure', vars['pressure']))
    if 'Q' in vars:
        if (vars['Q'] < 0) or (vars['Q'] > 1):
            raise InputError('{} must be <= 1 and >= 0. {} = {}'.format('Q', 'Q', vars['Q']))

def check_association(params):
    params = deepcopy(params)
    if ('e_assoc' in params) and ('vol_a' not in params):
        raise InputError('e_assoc was given, but not vol_a.')
    elif ('vol_a' in params) and ('e_assoc' not in params):
        raise InputError('vol_a was given, but not e_assoc.')

    if ('e_assoc' in params) and ('assoc_scheme' not in params):
        params['assoc_scheme'] = []
        for a in params['vol_a']:
            if a != 0:
                params['assoc_scheme'].append('2b')
            else:
                params['assoc_scheme'].append(None)

    if ('e_assoc' in params):
        params = create_assoc_matrix(params)

    return params

def create_assoc_matrix(params):
    charge = [] # whether the association site has a partial positive charge (i.e. hydrogen), negative charge, or elements of both (e.g. for acids modelled as type 1)

    scheme_charges = {
        '1': [0],
        '2a': [0, 0],
        '2b': [-1, 1],
        '3a': [0, 0, 0],
        '3b': [-1, -1, 1],
        '4a': [0, 0, 0, 0],
        '4b': [1, 1, 1, -1],
        '4c': [-1, -1, 1, 1]
    }

    assoc_num = []
    for comp in params['assoc_scheme']:
        if comp is None:
            assoc_num.append(0)
            pass
        elif type(comp) is list:
            num = 0
            for site in comp:
                if site.lower() not in scheme_charges:
                    raise InputError('{} is not a valid association type.'.format(site))
                charge.extend(scheme_charges[site.lower()])
                num += len(scheme_charges[site.lower()])
            assoc_num.append(num)
        else:
            if comp.lower() not in scheme_charges:
                raise InputError('{} is not a valid association type.'.format(comp))
            charge.extend(scheme_charges[comp.lower()])
            assoc_num.append(len(scheme_charges[comp.lower()]))
    params['assoc_num'] = np.asarray(assoc_num)

    params['assoc_matrix'] = np.zeros((len(charge)*len(charge)))
    ctr = 0
    for c1 in charge:
        for c2 in charge:
            if (c1 == 0 or c2 == 0):
                params['assoc_matrix'][ctr] = 1;
            elif (c1 == 1 and c2 == -1):
                params['assoc_matrix'][ctr] = 1;
            elif (c1 == -1 and c2 == 1):
                params['assoc_matrix'][ctr] = 1;
            else:
                params['assoc_matrix'][ctr] = 0;
            ctr += 1

    return params

def ensure_numpy_input(x, params):
    if np.isscalar(x):
        x = np.asarray([x], dtype=float)
    if np.isscalar(params['m']):
        params['m'] = np.asarray([params['m']], dtype=float)
    if np.isscalar(params['s']):
        params['s'] = np.asarray([params['s']], dtype=float)
    if np.isscalar(params['e']):
        params['e'] = np.asarray([params['e']], dtype=float)
    return x, params

def _safe_unit_log(values, floor=1e-300):
    vals = np.asarray(values, dtype=float)
    return np.log(np.maximum(vals, floor))


def _resolve_solvent_override(mixture, species, solvent):
    cdef ePCSAFTMixture mix = mixture
    z = np.asarray(mix._params.get("z", []), dtype=float).flatten()
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
            raise InputError("Unknown solvent label '{}'. Available species={}".format(token, list(species)))
        idx = species.index(token)
    if idx < 0 or idx >= z.size:
        raise InputError("solvent index out of bounds: {} (ncomp={})".format(idx, int(z.size)))
    if abs(float(z[idx])) > 1e-12:
        raise InputError("solvent override must reference a neutral species (z=0).")
    return True, idx

def np_to_vector_double(np_array):
    """Take a numpy array and return a C++ vector."""
    cdef vector[double] cpp_vector

    try:
        np_array = np_array.flatten()
        N = np_array.shape[0]
        for i in range(N):
            cpp_vector.push_back(np_array[i])
    except TypeError:
        cpp_vector.push_back(np_array)

    return cpp_vector

def np_to_vector_int(np_array):
    """Take a numpy array and return a C++ vector."""
    cdef vector[int] cpp_vector

    try:
        np_array = np_array.flatten()
        N = np_array.shape[0]
        for i in range(N):
            cpp_vector.push_back(np_array[i])
    except TypeError:
        cpp_vector.push_back(np_array)

    return cpp_vector

def create_struct(params):
    """Convert ePC-SAFT parameters to a C++ struct."""
    cdef add_args cppargs
    cdef int ncomp

    for removed_key in ('dipm', 'dip_num'):
        if removed_key in params:
            raise ValueError(
                'Removed polar parameter "{}" is not supported by the active ePC-SAFT package.'.format(removed_key)
            )

    cppargs.mixed_rel_perm_water_index = -1
    cppargs.m = np_to_vector_double(params['m'])
    ncomp = len(np.asarray(params['m']).flatten())
    cppargs.s = np_to_vector_double(params['s'])
    cppargs.e = np_to_vector_double(params['e'])
    if 'k_ij' in params:
        cppargs.k_ij = np_to_vector_double(params['k_ij'])
    if ('e_assoc' in params) and np.any(params['e_assoc']):
        cppargs.e_assoc = np_to_vector_double(params['e_assoc'])
    if ('vol_a' in params) and np.any(params['vol_a']):
        cppargs.vol_a = np_to_vector_double(params['vol_a'])
    z_arr = None
    if 'z' in params:
        z_arr = np.asarray(params['z'], dtype=float).flatten()
        if z_arr.size not in (0, ncomp):
            raise ValueError('params["z"] must have length {} (or be empty), got {}.'.format(ncomp, z_arr.size))
        if z_arr.size == ncomp:
            cppargs.z = np_to_vector_double(z_arr)
    if 'dielc' in params:
        dielc_arr = np.asarray(params['dielc'], dtype=float).flatten()
        if dielc_arr.size != ncomp:
            raise ValueError('params["dielc"] must have length {}, got {}.'.format(ncomp, dielc_arr.size))
        cppargs.dielc = np_to_vector_double(dielc_arr)
    if 'MW' in params:
        mw_arr = np.asarray(params['MW'], dtype=float).flatten()
        if mw_arr.size != ncomp:
            raise ValueError('params["MW"] must have length {}, got {}.'.format(ncomp, mw_arr.size))
        cppargs.mw = np_to_vector_double(mw_arr)
    if 'mixed_rel_perm_a' in params:
        mixed_a_arr = np.asarray(params['mixed_rel_perm_a'], dtype=float).flatten()
        if mixed_a_arr.size != ncomp:
            raise ValueError('params["mixed_rel_perm_a"] must have length {}, got {}.'.format(ncomp, mixed_a_arr.size))
        cppargs.mixed_rel_perm_a = np_to_vector_double(mixed_a_arr)
    if 'mixed_rel_perm_b' in params:
        mixed_b_arr = np.asarray(params['mixed_rel_perm_b'], dtype=float).flatten()
        if mixed_b_arr.size != ncomp:
            raise ValueError('params["mixed_rel_perm_b"] must have length {}, got {}.'.format(ncomp, mixed_b_arr.size))
        cppargs.mixed_rel_perm_b = np_to_vector_double(mixed_b_arr)
    if 'mixed_rel_perm_c' in params:
        mixed_c_arr = np.asarray(params['mixed_rel_perm_c'], dtype=float).flatten()
        if mixed_c_arr.size != ncomp:
            raise ValueError('params["mixed_rel_perm_c"] must have length {}, got {}.'.format(ncomp, mixed_c_arr.size))
        cppargs.mixed_rel_perm_c = np_to_vector_double(mixed_c_arr)
    if 'mixed_rel_perm_mask' in params:
        mixed_mask_arr = np.asarray(params['mixed_rel_perm_mask'], dtype=int).flatten()
        if mixed_mask_arr.size != ncomp:
            raise ValueError('params["mixed_rel_perm_mask"] must have length {}, got {}.'.format(ncomp, mixed_mask_arr.size))
        cppargs.mixed_rel_perm_mask = np_to_vector_int(mixed_mask_arr)
    if 'mixed_rel_perm_water_index' in params:
        cppargs.mixed_rel_perm_water_index = int(params['mixed_rel_perm_water_index'])
    if cppargs.z.size() > 0 and cppargs.dielc.size() == 0:
        raise ValueError('Electrolyte parameters require params["dielc"] as a per-species array.')
    d_born_arr = None
    if 'd_born' in params:
        d_born_arr = np.asarray(params['d_born'], dtype=float).flatten()
        if d_born_arr.size != ncomp:
            raise ValueError('params["d_born"] must have length {}, got {}.'.format(ncomp, d_born_arr.size))
        cppargs.d_born = np_to_vector_double(d_born_arr)
    if 'f_solv' in params:
        cppargs.f_solv = np_to_vector_double(np.asarray(params['f_solv'], dtype=float))

    unsupported_flat_elec_keys = {
        'dielc_rule', 'dielc_diff_mode', 'born_model', 'born_radius_model',
        'born_diff_mode', 'born_eps_mode', 'DH_model', 'bjeruum_treatment',
        'd_ion_mode', 'include_born_model', 'd_Born_mode',
        'born_solvation_shell_model', 'born_dielectric_saturation', 'born_bulk_mode',
        'mu_DH_diff_mode', 'mu_DH_comp_dep_rel_perm', 'mu_DH_include_sum_term',
        'mu_born_diff_mode', 'mu_born_comp_dep_rel_perm', 'mu_born_include_sum_term', 'mu_born_comp_dep_delta_d'
    }
    if any((k in params) for k in unsupported_flat_elec_keys):
        raise ValueError(
            'Flat electrostatic params are no longer supported; provide nested params["elec_model"] schema.'
        )

    elec_model = params.get('elec_model', None)
    if elec_model is None:
        if cppargs.z.size() > 0:
            # Apply the canonical electrolyte defaults when ionic parameters are present.
            elec_model = {
                'rel_perm': {'rule': 1, 'differential_mode': 'analytical'},
                'hc_model': {'dadx_differential_mode': 'analytical'},
                'disp_model': {'dadx_differential_mode': 'analytical'},
                'assoc_model': {'dadx_differential_mode': 'analytical'},
                'DH_model': {
                    'd_ion_mode': 1,
                    'bjeruum_treatment': False,
                    'mu_DH_model': {
                        'differential_mode': 'analytical',
                        'comp_dep_rel_perm': True,
                        'include_sum_term': True,
                    },
                },
                'include_born_model': True,
                'born_model': {
                    'd_Born_mode': 0,
                    'solvation_shell_model': False,
                    'dielectric_saturation': False,
                    'bulk_mode': 'mix',
                    'mu_born_model': {
                        'differential_mode': 'analytical',
                        'comp_dep_rel_perm': True,
                        'include_sum_term': True,
                        'comp_dep_delta_d': False,
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
            raise ValueError('{} contains unsupported key(s): {}.'.format(label, unknown))

    _reject_unknown_keys(
        elec_model,
        {'rel_perm', 'hc_model', 'disp_model', 'assoc_model', 'DH_model', 'include_born_model', 'born_model'},
        'params["elec_model"]'
    )

    def _as_bool(v):
        if isinstance(v, bool):
            return v
        if isinstance(v, (int, np.integer)):
            return bool(v)
        if isinstance(v, str):
            s = v.strip().lower()
            if s in {'1', 'true', 'yes', 'y', 'on'}:
                return True
            if s in {'0', 'false', 'no', 'n', 'off'}:
                return False
        raise ValueError('Could not coerce value to bool: {}'.format(v))

    def _as_int_alias(v, aliases):
        if isinstance(v, (int, np.integer)):
            return int(v)
        if isinstance(v, str):
            s = v.strip().lower()
            if s in aliases:
                return int(aliases[s])
            if s.isdigit() or (s.startswith('-') and s[1:].isdigit()):
                return int(s)
        raise ValueError('Unknown option value: {}'.format(v))

    rule_alias = {
        'constant': 0,
        'rule0': 0,
        'linear': 1,
        'linear-molefraction': 1,
        'linear-mixing-mole': 1,
        'rule1': 1,
        'linear-massfraction': 2,
        'linear-mixing-weight': 2,
        'rule2': 2,
        'combined': 3,
        'rule3': 3,
        'empirical': 4,
        'rule4': 4,
        'rule5': 5,
        'rule6': 6,
        'aqueous-organic': 8,
        'aqueous_organic': 8,
        'mixed-aqueous-organic': 8,
        'mixed_aqueous_organic': 8,
        'rule8': 8,
    }
    diff_alias = {
        'analytic': 0,
        'analytical': 0,
        'numeric': 1,
        'numerical': 1,
        'finite_difference': 1,
        'finite-difference': 1,
        'finite difference': 1,
        'fd': 1,
        'autodiff': 2,
        'automatic_differentiation': 2,
        'automatic-differentiation': 2,
        'automatic differentiation': 2,
    }
    d_ion_alias = {'t_indep': 0, 't_dep_1': 1, 't_dep_2': 2}
    d_born_alias = {'t_indep': 0, 't_dep_1': 1, 't_dep_2': 2, 'fitted_param': 3}
    bulk_alias = {'mix': 0, 'bulk': 0, 'solvent': 1}

    rel_perm = elec_model.get('rel_perm', {})
    if not isinstance(rel_perm, dict):
        raise ValueError('params["elec_model"]["rel_perm"] must be a dict.')
    _reject_unknown_keys(rel_perm, {'rule', 'differential_mode'}, 'params["elec_model"]["rel_perm"]')
    hc_model_dict = elec_model.get('hc_model', {})
    if not isinstance(hc_model_dict, dict):
        raise ValueError('params["elec_model"]["hc_model"] must be a dict.')
    _reject_unknown_keys(hc_model_dict, {'dadx_differential_mode'}, 'params["elec_model"]["hc_model"]')
    disp_model_dict = elec_model.get('disp_model', {})
    if not isinstance(disp_model_dict, dict):
        raise ValueError('params["elec_model"]["disp_model"] must be a dict.')
    _reject_unknown_keys(disp_model_dict, {'dadx_differential_mode'}, 'params["elec_model"]["disp_model"]')
    assoc_model_dict = elec_model.get('assoc_model', {})
    if not isinstance(assoc_model_dict, dict):
        raise ValueError('params["elec_model"]["assoc_model"] must be a dict.')
    _reject_unknown_keys(assoc_model_dict, {'dadx_differential_mode'}, 'params["elec_model"]["assoc_model"]')
    dh_model_dict = elec_model.get('DH_model', {})
    if not isinstance(dh_model_dict, dict):
        raise ValueError('params["elec_model"]["DH_model"] must be a dict.')
    _reject_unknown_keys(dh_model_dict, {'d_ion_mode', 'bjeruum_treatment', 'mu_DH_model'}, 'params["elec_model"]["DH_model"]')
    born_model_dict = elec_model.get('born_model', {})
    if not isinstance(born_model_dict, dict):
        raise ValueError('params["elec_model"]["born_model"] must be a dict.')
    _reject_unknown_keys(
        born_model_dict,
        {'d_Born_mode', 'solvation_shell_model', 'dielectric_saturation', 'bulk_mode', 'mu_born_model'},
        'params["elec_model"]["born_model"]'
    )
    mu_dh = dh_model_dict.get('mu_DH_model', {})
    if not isinstance(mu_dh, dict):
        raise ValueError('params["elec_model"]["DH_model"]["mu_DH_model"] must be a dict.')
    _reject_unknown_keys(
        mu_dh,
        {'differential_mode', 'comp_dep_rel_perm', 'include_sum_term'},
        'params["elec_model"]["DH_model"]["mu_DH_model"]'
    )
    mu_born = born_model_dict.get('mu_born_model', {})
    if not isinstance(mu_born, dict):
        raise ValueError('params["elec_model"]["born_model"]["mu_born_model"] must be a dict.')
    _reject_unknown_keys(
        mu_born,
        {'differential_mode', 'comp_dep_rel_perm', 'include_sum_term', 'comp_dep_delta_d'},
        'params["elec_model"]["born_model"]["mu_born_model"]'
    )

    cppargs.dielc_rule = _as_int_alias(rel_perm.get('rule', 1), rule_alias)
    cppargs.dielc_diff_mode = _as_int_alias(rel_perm.get('differential_mode', 'analytical'), diff_alias)
    if cppargs.dielc_diff_mode not in (0, 1, 2):
        raise ValueError('Unknown rel_perm differential_mode. Supported values are analytical/numerical/autodiff (0/1/2).')
    cppargs.hc_dadx_diff_mode = _as_int_alias(hc_model_dict.get('dadx_differential_mode', 'analytical'), diff_alias)
    if cppargs.hc_dadx_diff_mode not in (0, 1, 2):
        raise ValueError('Unknown hc_model dadx_differential_mode. Supported values are analytical/numerical/autodiff (0/1/2).')
    cppargs.disp_dadx_diff_mode = _as_int_alias(disp_model_dict.get('dadx_differential_mode', 'analytical'), diff_alias)
    if cppargs.disp_dadx_diff_mode not in (0, 1, 2):
        raise ValueError('Unknown disp_model dadx_differential_mode. Supported values are analytical/numerical/autodiff (0/1/2).')
    cppargs.assoc_dadx_diff_mode = _as_int_alias(assoc_model_dict.get('dadx_differential_mode', 'analytical'), diff_alias)
    if cppargs.assoc_dadx_diff_mode not in (0, 1, 2):
        raise ValueError('Unknown assoc_model dadx_differential_mode. Supported values are analytical/numerical/autodiff (0/1/2).')
    if cppargs.dielc_rule < 0 or cppargs.dielc_rule > 8:
        raise ValueError('Unknown rel_perm rule. Supported values are 0..8.')

    cppargs.d_ion_mode = _as_int_alias(dh_model_dict.get('d_ion_mode', 1), d_ion_alias)
    if cppargs.d_ion_mode not in (0, 1, 2):
        raise ValueError('Unknown d_ion_mode. Supported values are 0,1,2.')
    bjeruum = _as_bool(dh_model_dict.get('bjeruum_treatment', False))
    cppargs.mu_DH_diff_mode = _as_int_alias(mu_dh.get('differential_mode', 'analytical'), diff_alias)
    if cppargs.mu_DH_diff_mode not in (0, 1, 2):
        raise ValueError('Unknown mu_DH differential_mode. Supported values are analytical/numerical/autodiff (0/1/2).')
    cppargs.mu_DH_comp_dep_rel_perm = int(_as_bool(mu_dh.get('comp_dep_rel_perm', True)))
    cppargs.mu_DH_include_sum_term = int(_as_bool(mu_dh.get('include_sum_term', True)))

    cppargs.include_born_model = int(_as_bool(elec_model.get('include_born_model', True)))
    cppargs.d_born_mode = _as_int_alias(born_model_dict.get('d_Born_mode', 0), d_born_alias)
    if cppargs.d_born_mode not in (0, 1, 2, 3):
        raise ValueError('Unknown d_Born_mode. Supported values are 0,1,2,3.')
    cppargs.born_solvation_shell_model = int(_as_bool(born_model_dict.get('solvation_shell_model', False)))
    cppargs.born_dielectric_saturation = int(_as_bool(born_model_dict.get('dielectric_saturation', False)))
    cppargs.born_bulk_mode = _as_int_alias(born_model_dict.get('bulk_mode', 'mix'), bulk_alias)
    cppargs.mu_born_diff_mode = _as_int_alias(mu_born.get('differential_mode', 'analytical'), diff_alias)
    if cppargs.mu_born_diff_mode not in (0, 1, 2):
        raise ValueError('Unknown mu_born differential_mode. Supported values are analytical/numerical/autodiff (0/1/2).')
    cppargs.mu_born_comp_dep_rel_perm = int(_as_bool(mu_born.get('comp_dep_rel_perm', True)))
    cppargs.mu_born_include_sum_term = int(_as_bool(mu_born.get('include_sum_term', True)))
    cppargs.mu_born_comp_dep_delta_d = int(_as_bool(mu_born.get('comp_dep_delta_d', False)))

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
        elif cppargs.mu_born_comp_dep_rel_perm == 0:
            cppargs.born_diff_mode = 3
        elif cppargs.mu_born_include_sum_term == 0:
            cppargs.born_diff_mode = 2
        else:
            cppargs.born_diff_mode = 0

    if cppargs.born_model == 1 and cppargs.born_radius_model == 5:
        raise ValueError('d_Born_mode="fitted_param" requires SSM/DS Born path (include_born_model=true and SSM or DS true).')

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
        raise ValueError("Bjerrum treatment is reserved and not implemented (DH_model=2).")
    if cppargs.DH_model < 0 or cppargs.DH_model > 2:
        raise ValueError("Unknown DH_model. Supported values are 0, 1, and reserved 2.")
    if 'assoc_num' in params:
        cppargs.assoc_num = np_to_vector_int(params['assoc_num'])
    if 'assoc_matrix' in params:
        cppargs.assoc_matrix = np_to_vector_int(params['assoc_matrix'])
    if 'k_hb' in params:
        cppargs.k_hb = np_to_vector_double(params['k_hb'])
    if 'l_ij' in params:
        cppargs.l_ij = np_to_vector_double(params['l_ij'])

    return cppargs






