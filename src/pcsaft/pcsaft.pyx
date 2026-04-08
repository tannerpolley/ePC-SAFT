# -*- coding: utf-8 -*-
# setuptools: language=c++

import math
import numpy as np
from libc.math cimport cosh
from libc.math cimport sinh
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr
from copy import deepcopy
from . cimport pcsaft
from .parameters import get_prop_dict
from ._types import ActivityCoeffResult
from ._types import FlashResult
from ._types import InputError
from ._types import PhaseResult
from ._types import SolutionError
from ._types import VaporizationResult
from ._types import phase_to_int
from ._types import vector_to_array


cdef inline double _aly_lee(double t, double c0, double c1, double c2, double c3, double c4) noexcept:
    return (c0 + c1 * (c2 / t / sinh(c2 / t)) ** 2 + c3 * (c4 / t / cosh(c4 / t)) ** 2) / 1000.0


cdef class PCSAFTMixture:
    """Native-backed PC-SAFT parameter model and state factory."""
    cdef shared_ptr[PCSAFTMixtureNative] _native
    cdef object _params
    cdef object _species

    def __init__(self, params=None, species=None):
        """Create a mixture from a resolved parameter payload."""
        self._native = shared_ptr[PCSAFTMixtureNative]()
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
        params = get_prop_dict(dataset_name, species, x, T, user_options=user_options)
        return cls(params=params, species=species)

    cdef void _init_from_params(self, params, species=None):
        """Initialize the native mixture from a normalized parameter payload."""
        params = check_association(params)
        cppargs = create_struct(params)
        self._native = shared_ptr[PCSAFTMixtureNative](new PCSAFTMixtureNative(cppargs))
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

    def state(self, T, x, P=None, rho=None, phase="liq"):
        """Create an immutable thermodynamic state for the mixture."""
        if (P is None) == (rho is None):
            raise InputError("Provide exactly one of P or rho when constructing a state.")
        return PCSAFTState(self, T, x, P=P, rho=rho, phase=phase)

    def __repr__(self):
        """Return a short debugging representation of the mixture."""
        return f"PCSAFTMixture(ncomp={self.ncomp}, species={self._species})"


cdef class PCSAFTState:
    """Immutable thermodynamic state bound to one mixture."""
    cdef shared_ptr[PCSAFTStateNative] _native
    cdef object _mixture
    cdef object _x
    cdef double _T
    cdef object _P
    cdef object _rho
    cdef int _phase

    def __init__(self, mixture, T, x, P=None, rho=None, phase="liq"):
        """Create a state with exactly one intensive variable fixed."""
        if not isinstance(mixture, PCSAFTMixture):
            raise InputError("mixture must be a PCSAFTMixture instance.")
        cdef PCSAFTMixture mix = mixture
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
        self._native = shared_ptr[PCSAFTStateNative](new PCSAFTStateNative(
            mix._native,
            float(T),
            cpp_x,
            phase_num,
            has_p,
            float(P) if P is not None else 0.0,
            has_rho,
            float(rho) if rho is not None else 0.0,
        ))
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

    def pressure(self):
        """Return the pressure of the bound state."""
        return float(self._native.get().pressure())

    def density(self):
        """Return the density of the bound state."""
        return float(self._native.get().density())

    def Z(self):
        """Return the compressibility factor."""
        return float(self._native.get().Z())

    def a_res(self):
        """Return the residual Helmholtz energy."""
        return float(self._native.get().a_res())

    def dadt(self):
        """Return the temperature derivative of the residual Helmholtz energy."""
        return float(self._native.get().dadt())

    def h_res(self):
        """Return the residual enthalpy."""
        return float(self._native.get().h_res())

    def s_res(self):
        """Return the residual entropy."""
        return float(self._native.get().s_res())

    def g_res(self):
        """Return the residual Gibbs energy."""
        return float(self._native.get().g_res())

    def mu_res(self):
        """Return the residual chemical potentials."""
        return vector_to_array(self._native.get().mu_res())

    def gamma(self):
        """Return activity coefficients from the residual chemical potentials."""
        return vector_to_array(self._native.get().gamma())

    def lnfugcoef(self):
        """Return logarithmic fugacity coefficients."""
        return vector_to_array(self._native.get().lnfugcoef())

    def lnfugcoef_terms(self):
        """Return the term-by-term fugacity coefficient decomposition."""
        flat = vector_to_array(self._native.get().lnfugcoef_terms())
        ncomp = int(self._x.size)
        expected = 20 * ncomp + 25
        if flat.size != expected:
            raise SolutionError("Unexpected lnfug term payload size: expected {}, got {}.".format(expected, int(flat.size)))
        blocks = flat[:20 * ncomp].reshape((20, ncomp))
        scalars = flat[20 * ncomp:]
        return {
            "mu_hc": blocks[0],
            "mu_disp": blocks[1],
            "mu_polar": blocks[2],
            "mu_assoc": blocks[3],
            "mu_ion": blocks[4],
            "mu_born": blocks[5],
            "mu_total": blocks[6],
            "lnfugcoef_total": blocks[7],
            "lnfugcoef_hc": blocks[8],
            "lnfugcoef_disp": blocks[9],
            "lnfugcoef_polar": blocks[10],
            "lnfugcoef_assoc": blocks[11],
            "lnfugcoef_ion": blocks[12],
            "lnfugcoef_born": blocks[13],
            "dadx_hc": blocks[14],
            "dadx_disp": blocks[15],
            "dadx_polar": blocks[16],
            "dadx_assoc": blocks[17],
            "dadx_ion": blocks[18],
            "dadx_born": blocks[19],
            "a_hc": float(scalars[0]),
            "a_disp": float(scalars[1]),
            "a_polar": float(scalars[2]),
            "a_assoc": float(scalars[3]),
            "a_ion": float(scalars[4]),
            "a_born": float(scalars[5]),
            "sum_x_dadx_hc": float(scalars[6]),
            "sum_x_dadx_disp": float(scalars[7]),
            "sum_x_dadx_polar": float(scalars[8]),
            "sum_x_dadx_assoc": float(scalars[9]),
            "sum_x_dadx_ion": float(scalars[10]),
            "sum_x_dadx_born": float(scalars[11]),
            "z_raw_hc": float(scalars[12]),
            "z_raw_disp": float(scalars[13]),
            "z_raw_polar": float(scalars[14]),
            "z_raw_assoc": float(scalars[15]),
            "z_raw_ion": float(scalars[16]),
            "z_raw_born": float(scalars[17]),
            "z_hc": float(scalars[18]),
            "z_disp": float(scalars[19]),
            "z_polar": float(scalars[20]),
            "z_assoc": float(scalars[21]),
            "z_ion": float(scalars[22]),
            "z_born": float(scalars[23]),
            "z_total": float(scalars[24]),
        }

    def dielectric_eval(self):
        """Return the dielectric model evaluation for the current state."""
        flat = vector_to_array(self._native.get().dielectric_eval())
        return flat[0], flat[1:]

    def osmoticC(self):
        """Return the osmotic coefficient."""
        return np.asarray([self._native.get().osmoticC()], dtype=float)

    def actcoeff(self, species=None, solvent=None):
        """Return a bundled activity-coefficient result for the state."""
        species = self._mixture.species if species is None else [str(s) for s in species]
        if len(species) != self._x.size:
            raise InputError("species length ({}) must match composition length ({}).".format(len(species), self._x.size))
        has_solvent_override, solvent_index = _resolve_solvent_override(self._mixture, species, solvent)
        cdef bint has_solvent_override_c = bool(has_solvent_override)
        cdef int solvent_index_c = int(solvent_index)
        cdef ActivityCoeffNative out = self._native.get().actcoeff(has_solvent_override_c, solvent_index_c)
        pair_cat = np.asarray(out.pair_cation_indices, dtype=int)
        pair_an = np.asarray(out.pair_anion_indices, dtype=int)
        pair_labels = tuple(species[int(ic)] + species[int(ia)] for ic, ia in zip(pair_cat.tolist(), pair_an.tolist()))
        ion_idx = np.sort(np.unique(np.concatenate([np.asarray(out.cation_indices, dtype=int), np.asarray(out.anion_indices, dtype=int)])))
        ion_labels = tuple(species[int(i)] for i in ion_idx.tolist())
        return ActivityCoeffResult(
            species=tuple(species),
            component_gamma=np.asarray(out.gamma_components, dtype=float),
            gsolv_values=np.asarray(out.gsolv, dtype=float),
            gamma_mean_ionic_x_values=np.asarray(out.gamma_mean_ionic_x, dtype=float),
            gamma_mean_ionic_m_values=np.asarray(out.gamma_mean_ionic_m, dtype=float),
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
            osmotic_c=float(out.osmotic_c),
        )

    def breakdown(self, species=None):
        """Return a diagnostic dictionary of the main state properties."""
        cdef PCSAFTMixture mix = self._mixture
        species = self._mixture.species if species is None else species
        z = np.asarray(mix._params.get("z", []), dtype=float).flatten()
        has_ions = bool(np.any(np.abs(z) > 1e-12))
        terms = self.lnfugcoef_terms()
        lnfugcoef = self.lnfugcoef()
        gamma = np.exp(lnfugcoef)
        if has_ions:
            act = self.actcoeff(species=species)
            dielectric_eval = self.dielectric_eval()
            osmoticC = np.asarray([act.osmotic_c], dtype=float)
            miac_m = act.mean_ionic_m()
            miac = act.mean_ionic_x()
            gsolv = act.ion()
        else:
            dielectric_eval = None
            osmoticC = None
            miac_m = {}
            miac = {}
            gsolv = {}
        return {
            "T": self._T,
            "phase": self._phase,
            "x": np.asarray(self._x, dtype=float),
            "pressure": self.pressure(),
            "density": self.density(),
            "Z": self.Z(),
            "a_res": self.a_res(),
            "h_res": self.h_res(),
            "s_res": self.s_res(),
            "g_res": self.g_res(),
            "mu_res": self.mu_res(),
            "lnfugcoef": lnfugcoef,
            "gamma": gamma,
            "lnfugcoef_terms": terms,
            "dielectric_eval": dielectric_eval,
            "osmoticC": osmoticC,
            "miac_m": miac_m,
            "miac": miac,
            "gsolv": gsolv,
        }

    def miac_m(self, species=None):
        """Return mean-ionic activity coefficients on the molality basis."""
        return self.actcoeff(species=species).mean_ionic_m()

    def miac(self, species=None):
        """Return mean-ionic activity coefficients on the mole-fraction basis."""
        return self.actcoeff(species=species).mean_ionic_x()

    def gsolv(self, species=None):
        """Return ion solvation free-energy values keyed by species."""
        return self.actcoeff(species=species).ion()

    def flashTQ(self, q, p_guess=None):
        """Solve a temperature-quality flash at the state temperature."""
        if p_guess is None:
            out = self._native.get().flashTQ(float(q), False, 0.0)
        else:
            out = self._native.get().flashTQ(float(q), True, float(p_guess))
        xl = vector_to_array(out.xl)
        xv = vector_to_array(out.xv)
        pl = float(out.value)
        liq_state = PCSAFTState(self._mixture, self._T, xl, P=pl, phase="liq")
        vap_state = PCSAFTState(self._mixture, self._T, xv, P=pl, phase="vap")
        lnfug_l = liq_state.lnfugcoef() + np.log(np.maximum(liq_state.x, 1e-300)) + math.log(pl)
        lnfug_v = vap_state.lnfugcoef() + np.log(np.maximum(vap_state.x, 1e-300)) + math.log(pl)
        phases = [
            PhaseResult(1.0 - float(q), xl, liq_state.density(), liq_state.lnfugcoef(), lnfug_l),
            PhaseResult(float(q), xv, vap_state.density(), vap_state.lnfugcoef(), lnfug_v),
        ]
        return FlashResult(out.value, phases, "TQ")

    def flashPQ(self, p, q, t_guess=None):
        """Solve a pressure-quality flash with an optional temperature guess."""
        if t_guess is None:
            out = self._native.get().flashPQ(float(p), float(q), False, 0.0)
        else:
            out = self._native.get().flashPQ(float(p), float(q), True, float(t_guess))
        xl = vector_to_array(out.xl)
        xv = vector_to_array(out.xv)
        tt = float(out.value)
        liq_state = PCSAFTState(self._mixture, tt, xl, P=float(p), phase="liq")
        vap_state = PCSAFTState(self._mixture, tt, xv, P=float(p), phase="vap")
        lnfug_l = liq_state.lnfugcoef() + np.log(np.maximum(liq_state.x, 1e-300)) + math.log(float(p))
        lnfug_v = vap_state.lnfugcoef() + np.log(np.maximum(vap_state.x, 1e-300)) + math.log(float(p))
        phases = [
            PhaseResult(1.0 - float(q), xl, liq_state.density(), liq_state.lnfugcoef(), lnfug_l),
            PhaseResult(float(q), xv, vap_state.density(), vap_state.lnfugcoef(), lnfug_v),
        ]
        return FlashResult(out.value, phases, "PQ")

    def Hvap(self, p_guess=None):
        """Return the vaporization pressure for the current mixture state."""
        if p_guess is None:
            out = self._native.get().Hvap(False, 0.0)
        else:
            out = self._native.get().Hvap(True, float(p_guess))
        return VaporizationResult(out.value, out.pressure)

    def cp(self, aly_lee_params):
        """Estimate the isobaric heat capacity using an Aly-Lee fit."""
        cdef PCSAFTMixture mix = self._mixture
        aly_lee_params = np.asarray(aly_lee_params, dtype=float).flatten()
        if aly_lee_params.size != 5:
            raise InputError("aly_lee_params must have length 5.")
        cdef double c0 = float(aly_lee_params[0])
        cdef double c1 = float(aly_lee_params[1])
        cdef double c2 = float(aly_lee_params[2])
        cdef double c3 = float(aly_lee_params[3])
        cdef double c4 = float(aly_lee_params[4])
        cdef add_args cppargs = create_struct(check_association(mix._params))
        cdef object x = np.asarray(self._x, dtype=float)
        cdef double t = self._T
        cdef double rho = self.density()
        cdef int ph = 0 if rho > 900 else 1
        cdef double cp_ideal = _aly_lee(t, c0, c1, c2, c3, c4)
        cdef double p = p_cpp(t, rho, x, cppargs)
        cdef double rho0 = _pcsaft_den_checked(t - 0.001, p, x, ph, cppargs)
        cdef double hres0 = hres_cpp(t - 0.001, rho0, x, cppargs)
        cdef double rho1 = _pcsaft_den_checked(t + 0.001, p, x, ph, cppargs)
        cdef double hres1 = hres_cpp(t + 0.001, rho1, x, cppargs)
        cdef double dhdt = (hres1 - hres0) / 0.002
        return float(cp_ideal + dhdt)

    def __repr__(self):
        """Return a short debugging representation of the state."""
        return f"PCSAFTState(T={self._T}, phase={self._phase}, x={self._x})"


cdef double _pcsaft_den_checked(double t, double p, vector[double] x, int phase, add_args &cppargs) except *:
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
    cdef PCSAFTMixture mix = mixture
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
    """Convert PC-SAFT parameters to a C++ struct."""
    cdef add_args cppargs
    cdef int ncomp

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
    if ('dipm' in params) and np.any(params['dip_num']) and np.any(params['dipm']):
        cppargs.dipm = np_to_vector_double(params['dipm'])
    if ('dip_num' in params) and np.any(params['dip_num']):
        cppargs.dip_num = np_to_vector_double(params['dip_num'])
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

    legacy_elec_keys = {
        'dielc_rule', 'dielc_diff_mode', 'born_model', 'born_radius_model',
        'born_diff_mode', 'born_eps_mode', 'DH_model', 'bjeruum_treatment',
        'd_ion_mode', 'include_born_model', 'd_Born_mode',
        'born_solvation_shell_model', 'born_dielectric_saturation', 'born_bulk_mode',
        'mu_DH_diff_mode', 'mu_DH_comp_dep_rel_perm', 'mu_DH_include_sum_term',
        'mu_born_diff_mode', 'mu_born_comp_dep_rel_perm', 'mu_born_include_sum_term', 'mu_born_comp_dep_delta_d'
    }
    if any((k in params) for k in legacy_elec_keys):
        raise ValueError(
            'Flat electrostatic params are no longer supported; provide nested params["elec_model"] schema.'
        )

    elec_model = params.get('elec_model', None)
    if elec_model is None:
        if cppargs.z.size() > 0:
            # Backward-compatible electrolyte defaults when no explicit user options are provided.
            elec_model = {
                'rel_perm': {'rule': 1, 'differential_mode': 'analytical'},
                'hc_model': {'dadx_differential_mode': 'analytical'},
                'disp_model': {'dadx_differential_mode': 'analytical'},
                'assoc_model': {'dadx_differential_mode': 'analytical'},
                'polar_model': {'dadx_differential_mode': 'analytical'},
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
    diff_alias = {'analytic': 0, 'analytical': 0, 'numeric': 1, 'numerical': 1}
    d_ion_alias = {'t_indep': 0, 't_dep_1': 1, 't_dep_2': 2}
    d_born_alias = {'t_indep': 0, 't_dep_1': 1, 't_dep_2': 2, 'fitted_param': 3}
    bulk_alias = {'mix': 0, 'bulk': 0, 'solvent': 1}

    rel_perm = elec_model.get('rel_perm', {})
    if not isinstance(rel_perm, dict):
        raise ValueError('params["elec_model"]["rel_perm"] must be a dict.')
    hc_model_dict = elec_model.get('hc_model', {})
    if not isinstance(hc_model_dict, dict):
        raise ValueError('params["elec_model"]["hc_model"] must be a dict.')
    disp_model_dict = elec_model.get('disp_model', {})
    if not isinstance(disp_model_dict, dict):
        raise ValueError('params["elec_model"]["disp_model"] must be a dict.')
    assoc_model_dict = elec_model.get('assoc_model', {})
    if not isinstance(assoc_model_dict, dict):
        raise ValueError('params["elec_model"]["assoc_model"] must be a dict.')
    polar_model_dict = elec_model.get('polar_model', {})
    if not isinstance(polar_model_dict, dict):
        raise ValueError('params["elec_model"]["polar_model"] must be a dict.')
    dh_model_dict = elec_model.get('DH_model', {})
    if not isinstance(dh_model_dict, dict):
        raise ValueError('params["elec_model"]["DH_model"] must be a dict.')
    born_model_dict = elec_model.get('born_model', {})
    if not isinstance(born_model_dict, dict):
        raise ValueError('params["elec_model"]["born_model"] must be a dict.')
    mu_dh = dh_model_dict.get('mu_DH_model', {})
    if not isinstance(mu_dh, dict):
        raise ValueError('params["elec_model"]["DH_model"]["mu_DH_model"] must be a dict.')
    mu_born = born_model_dict.get('mu_born_model', {})
    if not isinstance(mu_born, dict):
        raise ValueError('params["elec_model"]["born_model"]["mu_born_model"] must be a dict.')

    cppargs.dielc_rule = _as_int_alias(rel_perm.get('rule', 1), rule_alias)
    cppargs.dielc_diff_mode = _as_int_alias(rel_perm.get('differential_mode', 'analytical'), diff_alias)
    if cppargs.dielc_diff_mode not in (0, 1):
        raise ValueError('Unknown rel_perm differential_mode. Supported values are analytical/numerical (0/1).')
    cppargs.hc_dadx_diff_mode = _as_int_alias(hc_model_dict.get('dadx_differential_mode', 'analytical'), diff_alias)
    if cppargs.hc_dadx_diff_mode not in (0, 1):
        raise ValueError('Unknown hc_model dadx_differential_mode. Supported values are analytical/numerical (0/1).')
    cppargs.disp_dadx_diff_mode = _as_int_alias(disp_model_dict.get('dadx_differential_mode', 'analytical'), diff_alias)
    if cppargs.disp_dadx_diff_mode not in (0, 1):
        raise ValueError('Unknown disp_model dadx_differential_mode. Supported values are analytical/numerical (0/1).')
    cppargs.assoc_dadx_diff_mode = _as_int_alias(assoc_model_dict.get('dadx_differential_mode', 'analytical'), diff_alias)
    if cppargs.assoc_dadx_diff_mode not in (0, 1):
        raise ValueError('Unknown assoc_model dadx_differential_mode. Supported values are analytical/numerical (0/1).')
    cppargs.polar_dadx_diff_mode = _as_int_alias(polar_model_dict.get('dadx_differential_mode', 'analytical'), diff_alias)
    if cppargs.polar_dadx_diff_mode not in (0, 1):
        raise ValueError('Unknown polar_model dadx_differential_mode. Supported values are analytical/numerical (0/1).')
    if cppargs.dielc_rule < 0 or cppargs.dielc_rule > 8:
        raise ValueError('Unknown rel_perm rule. Supported values are 0..8.')

    cppargs.d_ion_mode = _as_int_alias(dh_model_dict.get('d_ion_mode', 1), d_ion_alias)
    if cppargs.d_ion_mode not in (0, 1, 2):
        raise ValueError('Unknown d_ion_mode. Supported values are 0,1,2.')
    bjeruum = _as_bool(dh_model_dict.get('bjeruum_treatment', False))
    cppargs.mu_DH_diff_mode = _as_int_alias(mu_dh.get('differential_mode', 'analytical'), diff_alias)
    if cppargs.mu_DH_diff_mode not in (0, 1):
        raise ValueError('Unknown mu_DH differential_mode. Supported values are analytical/numerical (0/1).')
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

        # Legacy radius projection kept for C++ internals.
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
    cppargs.debug = int(bool(params['debug'])) if 'debug' in params else 0
    if 'assoc_num' in params:
        cppargs.assoc_num = np_to_vector_int(params['assoc_num'])
    if 'assoc_matrix' in params:
        cppargs.assoc_matrix = np_to_vector_int(params['assoc_matrix'])
    if 'k_hb' in params:
        cppargs.k_hb = np_to_vector_double(params['k_hb'])
    if 'l_ij' in params:
        cppargs.l_ij = np_to_vector_double(params['l_ij'])

    return cppargs




