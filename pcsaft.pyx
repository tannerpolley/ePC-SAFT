# -*- coding: utf-8 -*-
# setuptools: language=c++

import math
import numpy as np
from libcpp.vector cimport vector
from copy import deepcopy
cimport pcsaft

class InputError(Exception):
    # Exception raised for errors in the input.
    def __init__(self, message):
        self.message = message

class SolutionError(Exception):
    # Exception raised when a solver does not return a value.
    def __init__(self, message):
        self.message = message

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

def pcsaft_p(t, rho, x, params):
    """
    Calculate pressure.

    Parameters
    ----------
    t : float
        Temperature (K)
    rho : float
        Molar density (mol m\ :sup:`-3`)
    x : ndarray, shape (n,)
        Mole fractions of each component. It has a length of n, where n is
        the number of components in the system.
    params : dict
        A dictionary containing PC-SAFT parameters that can be passed for
        use in PC-SAFT:

        m : ndarray, shape (n,)
            Segment number for each component.
        s : ndarray, shape (n,)
            Segment diameter for each component. For ions this is the diameter of
            the hydrated ion. Units of Angstrom.
        e : ndarray, shape (n,)
            Dispersion energy of each component. For ions this is the dispersion
            energy of the hydrated ion. Units of K.
        k_ij : ndarray, shape (n,n)
            Binary interaction parameters between components in the mixture.
            (dimensions: ncomp x ncomp)
        e_assoc : ndarray, shape (n,)
            Association energy of the associating components. For non associating
            compounds this is set to 0. Units of K.
        vol_a : ndarray, shape (n,)
            Effective association volume of the associating components. For non
            associating compounds this is set to 0.
        dipm : ndarray, shape (n,)
            Dipole moment of the polar components. For components where the dipole
            term is not used this is set to 0. Units of Debye.
        dip_num : ndarray, shape (n,)
            The effective number of dipole functional groups on each component
            molecule. Generally this is set to 1, but some implementations use this
            as an adjustable parameter that is fit to data.
        z : ndarray, shape (n,)
            Charge number of the ions
        dielc : ndarray, shape (n,)
            Component dielectric constants used for electrolyte calculations.
        MW : ndarray, shape (n,)
            Molecular weights in kg/mol; required for non-water solvent
            molality conversion.
        osmotic_solvent_index : int, optional
            Explicit solvent index for the molality basis. If omitted, the
            method chooses a neutral component (z=0) automatically.
        assoc_scheme : list, shape (n,)
            The types of association sites for each component. Use `None` for molecules
            without association sites. If a molecule has multiple association sites,
            use a nested list for that component to specify the association scheme for
            each site. The accepted association schemes are those given by Huang and
            Radosz (1990): 1, 2A, 2B, 3A, 3B, 4A, 4B, 4C. If `e_assoc` and `vol_a` are
            given but `assoc_scheme` is not, the 2B association scheme is assumed (which
            would, for example, correspond to one hydroxyl functional group).

    Returns
    -------
    P : float
        Pressure (Pa)
    """
    x, params = ensure_numpy_input(x, params)
    check_input(x, {'density':rho, 'temperature':t})
    params = check_association(params)
    cppargs = create_struct(params)
    return pcsaft_p_cpp(t, rho, x, cppargs)


def pcsaft_lnfugcoef(t, rho, x, params):
    """
    Calculate the natural logarithm of the fugacity coefficients for one phase of the system.

    Parameters
    ----------
    t : float
        Temperature (K)
    rho : float
        Molar density (mol m\ :sup:`-3`)
    x : ndarray, shape (n,)
        Mole fractions of each component. It has a length of n, where n is
        the number of components in the system.
    params : dict
        A dictionary containing PC-SAFT parameters that can be passed for
        use in PC-SAFT:

        m : ndarray, shape (n,)
            Segment number for each component.
        s : ndarray, shape (n,)
            Segment diameter for each component. For ions this is the diameter of
            the hydrated ion. Units of Angstrom.
        e : ndarray, shape (n,)
            Dispersion energy of each component. For ions this is the dispersion
            energy of the hydrated ion. Units of K.
        k_ij : ndarray, shape (n,n)
            Binary interaction parameters between components in the mixture.
            (dimensions: ncomp x ncomp)
        e_assoc : ndarray, shape (n,)
            Association energy of the associating components. For non associating
            compounds this is set to 0. Units of K.
        vol_a : ndarray, shape (n,)
            Effective association volume of the associating components. For non
            associating compounds this is set to 0.
        dipm : ndarray, shape (n,)
            Dipole moment of the polar components. For components where the dipole
            term is not used this is set to 0. Units of Debye.
        dip_num : ndarray, shape (n,)
            The effective number of dipole functional groups on each component
            molecule. Generally this is set to 1, but some implementations use this
            as an adjustable parameter that is fit to data.
        z : ndarray, shape (n,)
            Charge number of the ions
        dielc : float
            Dielectric constant of the medium to be used for electrolyte
            calculations.
        assoc_scheme : list, shape (n,)
            The types of association sites for each component. Use `None` for molecules
            without association sites. If a molecule has multiple association sites,
            use a nested list for that component to specify the association scheme for
            each site. The accepted association schemes are those given by Huang and
            Radosz (1990): 1, 2A, 2B, 3A, 3B, 4A, 4B, 4C. If `e_assoc` and `vol_a` are
            given but `assoc_scheme` is not, the 2B association scheme is assumed (which
            would, for example, correspond to one hydroxyl functional group).

    Returns
    -------
    lnfugcoef : ndarray, shape (n,)
        Natural logarithm of the fugacity coefficients for each component.
    """
    x, params = ensure_numpy_input(x, params)
    check_input(x, {'density':rho, 'temperature':t})
    params = check_association(params)
    cppargs = create_struct(params)
    return np.asarray(pcsaft_lnfug_cpp(t, rho, x, cppargs))


def pcsaft_fugcoef(t, rho, x, params):
    """
    Calculate the fugacity coefficients for one phase of the system.

    Parameters
    ----------
    t : float
        Temperature (K)
    rho : float
        Molar density (mol m\ :sup:`-3`)
    x : ndarray, shape (n,)
        Mole fractions of each component. It has a length of n, where n is
        the number of components in the system.
    params : dict
        A dictionary containing PC-SAFT parameters that can be passed for
        use in PC-SAFT:

        m : ndarray, shape (n,)
            Segment number for each component.
        s : ndarray, shape (n,)
            Segment diameter for each component. For ions this is the diameter of
            the hydrated ion. Units of Angstrom.
        e : ndarray, shape (n,)
            Dispersion energy of each component. For ions this is the dispersion
            energy of the hydrated ion. Units of K.
        k_ij : ndarray, shape (n,n)
            Binary interaction parameters between components in the mixture.
            (dimensions: ncomp x ncomp)
        e_assoc : ndarray, shape (n,)
            Association energy of the associating components. For non associating
            compounds this is set to 0. Units of K.
        vol_a : ndarray, shape (n,)
            Effective association volume of the associating components. For non
            associating compounds this is set to 0.
        dipm : ndarray, shape (n,)
            Dipole moment of the polar components. For components where the dipole
            term is not used this is set to 0. Units of Debye.
        dip_num : ndarray, shape (n,)
            The effective number of dipole functional groups on each component
            molecule. Generally this is set to 1, but some implementations use this
            as an adjustable parameter that is fit to data.
        z : ndarray, shape (n,)
            Charge number of the ions
        dielc : float
            Dielectric constant of the medium to be used for electrolyte
            calculations.
        assoc_scheme : list, shape (n,)
            The types of association sites for each component. Use `None` for molecules
            without association sites. If a molecule has multiple association sites,
            use a nested list for that component to specify the association scheme for
            each site. The accepted association schemes are those given by Huang and
            Radosz (1990): 1, 2A, 2B, 3A, 3B, 4A, 4B, 4C. If `e_assoc` and `vol_a` are
            given but `assoc_scheme` is not, the 2B association scheme is assumed (which
            would, for example, correspond to one hydroxyl functional group).

    Returns
    -------
    fugcoef : ndarray, shape (n,)
        Fugacity coefficients of each component.
    """
    x, params = ensure_numpy_input(x, params)
    check_input(x, {'density':rho, 'temperature':t})
    params = check_association(params)
    cppargs = create_struct(params)
    return np.asarray(pcsaft_fugcoef_cpp(t, rho, x, cppargs))


def pcsaft_Z(t, rho, x, params):
    """
    Calculate the compressibility factor.

    Parameters
    ----------
    t : float
        Temperature (K)
    rho : float
        Molar density (mol m\ :sup:`-3`)
    x : ndarray, shape (n,)
        Mole fractions of each component. It has a length of n, where n is
        the number of components in the system.
    params : dict
        A dictionary containing PC-SAFT parameters that can be passed for
        use in PC-SAFT:

        m : ndarray, shape (n,)
            Segment number for each component.
        s : ndarray, shape (n,)
            Segment diameter for each component. For ions this is the diameter of
            the hydrated ion. Units of Angstrom.
        e : ndarray, shape (n,)
            Dispersion energy of each component. For ions this is the dispersion
            energy of the hydrated ion. Units of K.
        k_ij : ndarray, shape (n,n)
            Binary interaction parameters between components in the mixture.
            (dimensions: ncomp x ncomp)
        e_assoc : ndarray, shape (n,)
            Association energy of the associating components. For non associating
            compounds this is set to 0. Units of K.
        vol_a : ndarray, shape (n,)
            Effective association volume of the associating components. For non
            associating compounds this is set to 0.
        dipm : ndarray, shape (n,)
            Dipole moment of the polar components. For components where the dipole
            term is not used this is set to 0. Units of Debye.
        dip_num : ndarray, shape (n,)
            The effective number of dipole functional groups on each component
            molecule. Generally this is set to 1, but some implementations use this
            as an adjustable parameter that is fit to data.
        z : ndarray, shape (n,)
            Charge number of the ions
        dielc : float
            Dielectric constant of the medium to be used for electrolyte
            calculations.
        assoc_scheme : list, shape (n,)
            The types of association sites for each component. Use `None` for molecules
            without association sites. If a molecule has multiple association sites,
            use a nested list for that component to specify the association scheme for
            each site. The accepted association schemes are those given by Huang and
            Radosz (1990): 1, 2A, 2B, 3A, 3B, 4A, 4B, 4C. If `e_assoc` and `vol_a` are
            given but `assoc_scheme` is not, the 2B association scheme is assumed (which
            would, for example, correspond to one hydroxyl functional group).

    Returns
    -------
    Z : float
        Compressibility factor
    """
    x, params = ensure_numpy_input(x, params)
    check_input(x, {'density':rho, 'temperature':t})
    params = check_association(params)
    cppargs = create_struct(params)
    return pcsaft_Z_cpp(t, rho, x, cppargs)


def flashPQ(p, q, x, params, t_guess=None):
    """
    Calculate the temperature of the system where vapor and liquid phases are in equilibrium.

    Parameters
    ----------
    p : float
        Pressure (Pa)
    q : float
        Mole fraction of the fluid in the vapor phase
    x : ndarray, shape (n,)
        Mole fractions of each component. It has a length of n, where n is
        the number of components in the system.
    params : dict
        A dictionary containing PC-SAFT parameters that can be passed for
        use in PC-SAFT:

        m : ndarray, shape (n,)
            Segment number for each component.
        s : ndarray, shape (n,)
            Segment diameter for each component. For ions this is the diameter of
            the hydrated ion. Units of Angstrom.
        e : ndarray, shape (n,)
            Dispersion energy of each component. For ions this is the dispersion
            energy of the hydrated ion. Units of K.
        k_ij : ndarray, shape (n,n)
            Binary interaction parameters between components in the mixture.
            (dimensions: ncomp x ncomp)
        e_assoc : ndarray, shape (n,)
            Association energy of the associating components. For non associating
            compounds this is set to 0. Units of K.
        vol_a : ndarray, shape (n,)
            Effective association volume of the associating components. For non
            associating compounds this is set to 0.
        dipm : ndarray, shape (n,)
            Dipole moment of the polar components. For components where the dipole
            term is not used this is set to 0. Units of Debye.
        dip_num : ndarray, shape (n,)
            The effective number of dipole functional groups on each component
            molecule. Generally this is set to 1, but some implementations use this
            as an adjustable parameter that is fit to data.
        z : ndarray, shape (n,)
            Charge number of the ions
        dielc : float
            Dielectric constant of the medium to be used for electrolyte
            calculations.
        assoc_scheme : list, shape (n,)
            The types of association sites for each component. Use `None` for molecules
            without association sites. If a molecule has multiple association sites,
            use a nested list for that component to specify the association scheme for
            each site. The accepted association schemes are those given by Huang and
            Radosz (1990): 1, 2A, 2B, 3A, 3B, 4A, 4B, 4C. If `e_assoc` and `vol_a` are
            given but `assoc_scheme` is not, the 2B association scheme is assumed (which
            would, for example, correspond to one hydroxyl functional group).

    t_guess : float
        Initial guess for the temperature (K) (optional)

    Returns
    -------
    t : float
        Temperature (K)
    xl : ndarray, shape (n,)
        Liquid mole fractions after flash
    xv : ndarray, shape (n,)
        Vapor mole fractions after flash

    Notes
    -----
    To solve the PQ flash the temperature must be varied. This adds additional complexity
    for water and electrolyte mixtures. For water, a temperature dependent sigma is often
    used. However, there does not appear to be a way to pass a Python function to the C++
    code without requiring the user to compile it using Cython. To avoid this, the `flashPQ`
    function uses the following relationship internally to calculate sigma for water as a
    function of temperature: ::

        3.8395 + 1.2828 * exp(-0.0074944 * t) - 1.3939 * exp(-0.00056029 * t);

    For electrolyte solutions the dielectric constant is calculated using the `dielc_water`
    function. This means that the sigma value for water and the dielectric constant given by
    the user are not used by the `flashPQ` function.

    The code identifies which component is water by the epsilon/k value. Therefore, when
    using `flashPQ` with water `e` must be exactly 353.9449, if you want the temperature
    dependence of sigma to be accounted for.

    If you want to use different functions for temperature dependent parameters with `flashPQ`
    then you will need to modify the source code and recompile it.
    """
    x, params = ensure_numpy_input(x, params)
    check_input(x, {'pressure':p, 'Q':q})
    params = check_association(params)
    cppargs = create_struct(params)
    try:
        if t_guess is not None:
            result = flashPQ_cpp(p, q, x, cppargs, t_guess)
        else:
            result = flashPQ_cpp(p, q, x, cppargs)
    except:
        raise SolutionError('A solution was not found for flashPQ. P={}'.format(p))

    t = result[0]
    xl = np.asarray(result[1:])
    xl, xv = np.split(xl, 2)
    return t, xl, xv


def flashTQ(t, q, x, params, p_guess=None):
    """
    Calculate the pressure of the system where vapor and liquid phases are in equilibrium.

    Parameters
    ----------
    t : float
        Temperature (K)
    q : float
        Mole fraction of the fluid in the vapor phase
    x : ndarray, shape (n,)
        Mole fractions of each component. It has a length of n, where n is
        the number of components in the system.
    params : dict
        A dictionary containing PC-SAFT parameters that can be passed for
        use in PC-SAFT:

        m : ndarray, shape (n,)
            Segment number for each component.
        s : ndarray, shape (n,)
            Segment diameter for each component. For ions this is the diameter of
            the hydrated ion. Units of Angstrom.
        e : ndarray, shape (n,)
            Dispersion energy of each component. For ions this is the dispersion
            energy of the hydrated ion. Units of K.
        k_ij : ndarray, shape (n,n)
            Binary interaction parameters between components in the mixture.
            (dimensions: ncomp x ncomp)
        e_assoc : ndarray, shape (n,)
            Association energy of the associating components. For non associating
            compounds this is set to 0. Units of K.
        vol_a : ndarray, shape (n,)
            Effective association volume of the associating components. For non
            associating compounds this is set to 0.
        dipm : ndarray, shape (n,)
            Dipole moment of the polar components. For components where the dipole
            term is not used this is set to 0. Units of Debye.
        dip_num : ndarray, shape (n,)
            The effective number of dipole functional groups on each component
            molecule. Generally this is set to 1, but some implementations use this
            as an adjustable parameter that is fit to data.
        z : ndarray, shape (n,)
            Charge number of the ions
        dielc : float
            Dielectric constant of the medium to be used for electrolyte
            calculations.
        assoc_scheme : list, shape (n,)
            The types of association sites for each component. Use `None` for molecules
            without association sites. If a molecule has multiple association sites,
            use a nested list for that component to specify the association scheme for
            each site. The accepted association schemes are those given by Huang and
            Radosz (1990): 1, 2A, 2B, 3A, 3B, 4A, 4B, 4C. If `e_assoc` and `vol_a` are
            given but `assoc_scheme` is not, the 2B association scheme is assumed (which
            would, for example, correspond to one hydroxyl functional group).

    p_guess : float
        Initial guess for the pressure (Pa) (optional)

    Returns
    -------
    p : float
        Pressure (Pa)
    xl : ndarray, shape (n,)
        Liquid mole fractions after flash
    xv : ndarray, shape (n,)
        Vapor mole fractions after flash
    """
    x, params = ensure_numpy_input(x, params)
    check_input(x, {'temperature':t, 'Q':q})
    params = check_association(params)
    cppargs = create_struct(params)
    try:
        if p_guess is not None:
            result = flashTQ_cpp(t, q, x, cppargs, p_guess)
        else:
            result = flashTQ_cpp(t, q, x, cppargs)
    except:
        raise SolutionError('A solution was not found for flashTQ. T={}'.format(t))

    p = result[0]
    xl = np.asarray(result[1:])
    xl, xv = np.split(xl, 2)
    return p, xl, xv

def pcsaft_Hvap(t, x, params, p_guess=None):
    """
    Calculate the enthalpy of vaporization.

    Parameters
    ----------
    t : float
        Temperature (K)
    x : ndarray, shape (n,)
        Mole fractions of each component. It has a length of n, where n is
        the number of components in the system.
    params : dict
        A dictionary containing PC-SAFT parameters that can be passed for
        use in PC-SAFT:

        m : ndarray, shape (n,)
            Segment number for each component.
        s : ndarray, shape (n,)
            Segment diameter for each component. For ions this is the diameter of
            the hydrated ion. Units of Angstrom.
        e : ndarray, shape (n,)
            Dispersion energy of each component. For ions this is the dispersion
            energy of the hydrated ion. Units of K.
        k_ij : ndarray, shape (n,n)
            Binary interaction parameters between components in the mixture.
            (dimensions: ncomp x ncomp)
        e_assoc : ndarray, shape (n,)
            Association energy of the associating components. For non associating
            compounds this is set to 0. Units of K.
        vol_a : ndarray, shape (n,)
            Effective association volume of the associating components. For non
            associating compounds this is set to 0.
        dipm : ndarray, shape (n,)
            Dipole moment of the polar components. For components where the dipole
            term is not used this is set to 0. Units of Debye.
        dip_num : ndarray, shape (n,)
            The effective number of dipole functional groups on each component
            molecule. Generally this is set to 1, but some implementations use this
            as an adjustable parameter that is fit to data.
        z : ndarray, shape (n,)
            Charge number of the ions
        dielc : float
            Dielectric constant of the medium to be used for electrolyte
            calculations.
        assoc_scheme : list, shape (n,)
            The types of association sites for each component. Use `None` for molecules
            without association sites. If a molecule has multiple association sites,
            use a nested list for that component to specify the association scheme for
            each site. The accepted association schemes are those given by Huang and
            Radosz (1990): 1, 2A, 2B, 3A, 3B, 4A, 4B, 4C. If `e_assoc` and `vol_a` are
            given but `assoc_scheme` is not, the 2B association scheme is assumed (which
            would, for example, correspond to one hydroxyl functional group).

    p_guess : float
        Guess for the vapor pressure (Pa) (optional)

    Returns
    -------
    output : list
        A list containing the following results:
            0 : enthalpy of vaporization (J/mol), float
            1 : vapor pressure (Pa), float
    """
    x, params = ensure_numpy_input(x, params)
    check_input(x, {'temperature': t})
    params = check_association(params)
    cppargs = create_struct(params)

    q = 0
    try:
        if p_guess is not None:
            result = np.asarray(flashTQ_cpp(t, q, x, cppargs, p_guess))
            Pvap = result[0]
        else:
            result = np.asarray(flashTQ_cpp(t, q, x, cppargs))
            Pvap = result[0]
    except:
        raise SolutionError('A solution was not found for flashTQ. T={}'.format(t))

    rho = pcsaft_den_cpp(t, Pvap, x, 0, cppargs)
    hres_l = pcsaft_hres_cpp(t, rho, x, cppargs)
    rho = pcsaft_den_cpp(t, Pvap, x, 1, cppargs)
    hres_v = pcsaft_hres_cpp(t, rho, x, cppargs)
    Hvap = hres_v - hres_l

    output = [Hvap, Pvap]
    return output


def pcsaft_osmoticC(t, rho, x, params):
    """
    Calculate the osmotic coefficient.

    Parameters
    ----------
    t : float
        Temperature (K)
    rho : float
        Molar density (mol m\ :sup:`-3`)
    x : ndarray, shape (n,)
        Mole fractions of each component. It has a length of n, where n is
        the number of components in the system.
    params : dict
        A dictionary containing PC-SAFT parameters that can be passed for
        use in PC-SAFT:

        m : ndarray, shape (n,)
            Segment number for each component.
        s : ndarray, shape (n,)
            Segment diameter for each component. For ions this is the diameter of
            the hydrated ion. Units of Angstrom.
        e : ndarray, shape (n,)
            Dispersion energy of each component. For ions this is the dispersion
            energy of the hydrated ion. Units of K.
        k_ij : ndarray, shape (n,n)
            Binary interaction parameters between components in the mixture.
            (dimensions: ncomp x ncomp)
        e_assoc : ndarray, shape (n,)
            Association energy of the associating components. For non associating
            compounds this is set to 0. Units of K.
        vol_a : ndarray, shape (n,)
            Effective association volume of the associating components. For non
            associating compounds this is set to 0.
        dipm : ndarray, shape (n,)
            Dipole moment of the polar components. For components where the dipole
            term is not used this is set to 0. Units of Debye.
        dip_num : ndarray, shape (n,)
            The effective number of dipole functional groups on each component
            molecule. Generally this is set to 1, but some implementations use this
            as an adjustable parameter that is fit to data.
        z : ndarray, shape (n,)
            Charge number of the ions
        dielc : float
            Dielectric constant of the medium to be used for electrolyte
            calculations.
        assoc_scheme : list, shape (n,)
            The types of association sites for each component. Use `None` for molecules
            without association sites. If a molecule has multiple association sites,
            use a nested list for that component to specify the association scheme for
            each site. The accepted association schemes are those given by Huang and
            Radosz (1990): 1, 2A, 2B, 3A, 3B, 4A, 4B, 4C. If `e_assoc` and `vol_a` are
            given but `assoc_scheme` is not, the 2B association scheme is assumed (which
            would, for example, correspond to one hydroxyl functional group).

    Returns
    -------
    osmC : float
        Molal osmotic coefficient
    """
    x, params = ensure_numpy_input(x, params)
    check_input(x, {'density':rho, 'temperature':t})
    params = check_association(params)
    cppargs = create_struct(params)

    # Select solvent component for molality basis.
    if 'osmotic_solvent_index' in params:
        indx_solvent = int(params['osmotic_solvent_index'])
        if indx_solvent < 0 or indx_solvent >= len(x):
            raise InputError('params["osmotic_solvent_index"] is out of bounds.')
    else:
        indx_solvent = -1
        z = np.asarray(params.get('z', []), dtype=float).flatten()
        if z.size == len(x):
            idx_sol = np.where(np.abs(z) < 1e-12)[0]
            if idx_sol.size == 1:
                indx_solvent = int(idx_sol[0])
            elif idx_sol.size > 1:
                # For mixed neutral solvents, default to the dominant neutral component.
                indx_solvent = int(idx_sol[np.argmax(x[idx_sol])])
        # Backward-compatible fallback for legacy water-only parameter sets.
        if indx_solvent < 0 and 'e' in params:
            idx_water = np.where(np.isclose(np.asarray(params['e'], dtype=float), 353.9449, atol=1e-6))[0]
            if idx_water.size == 1:
                indx_solvent = int(idx_water[0])
        if indx_solvent < 0:
            raise InputError('pcsaft_osmoticC requires a solvent component. Provide params["osmotic_solvent_index"] or one neutral species (z=0).')

    mw = np.asarray(params.get('MW', []), dtype=float).flatten()
    if mw.size == len(x):
        mw_solvent = float(mw[indx_solvent])
    elif 'e' in params and np.isclose(np.asarray(params['e'], dtype=float)[indx_solvent], 353.9449, atol=1e-6):
        # Legacy fallback for water if MW is omitted.
        mw_solvent = 18.0153/1000.
    else:
        raise InputError('pcsaft_osmoticC requires params["MW"] to compute molality for non-water solvent.')
    # If a molar mass in g/mol was passed by mistake, convert to kg/mol.
    if mw_solvent > 1.0:
        mw_solvent = mw_solvent/1000.0
    if mw_solvent <= 0.0:
        raise InputError('Solvent molecular weight must be positive.')
    if x[indx_solvent] <= 0.0:
        raise InputError('Solvent mole fraction must be positive.')

    molality = x/(x[indx_solvent]*mw_solvent)
    molality[indx_solvent] = 0
    molality_sum = float(np.sum(molality))
    if molality_sum <= 0.0:
        raise InputError('Total molality is zero; osmotic coefficient is undefined.')
    x0 = np.zeros_like(x)
    x0[indx_solvent] = 1.

    fugcoef = np.asarray(pcsaft_fugcoef_cpp(t, rho, x, cppargs))
    p = pcsaft_p_cpp(t, rho, x, cppargs)
    if rho < 900:
        ph = 1
    else:
        ph = 0
    rho0 = pcsaft_den_cpp(t, p, x0, ph, cppargs)
    fugcoef0 = np.asarray(pcsaft_fugcoef_cpp(t, rho0, x0, cppargs))
    gamma = float(fugcoef[indx_solvent]/fugcoef0[indx_solvent])

    osmC = -np.log(x[indx_solvent]*gamma)/(mw_solvent*molality_sum)
    # Keep legacy return shape for callers/tests that index result[0].
    return np.asarray([osmC], dtype=float)

def pcsaft_miac_m(t, rho, x, params, species=None):
    """
    Molality-scale mean ionic activity coefficient (MIAC_m).
    """
    x, params = ensure_numpy_input(x, params)
    check_input(x, {'density':rho, 'temperature':t})
    params = check_association(params)
    cppargs = create_struct(params)

    z = np.asarray(params.get('z', []), dtype=float)
    if z.size == 0 or np.allclose(z, 0):
        raise InputError('pcsaft_miac_m requires ionic species (non-zero z).')
    idx_cat = np.where(z > 0)[0]
    idx_an = np.where(z < 0)[0]
    idx_sol = np.where(np.abs(z) < 1e-12)[0]
    if len(idx_cat) == 0 or len(idx_an) == 0:
        raise InputError('pcsaft_miac_m needs at least one cation and one anion.')
    if len(idx_sol) == 0:
        raise InputError('pcsaft_miac_m needs a neutral solvent to define molality.')

    mw = np.asarray(params['MW'], dtype=float)
    mass_solvent = float(np.sum(x[idx_sol] * mw[idx_sol]))
    if mass_solvent <= 0:
        raise InputError('Solvent mass is zero; check solvent mole fraction and MW.')

    fugcoef = np.asarray(pcsaft_fugcoef_cpp(t, rho, x, cppargs), dtype=float)

    eps = 1e-12
    x_inf = np.full_like(x, eps)
    x_inf[idx_sol[0]] = max(1.0 - eps * (len(x) - 1), eps)
    x_inf /= np.sum(x_inf)
    rho_inf = pcsaft_den_cpp(t, pcsaft_p_cpp(t, rho, x, cppargs), x_inf, 0, cppargs)
    fugcoef_inf = np.asarray(pcsaft_fugcoef_cpp(t, rho_inf, x_inf, cppargs), dtype=float)
    if np.any(fugcoef_inf <= 0):
        raise SolutionError('Non-positive fugacity at infinite dilution.')
    gamma_i = fugcoef / fugcoef_inf

    mass_neutral = x[idx_sol] * mw[idx_sol]
    w_sf = mass_neutral / mass_neutral.sum()
    M_solvent_mix = 1.0 / np.sum(w_sf / mw[idx_sol])

    if species is None or len(species) != len(x):
        raise InputError('species list (matching x order) is required to label salts.')

    result = {}
    for ic in idx_cat:
        for ia in idx_an:
            zc = int(round(abs(z[ic])))
            za = int(round(abs(z[ia])))
            g = math.gcd(zc, za)
            nu_cat = za // g
            nu_an = zc // g
            n_salt = 0.5 * (x[ic] / nu_cat + x[ia] / nu_an)
            m_salt = n_salt / mass_solvent
            sum_nu = float(nu_cat + nu_an)
            ln_gamma_pm = (nu_cat * math.log(gamma_i[ic]) + nu_an * math.log(gamma_i[ia])) / sum_nu
            gamma_pm_x = math.exp(ln_gamma_pm)
            gamma_pm_m = gamma_pm_x / (1.0 + M_solvent_mix * m_salt * sum_nu)
            salt_name = species[ic] + species[ia]
            result[salt_name] = gamma_pm_m

    return result

def pcsaft_cp(t, rho, aly_lee_params, x, params):
    """
    Calculate the specific molar isobaric heat capacity.

    Parameters
    ----------
    t : float
        Temperature (K)
    rho : float
        Molar density (mol m\ :sup:`-3`)
    aly_lee_params : ndarray, shape (5,)
        Constants for the Aly-Lee equation. Can be substituted with parameters for
        another equation if the ideal gas heat capacity is given using a different
        equation.
    x : ndarray, shape (n,)
        Mole fractions of each component. It has a length of n, where n is
        the number of components in the system.
    params : dict
        A dictionary containing PC-SAFT parameters that can be passed for
        use in PC-SAFT:

        m : ndarray, shape (n,)
            Segment number for each component.
        s : ndarray, shape (n,)
            Segment diameter for each compopynent. For ions this is the diameter of
            the hydrated ion. Units of Angstrom.
        e : ndarray, shape (n,)
            Dispersion energy of each component. For ions this is the dispersion
            energy of the hydrated ion. Units of K.
        k_ij : ndarray, shape (n,n)
            Binary interaction parameters between components in the mixture.
            (dimensions: ncomp x ncomp)
        e_assoc : ndarray, shape (n,)
            Association energy of the associating components. For non associating
            compounds this is set to 0. Units of K.
        vol_a : ndarray, shape (n,)
            Effective association volume of the associating components. For non
            associating compounds this is set to 0.
        dipm : ndarray, shape (n,)
            Dipole moment of the polar components. For components where the dipole
            term is not used this is set to 0. Units of Debye.
        dip_num : ndarray, shape (n,)
            The effective number of dipole functional groups on each component
            molecule. Generally this is set to 1, but some implementations use this
            as an adjustable parameter that is fit to data.
        z : ndarray, shape (n,)
            Charge number of the ions
        dielc : float
            Dielectric constant of the medium to be used for electrolyte
            calculations.
        assoc_scheme : list, shape (n,)
            The types of association sites for each component. Use `None` for molecules
            without association sites. If a molecule has multiple association sites,
            use a nested list for that component to specify the association scheme for
            each site. The accepted association schemes are those given by Huang and
            Radosz (1990): 1, 2A, 2B, 3A, 3B, 4A, 4B, 4C. If `e_assoc` and `vol_a` are
            given but `assoc_scheme` is not, the 2B association scheme is assumed (which
            would, for example, correspond to one hydroxyl functional group).

    Returns
    -------
    cp : float
        Specific molar isobaric heat capacity (J mol\ :sup:`-1` K\ :sup:`-1`)
    """
    x, params = ensure_numpy_input(x, params)
    check_input(x, {'density':rho, 'temperature':t})
    params = check_association(params)

    if rho > 900:
        ph = 0
    else:
        ph = 1

    cppargs = create_struct(params)

    cp_ideal = aly_lee(t, aly_lee_params)
    p = pcsaft_p_cpp(t, rho, x, cppargs)
    rho0 = pcsaft_den_cpp(t-0.001, p, x, ph, cppargs)
    hres0 = pcsaft_hres_cpp(t-0.001, rho0, x, cppargs)
    rho1 = pcsaft_den_cpp(t+0.001, p, x, ph, cppargs)
    hres1 = pcsaft_hres_cpp(t+0.001, rho1, x, cppargs)
    dhdt = (hres1-hres0)/0.002 # a numerical derivative is used for now until analytical derivatives are ready
    return cp_ideal + dhdt


def pcsaft_den(t, p, x, params, phase='liq'):
    """
    Calculate the molar density.

    Parameters
    ----------
    t : float
        Temperature (K)
    p : float
        Pressure (Pa)
    x : ndarray, shape (n,)
        Mole fractions of each component. It has a length of n, where n is
        the number of components in the system.
    params : dict
        A dictionary containing PC-SAFT parameters that can be passed for
        use in PC-SAFT:

        m : ndarray, shape (n,)
            Segment number for each component.
        s : ndarray, shape (n,)
            Segment diameter for each component. For ions this is the diameter of
            the hydrated ion. Units of Angstrom.
        e : ndarray, shape (n,)
            Dispersion energy of each component. For ions this is the dispersion
            energy of the hydrated ion. Units of K.
        k_ij : ndarray, shape (n,n)
            Binary interaction parameters between components in the mixture.
            (dimensions: ncomp x ncomp)
        e_assoc : ndarray, shape (n,)
            Association energy of the associating components. For non associating
            compounds this is set to 0. Units of K.
        vol_a : ndarray, shape (n,)
            Effective association volume of the associating components. For non
            associating compounds this is set to 0.
        dipm : ndarray, shape (n,)
            Dipole moment of the polar components. For components where the dipole
            term is not used this is set to 0. Units of Debye.
        dip_num : ndarray, shape (n,)
            The effective number of dipole functional groups on each component
            molecule. Generally this is set to 1, but some implementations use this
            as an adjustable parameter that is fit to data.
        z : ndarray, shape (n,)
            Charge number of the ions
        dielc : float
            Dielectric constant of the medium to be used for electrolyte
            calculations.
        assoc_scheme : list, shape (n,)
            The types of association sites for each component. Use `None` for molecules
            without association sites. If a molecule has multiple association sites,
            use a nested list for that component to specify the association scheme for
            each site. The accepted association schemes are those given by Huang and
            Radosz (1990): 1, 2A, 2B, 3A, 3B, 4A, 4B, 4C. If `e_assoc` and `vol_a` are
            given but `assoc_scheme` is not, the 2B association scheme is assumed (which
            would, for example, correspond to one hydroxyl functional group).

    phase : string
        The phase for which the calculation is performed. Options: "liq" (liquid),
        "vap" (vapor).

    Returns
    -------
    rho : float
        Molar density (mol m\ :sup:`-3`)
    """
    x, params = ensure_numpy_input(x, params)
    check_input(x, {'pressure':p, 'temperature':t})
    params = check_association(params)
    cppargs = create_struct(params)
    if phase == 'liq':
        phase_num = 0
    else:
        phase_num = 1

    return pcsaft_den_cpp(t, p, x, phase_num, cppargs)


def pcsaft_hres(t, rho, x, params):
    """
    Calculate the residual enthalpy for one phase of the system.

    Parameters
    ----------
    t : float
        Temperature (K)
    rho : float
        Molar density (mol m\ :sup:`-3`)
    x : ndarray, shape (n,)
        Mole fractions of each component. It has a length of n, where n is
        the number of components in the system.
    params : dict
        A dictionary containing PC-SAFT parameters that can be passed for
        use in PC-SAFT:

        m : ndarray, shape (n,)
            Segment number for each component.
        s : ndarray, shape (n,)
            Segment diameter for each component. For ions this is the diameter of
            the hydrated ion. Units of Angstrom.
        e : ndarray, shape (n,)
            Dispersion energy of each component. For ions this is the dispersion
            energy of the hydrated ion. Units of K.
        k_ij : ndarray, shape (n,n)
            Binary interaction parameters between components in the mixture.
            (dimensions: ncomp x ncomp)
        e_assoc : ndarray, shape (n,)
            Association energy of the associating components. For non associating
            compounds this is set to 0. Units of K.
        vol_a : ndarray, shape (n,)
            Effective association volume of the associating components. For non
            associating compounds this is set to 0.
        dipm : ndarray, shape (n,)
            Dipole moment of the polar components. For components where the dipole
            term is not used this is set to 0. Units of Debye.
        dip_num : ndarray, shape (n,)
            The effective number of dipole functional groups on each component
            molecule. Generally this is set to 1, but some implementations use this
            as an adjustable parameter that is fit to data.
        z : ndarray, shape (n,)
            Charge number of the ions
        dielc : float
            Dielectric constant of the medium to be used for electrolyte
            calculations.
        assoc_scheme : list, shape (n,)
            The types of association sites for each component. Use `None` for molecules
            without association sites. If a molecule has multiple association sites,
            use a nested list for that component to specify the association scheme for
            each site. The accepted association schemes are those given by Huang and
            Radosz (1990): 1, 2A, 2B, 3A, 3B, 4A, 4B, 4C. If `e_assoc` and `vol_a` are
            given but `assoc_scheme` is not, the 2B association scheme is assumed (which
            would, for example, correspond to one hydroxyl functional group).

    Returns
    -------
    hres : float
        Residual enthalpy (J mol\ :sup:`-1`)
    """
    x, params = ensure_numpy_input(x, params)
    check_input(x, {'density':rho, 'temperature':t})
    params = check_association(params)
    cppargs = create_struct(params)
    return pcsaft_hres_cpp(t, rho, x, cppargs)

def pcsaft_sres(t, rho, x, params):
    """
    Calculate the residual entropy (constant volume) for one phase of the system.

    Parameters
    ----------
    t : float
        Temperature (K)
    rho : float
        Molar density (mol m\ :sup:`-3`)
    x : ndarray, shape (n,)
        Mole fractions of each component. It has a length of n, where n is
        the number of components in the system.
    params : dict
        A dictionary containing PC-SAFT parameters that can be passed for
        use in PC-SAFT:

        m : ndarray, shape (n,)
            Segment number for each component.
        s : ndarray, shape (n,)
            Segment diameter for each component. For ions this is the diameter of
            the hydrated ion. Units of Angstrom.
        e : ndarray, shape (n,)
            Dispersion energy of each component. For ions this is the dispersion
            energy of the hydrated ion. Units of K.
        k_ij : ndarray, shape (n,n)
            Binary interaction parameters between components in the mixture.
            (dimensions: ncomp x ncomp)
        e_assoc : ndarray, shape (n,)
            Association energy of the associating components. For non associating
            compounds this is set to 0. Units of K.
        vol_a : ndarray, shape (n,)
            Effective association volume of the associating components. For non
            associating compounds this is set to 0.
        dipm : ndarray, shape (n,)
            Dipole moment of the polar components. For components where the dipole
            term is not used this is set to 0. Units of Debye.
        dip_num : ndarray, shape (n,)
            The effective number of dipole functional groups on each component
            molecule. Generally this is set to 1, but some implementations use this
            as an adjustable parameter that is fit to data.
        z : ndarray, shape (n,)
            Charge number of the ions
        dielc : float
            Dielectric constant of the medium to be used for electrolyte
            calculations.
        assoc_scheme : list, shape (n,)
            The types of association sites for each component. Use `None` for molecules
            without association sites. If a molecule has multiple association sites,
            use a nested list for that component to specify the association scheme for
            each site. The accepted association schemes are those given by Huang and
            Radosz (1990): 1, 2A, 2B, 3A, 3B, 4A, 4B, 4C. If `e_assoc` and `vol_a` are
            given but `assoc_scheme` is not, the 2B association scheme is assumed (which
            would, for example, correspond to one hydroxyl functional group).

    Returns
    -------
    sres : float
        Residual entropy (J mol\ :sup:`-1` K\ :sup:`-1`)
    """
    x, params = ensure_numpy_input(x, params)
    check_input(x, {'density':rho, 'temperature':t})
    params = check_association(params)
    cppargs = create_struct(params)
    return pcsaft_sres_cpp(t, rho, x, cppargs)

def pcsaft_gres(t, rho, x, params):
    """
    Calculate the residual Gibbs energy for one phase of the system.

    Parameters
    ----------
    t : float
        Temperature (K)
    rho : float
        Molar density (mol m\ :sup:`-3`)
    x : ndarray, shape (n,)
        Mole fractions of each component. It has a length of n, where n is
        the number of components in the system.
    params : dict
        A dictionary containing PC-SAFT parameters that can be passed for
        use in PC-SAFT:

        m : ndarray, shape (n,)
            Segment number for each component.
        s : ndarray, shape (n,)
            Segment diameter for each component. For ions this is the diameter of
            the hydrated ion. Units of Angstrom.
        e : ndarray, shape (n,)
            Dispersion energy of each component. For ions this is the dispersion
            energy of the hydrated ion. Units of K.
        k_ij : ndarray, shape (n,n)
            Binary interaction parameters between components in the mixture.
            (dimensions: ncomp x ncomp)
        e_assoc : ndarray, shape (n,)
            Association energy of the associating components. For non associating
            compounds this is set to 0. Units of K.
        vol_a : ndarray, shape (n,)
            Effective association volume of the associating components. For non
            associating compounds this is set to 0.
        dipm : ndarray, shape (n,)
            Dipole moment of the polar components. For components where the dipole
            term is not used this is set to 0. Units of Debye.
        dip_num : ndarray, shape (n,)
            The effective number of dipole functional groups on each component
            molecule. Generally this is set to 1, but some implementations use this
            as an adjustable parameter that is fit to data.
        z : ndarray, shape (n,)
            Charge number of the ions
        dielc : float
            Dielectric constant of the medium to be used for electrolyte
            calculations.
        assoc_scheme : list, shape (n,)
            The types of association sites for each component. Use `None` for molecules
            without association sites. If a molecule has multiple association sites,
            use a nested list for that component to specify the association scheme for
            each site. The accepted association schemes are those given by Huang and
            Radosz (1990): 1, 2A, 2B, 3A, 3B, 4A, 4B, 4C. If `e_assoc` and `vol_a` are
            given but `assoc_scheme` is not, the 2B association scheme is assumed (which
            would, for example, correspond to one hydroxyl functional group).

    Returns
    -------
    gres : float
        Residual Gibbs energy (J mol\ :sup:`-1`)
    """
    x, params = ensure_numpy_input(x, params)
    check_input(x, {'density':rho, 'temperature':t})
    params = check_association(params)
    cppargs = create_struct(params)
    return pcsaft_gres_cpp(t, rho, x, cppargs)

def pcsaft_gsolv(t, rho, x, params, species=None):
    """
    Gibbs solvation energy at infinite dilution on the mole-fraction scale for ions.
    """
    x, params = ensure_numpy_input(x, params)
    check_input(x, {'density':rho, 'temperature':t})
    params = check_association(params)
    cppargs = create_struct(params)

    z = np.asarray(params.get('z', []), dtype=float)
    idx_ion = np.where(np.abs(z) > 1e-12)[0]
    if len(idx_ion) == 0:
        raise InputError('pcsaft_gsolv requires ionic species in params["z"].')
    if species is None or len(species) != len(x):
        raise InputError('species list (matching x order) is required to label ions.')

    idx_solv = np.where(np.abs(z) <= 1e-12)[0]
    if len(idx_solv) == 0:
        raise InputError('pcsaft_gsolv requires at least one solvent species (z=0).')

    x_ref = np.asarray(x, dtype=float).copy()
    x_ref[idx_ion] = 0.0
    solv_sum = np.sum(x_ref[idx_solv])
    if solv_sum > 0:
        x_ref[idx_solv] = x_ref[idx_solv] / solv_sum
    else:
        x_ref[idx_solv] = 1.0 / len(idx_solv)

    eps = 1e-12
    p = pcsaft_p_cpp(t, rho, x_ref, cppargs)
    phase = 1 if rho < 900 else 0
    result = {}
    for i in idx_ion:
        x_inf = x_ref.copy()
        x_inf[i] = eps
        x_inf /= np.sum(x_inf)
        rho_inf = pcsaft_den_cpp(t, p, x_inf, phase, cppargs)
        lnfug_inf = float(pcsaft_lnfug_cpp(t, rho_inf, x_inf, cppargs)[i])
        if not np.isfinite(lnfug_inf):
            raise SolutionError('Non-finite ln(fugacity coefficient) at infinite dilution.')
        result[species[i]] = 8.31446261815324 * t * lnfug_inf
    return result


def pcsaft_ares(t, rho, x, params):
    """
    Calculate the residual Helmholtz energy.

    Parameters
    ----------
    t : float
        Temperature (K)
    rho : float
        Molar density (mol m\ :sup:`-3`)
    x : ndarray, shape (n,)
        Mole fractions of each component. It has a length of n, where n is
        the number of components in the system.
    params : dict
        A dictionary containing PC-SAFT parameters that can be passed for
        use in PC-SAFT:

        m : ndarray, shape (n,)
            Segment number for each component.
        s : ndarray, shape (n,)
            Segment diameter for each component. For ions this is the diameter of
            the hydrated ion. Units of Angstrom.
        e : ndarray, shape (n,)
            Dispersion energy of each component. For ions this is the dispersion
            energy of the hydrated ion. Units of K.
        k_ij : ndarray, shape (n,n)
            Binary interaction parameters between components in the mixture.
            (dimensions: ncomp x ncomp)
        e_assoc : ndarray, shape (n,)
            Association energy of the associating components. For non associating
            compounds this is set to 0. Units of K.
        vol_a : ndarray, shape (n,)
            Effective association volume of the associating components. For non
            associating compounds this is set to 0.
        dipm : ndarray, shape (n,)
            Dipole moment of the polar components. For components where the dipole
            term is not used this is set to 0. Units of Debye.
        dip_num : ndarray, shape (n,)
            The effective number of dipole functional groups on each component
            molecule. Generally this is set to 1, but some implementations use this
            as an adjustable parameter that is fit to data.
        z : ndarray, shape (n,)
            Charge number of the ions
        dielc : float
            Dielectric constant of the medium to be used for electrolyte
            calculations.
        assoc_scheme : list, shape (n,)
            The types of association sites for each component. Use `None` for molecules
            without association sites. If a molecule has multiple association sites,
            use a nested list for that component to specify the association scheme for
            each site. The accepted association schemes are those given by Huang and
            Radosz (1990): 1, 2A, 2B, 3A, 3B, 4A, 4B, 4C. If `e_assoc` and `vol_a` are
            given but `assoc_scheme` is not, the 2B association scheme is assumed (which
            would, for example, correspond to one hydroxyl functional group).

    Returns
    -------
    ares : float
        Residual Helmholtz energy (J mol\ :sup:`-1`)
    """
    x, params = ensure_numpy_input(x, params)
    check_input(x, {'density':rho, 'temperature':t})
    params = check_association(params)
    cppargs = create_struct(params)
    return pcsaft_ares_cpp(t, rho, x, cppargs)


def pcsaft_dadt(t, rho, x, params):
    """
    Calculate the temperature derivative of the residual Helmholtz energy.

    Parameters
    ----------
    t : float
        Temperature (K)
    rho : float
        Molar density (mol m\ :sup:`-3`)
    x : ndarray, shape (n,)
        Mole fractions of each component. It has a length of n, where n is
        the number of components in the system.
    params : dict
        A dictionary containing PC-SAFT parameters that can be passed for
        use in PC-SAFT:

        m : ndarray, shape (n,)
            Segment number for each component.
        s : ndarray, shape (n,)
            Segment diameter for each component. For ions this is the diameter of
            the hydrated ion. Units of Angstrom.
        e : ndarray, shape (n,)
            Dispersion energy of each component. For ions this is the dispersion
            energy of the hydrated ion. Units of K.
        k_ij : ndarray, shape (n,n)
            Binary interaction parameters between components in the mixture.
            (dimensions: ncomp x ncomp)
        e_assoc : ndarray, shape (n,)
            Association energy of the associating components. For non associating
            compounds this is set to 0. Units of K.
        vol_a : ndarray, shape (n,)
            Effective association volume of the associating components. For non
            associating compounds this is set to 0.
        dipm : ndarray, shape (n,)
            Dipole moment of the polar components. For components where the dipole
            term is not used this is set to 0. Units of Debye.
        dip_num : ndarray, shape (n,)
            The effective number of dipole functional groups on each component
            molecule. Generally this is set to 1, but some implementations use this
            as an adjustable parameter that is fit to data.
        z : ndarray, shape (n,)
            Charge number of the ions
        dielc : float
            Dielectric constant of the medium to be used for electrolyte
            calculations.
        assoc_scheme : list, shape (n,)
            The types of association sites for each component. Use `None` for molecules
            without association sites. If a molecule has multiple association sites,
            use a nested list for that component to specify the association scheme for
            each site. The accepted association schemes are those given by Huang and
            Radosz (1990): 1, 2A, 2B, 3A, 3B, 4A, 4B, 4C. If `e_assoc` and `vol_a` are
            given but `assoc_scheme` is not, the 2B association scheme is assumed (which
            would, for example, correspond to one hydroxyl functional group).

    Returns
    -------
    dadt : float
        Temperature derivative of the residual Helmholtz energy (J mol\ :sup:`-1`)
    """
    x, params = ensure_numpy_input(x, params)
    check_input(x, {'density':rho, 'temperature':t})
    params = check_association(params)
    cppargs = create_struct(params)
    return pcsaft_dadt_cpp(t, rho, x, cppargs)


def aly_lee(t, c):
    """
    Calculate the ideal gas isobaric heat capacity using the Aly-Lee equation.

    Parameters
    ----------
    t : float
        Temperature (K)
    c : ndarray, shape (5,)
        Constants for the Aly-Lee equation

    Returns
    -------
    cp_ideal : float
        Ideal gas isobaric heat capacity (J mol\ :sup:`-1` K\ :sup:`-1`)

    References
    ----------
    - F. A. Aly and L. L. Lee, “Self-consistent equations for calculating the ideal gas heat capacity, enthalpy, and entropy,” Fluid Phase Equilibria, vol. 6, no. 3–4, pp. 169–179, 1981.
    """
    cp_ideal = (c[0] + c[1]*(c[2]/t/np.sinh(c[2]/t))**2 + c[3]*(c[4]/t/np.cosh(c[4]/t))**2)/1000.
    return cp_ideal

def dielc_water(t):
    """
    Return the dielectric constant of water at the given temperature.

    This equation was fit to values given in the reference. For temperatures from
    263.15 to 368.15 K values at 1 bar were used. For temperatures from 368.15 to
    443.15 K values at 10 bar were used. Below 263.15 K and above 443.15 K an
    error is raised.

    Parameters
    ----------
    t : float
        Temperature (K)

    Returns
    -------
    dielc : float
        Dielectric constant of water

    References
    ----------
    - D. G. Archer and P. Wang, “The Dielectric Constant of Water and Debye‐Hückel Limiting Law Slopes,” J. Phys. Chem. Ref. Data, vol. 19, no. 2, pp. 371–411, Mar. 1990.
    """
    if t < 263.15:
        raise ValueError('For dielc_water t must be greater than 263.15 K.')
    elif t > 443.15:
        raise ValueError('For dielc_water t must be less than 443.15 K.')

    if t <= 368.15:
        dielc = 7.6555618295E-04*t**2 - 8.1783881423E-01*t + 2.5419616803E+02
    else:
        dielc = 0.0005003272124*t**2 - 0.6285556029*t + 220.4467027
    return dielc


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
    if 'z' in params:
        cppargs.z = np_to_vector_double(params['z'])
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
    cppargs.dielc_rule = int(params['dielc_rule']) if 'dielc_rule' in params else 1
    cppargs.dielc_diff_mode = int(params['dielc_diff_mode']) if 'dielc_diff_mode' in params else 0
    if cppargs.dielc_diff_mode not in (0, 1):
        raise ValueError("Unknown dielc_diff_mode. Supported values are 0 (analytic) and 1 (finite-diff).")
    if cppargs.z.size() > 0 and cppargs.dielc.size() == 0:
        raise ValueError('Electrolyte parameters require params["dielc"] as a per-species array.')
    if 'd_born' in params:
        cppargs.d_born = np_to_vector_double(np.asarray(params['d_born'], dtype=float))
    if 'f_solv' in params:
        cppargs.f_solv = np_to_vector_double(np.asarray(params['f_solv'], dtype=float))
    cppargs.born_model = int(params['born_model']) if 'born_model' in params else 1
    cppargs.born_diff_mode = int(params['born_diff_mode']) if 'born_diff_mode' in params else 0
    if cppargs.born_diff_mode not in (0, 1, 2):
        raise ValueError("Unknown born_diff_mode. Supported values are 0 (analytic), 1 (finite-diff), and 2 (Eq.133-style).")
    cppargs.DH_model = int(params['DH_model']) if 'DH_model' in params else 1
    if cppargs.DH_model == 2:
        raise ValueError("DH_model=2 (Bjerrum treatment) is reserved and not implemented.")
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


def pcsaft_dielc_eval(x, params):
    """
    Evaluate mixed dielectric constant and composition derivatives using the C++ dielectric engine.
    """
    x, params = ensure_numpy_input(x, params)
    check_input(x, {})
    params = check_association(params)
    cppargs = create_struct(params)
    eps = pcsaft_dielc_eps_cpp(x, cppargs)
    deps = np.asarray(pcsaft_dielc_diff_cpp(x, cppargs))
    return eps, deps
