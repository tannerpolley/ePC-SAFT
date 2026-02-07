import numpy as np
import types
import math

pcsaft_prop = {
    'CO2': {
        'MW': 44.01e-3,  # kg/mol
        'm': 2.079, 's': 2.7852, 'e': 169.21,
        'e_assoc': 0., 'vol_a': 0., 'assoc_scheme': None,
        'dipm': 0., 'dip_num': 1,
        'z': 0., 'dielc': 1.4122 # Schick 2022
    },
    'MEA-2B': {
        'MW': 61.08e-3,  # kg/mol
        'm': 3.0353, 's': 3.0435, 'e': 277.174,
        'e_assoc': 2586.3, 'vol_a': 0.037470, 'assoc_scheme': '2B',
        'dipm': 0., 'dip_num': 1,
        'z': 0., 'dielc': 32.
    },
    'MEA-4C': {
        'MW': 61.08e-3,  # kg/mol
        'm': 4.5208, 's': 2.6574, 'e': 237.6864,
        'e_assoc': 989.8984, 'vol_a': 0.187533, 'assoc_scheme': '4C',
        'dipm': 0., 'dip_num': 1,
        'z': 0., 'dielc': 0.
    },
    'H2O-2B-CC': {
        'MW': 18.01528e-3,  # kg/mol
        'm': 1.9599, 's': 2.362, 'e': 279.42,
        'e_assoc': 2059.28, 'vol_a': 0.1750, 'assoc_scheme': '2B',
        'dipm': 0., 'dip_num': 1,
        'z': 0.,
        'dielc': 78.09,
    },
    'H2O-4C-CC': {
        'MW': 18.01528e-3,  # kg/mol
        'm': 2.1945, 's': 2.229, 'e': 141.66,
        'e_assoc': 1804.17, 'vol_a': 0.2039, 'assoc_scheme': '4C',
        'dipm': 0., 'dip_num': 1,
        'z': 0.,
        'dielc': 78.09,
    },
    'MEAH+': {
        'MW': 62.09e-3,  # kg/mol
        'm': 1., 's': 3.0435, 'e': 277.174,
        'e_assoc': 0., 'vol_a': 0., 'assoc_scheme': None,
        'dipm': 0., 'dip_num': 1,
        'z': 1., 'dielc': 8.
    },
    'MEACOO-': {
        'MW': 75.07e-3,  # kg/mol
        'm': 1., 's': 3.0435, 'e': 277.174,
        'e_assoc': 0., 'vol_a': 0., 'assoc_scheme': None,
        'dipm': 0., 'dip_num': 1,
        'z': -1., 'dielc': 8.
    },
    'HCO3-': {
        'MW': 61.0168e-3,  # kg/mol
        'm': 1., 's': 3., 'e': 300.,
        'e_assoc': 0., 'vol_a': 0., 'assoc_scheme': None,
        'dipm': 0., 'dip_num': 1,
        'z': -1., 'dielc': 8.
    },
    'CO32-': {
        'MW': 60.01e-3,  # kg/mol
        'm': 1., 's': 3., 'e': 300.,
        'e_assoc': 0., 'vol_a': 0., 'assoc_scheme': None,
        'dipm': 0., 'dip_num': 1,
        'z': -1., 'dielc': 8.
    },
    'H3O+': {
        'MW': 19.02e-3,  # kg/mol
        'm': 1., 's': 3., 'e': 300.,
        'e_assoc': 0., 'vol_a': 0., 'assoc_scheme': None,
        'dipm': 0., 'dip_num': 1,
        'z': 1., 'dielc': 8.
    },
    'OH-': {
        'MW': 17.01e-3,  # kg/mol
        'm': 1., 's': 3., 'e': 300.,
        'e_assoc': 0., 'vol_a': 0., 'assoc_scheme': None,
        'dipm': 0., 'dip_num': 1,
        'z': -1., 'dielc': 8.
    },

    'Hexane': {
        'MW': 86.17848e-3,  # kg/mol
        'm': 3.0576, 's': 3.7983, 'e': 236.77,
        'e_assoc': 0.0, 'vol_a': 0.0, 'assoc_scheme': None,
        'dipm': 0., 'dip_num': 1,
        'z': 0., 'dielc': 0.
    },
    'Methane': {
        'MW': 16.04e-3,  # kg/mol
        'm': 1.0, 's': 3.7039, 'e': 150.03,
        'e_assoc': 0., 'vol_a': 0., 'assoc_scheme': None,
        'dipm': 0., 'dip_num': 1,
        'z': 0., 'dielc': 0.},
    'Ethane': {
        'MW': 30.07e-3,  # kg/mol
        'm': 1.6069, 's': 3.5206, 'e': 191.42,
        'e_assoc': 0., 'vol_a': 0., 'assoc_scheme': None,
        'dipm': 0., 'dip_num': 1,
        'z': 0., 'dielc': 0.},
    'Propane': {
        'MW': 44.10e-3,  # kg/mol
        'm': 2.0020, 's': 3.6184, 'e': 208.11,
        'e_assoc': 0., 'vol_a': 0., 'assoc_scheme': None,
        'dipm': 0., 'dip_num': 1,
        'z': 0., 'dielc': 0.},
    'Methanol': {
        'MW': 32.04e-3,  # kg/mol
        'm': 1.5255, 's': 3.2300, 'e': 188.90,
        'e_assoc': 2899.5, 'vol_a': 0.03518, 'assoc_scheme': '2B',
        'dipm': 0., 'dip_num': 1,
        'z': 0., 'dielc': 33.05,
        'f_solv': 1.4,
    },
    'Ethanol': {
        'MW': 46.068e-3,  # kg/mol
        'm': 2.3827, 's': 3.1771, 'e': 198.24,
        'e_assoc': 2653.4, 'vol_a': 0.03238, 'assoc_scheme': '2B',
        'dipm': 0., 'dip_num': 1,
        'z': 0., 'dielc': 24.88,
        'f_solv': 1.6,
    },
    'Butanol': {
        'MW': 74.12e-3,  # kg/mol
        'm': 2.7510, 's': 3.6139, 'e': 259.59,
        'e_assoc': 2544.56, 'vol_a': 0.00669, 'assoc_scheme': '2B',
        'dipm': 0., 'dip_num': 1,
        'z': 0., 'dielc': 20.47},

    'H2O-2B-Li': {
        'MW': 18.01528e-3,  # kg/mol
        'm': 1.2047, 's': lambda T: 2.7927 + (10.11 * np.exp(-.01775 * T) - 1.417 * np.exp(-.01146 * T)), 'e': 353.9449,
        'e_assoc': 2425.7, 'vol_a': .04509, 'assoc_scheme': '2B',
        'dipm': 0., 'dip_num': 1,
        'z': 0.,
        # 'dielc': lambda T: -105.2*np.log(T) + 677.480,
        'dielc': 78.09,
        'f_solv': 1.5,
    },

    'H2O-2B-NaCl': {
        'MW': 18.01528e-3,  # kg/mol
        'm': 1.2047, 's': lambda T: 2.7927 + (10.11 * np.exp(-.01775 * T) - 1.417 * np.exp(-.01146 * T)), 'e': 353.9449,
        'e_assoc': 2425.67, 'vol_a': .0451, 'assoc_scheme': '2B',
        'dipm': 0., 'dip_num': 1,
        'z': 0.,
        'dielc': 78.09,
        'f_solv': 1.5,
    },

    'TOP': {
        'MW': 434.63e-3,  # kg/mol
        'm': 4.2032, 's': 5.4506, 'e': 280.4777,
        'e_assoc': 6393.5, 'vol_a': .0001, 'assoc_scheme': '2B',
        'dipm': 0., 'dip_num': 1,
        'z': 0., 'dielc': 11
    },
    'IL': {
        'MW': 407.31e-3,  # kg/mol
        'm': 4.073, 's': 4.6432, 'e': 434.6130,
        'e_assoc': 5000, 'vol_a': .1, 'assoc_scheme': '2B',
        'dipm': 0., 'dip_num': 1,
        'z': 0., 'dielc': 11
    },

    'Li+': {
        'MW': 6.94e-3,  # kg/mol
        'm': 1., 's': 2.8449, 'e': 360.00,
        'e_assoc': 0.0, 'vol_a': 100, 'assoc_scheme': '2B',
        'dipm': 0., 'dip_num': 1,
        'z': 1, 'dielc': 8,
        'd_born': 2.784,
        'f_solv': 1.0
    },
    'Na+': {
        'MW': 22.98e-3,  # kg/mol
        'm': 1., 's': 2.8232, 'e': 230.00,
        'e_assoc': 0.0, 'vol_a': 0.0, 'assoc_scheme': None,
        'dipm': 0., 'dip_num': 1,
        'z': 1, 'dielc': 8,
        'd_born': 3.445,
        'f_solv': 1.0
    },
    # 'Na+': {
    #     'MW': 22.98e-3,  # kg/mol
    #     'm': 1., 's': 2.4122, 'e': 230.00,
    #     'e_assoc': 0.0, 'vol_a': 0.0, 'assoc_scheme': None,
    #     'dipm': 0., 'dip_num': 1,
    #     'z': 1, 'dielc': 8,
    #     'd_born': 3.445,
    #     'f_solv': 1.0
    # },
    'K+': {
        'MW': 39.0983e-3,  # kg/mol
        'm': 1., 's': 3.3417, 'e': 200.00,
        'e_assoc': 0.0, 'vol_a': 0.0, 'assoc_scheme': None,
        'dipm': 0., 'dip_num': 1,
        'z': 1, 'dielc': 8,
        'd_born': 4.150,
        'f_solv': 1.0
    },
    'Mg2+': {
        'MW': 24.31e-3,  # kg/mol
        'm': 1., 's': 3.1327, 'e': 1500,
        'e_assoc': 0., 'vol_a': 0, 'assoc_scheme': None,
        'dipm': 0., 'dip_num': 1,
        'z': 2, 'dielc': 8
    },
    'F-': {
        'MW': 18.998e-3,  # kg/mol
        'm': 1., 's': 1.7712, 'e': 275.00,
        'e_assoc': 0., 'vol_a': 0, 'assoc_scheme': None,
        'dipm': 0., 'dip_num': 1,
        'z': -1, 'dielc': 8
    },
    'Cl-': {
        'MW': 35.45e-3,  # kg/mol
        'm': 1., 's': 2.7560, 'e': 170.00,
        'e_assoc': 0., 'vol_a': 0, 'assoc_scheme': None,
        'dipm': 0., 'dip_num': 1,
        'z': -1, 'dielc': 8,
        'd_born': 4.100,
        'f_solv': 1.0
    },
    'Br-': {
        'MW': 79.904e-3,  # kg/mol
        'm': 1., 's': 3.0707, 'e': 190.00,
        'e_assoc': 0., 'vol_a': 0, 'assoc_scheme': None,
        'dipm': 0., 'dip_num': 1,
        'z': -1, 'dielc': 8,
        'd_born': 4.480,
        'f_solv': 1.0
    },
    'I-': {
        'MW': 126.90447e-3,  # kg/mol
        'm': 1., 's': 3.6672, 'e': 200.00,
        'e_assoc': 0., 'vol_a': 0, 'assoc_scheme': None,
        'dipm': 0., 'dip_num': 1,
        'z': -1, 'dielc': 8,
        'd_born': 4.985,
        'f_solv': 1.0
    },

}

# Defaults for components without explicit Born/solvation-shell parameters
for _sp, _props in pcsaft_prop.items():
    _props.setdefault('d_born', 0.0)
    _props.setdefault('f_solv', 1.0)


# Create the binary interaction parameter dictionary for dispersion forces

k_ij_dict = {

    # CO2-MEA-H2O System
    ("CO2", "H2O"): lambda T: -2.2e-2 + 4.2e-4 * (T - 298) - 1.7e-6 * (T - 298),
    ("CO2", "MEA"): 0.0,
    # ("MEA-2B", "H2O-2B-CC"): -0.0420, # Baygi 2015
    ("MEA-2B", "H2O-2B-CC"): 0.250,  # Baygi 2015
    ("MEAH+", "MEACOO-"): 0.0,

    # Example System for hydrocarbons from LearnChemE
    ("Methane", "Ethane"): 3e-4,
    ("Methane", "Propane"): 1.15e-2,
    ("Ethane", "Propane"): 5.10e-3,

    # Lithium Extraction with Ionic Liquids from Yu 2024 10.1016/j.ces.2023.119682
    ("H2O-2B-Li", "IL"): .007,
    ("Li+", "TOP"): .3,
    ("H2O-2B-Li", "TOP"): 1,
    ("TOP", "IL"): 1,
    ("Li+", "IL"): 1,

    # Testing from pubs.acs.org/IECR Article
    # Predicting Thermodynamic Properties of Ions in Single Solvents
    # and in Mixed Solvents Using a Modified Born Term within the ePC-SAFT Framework
    # 10.1021/acs.iecr.5c00475
    # ("Li+", "H2O"): -.2500,
    # ("Li+", "H2O-2B-Li"): -.4,
    # ("Na+", "H2O-2B-Li"): .0045,
    # ("K+", "H2O-2B-Li"): .1997,
    # ("F-", "H2O-2B-Li"): 0.000,
    # ("Cl-", "H2O-2B-Li"): -.250,
    # ("Br-", "H2O-2B-Li"): -.250,
    # ("I-", "H2O-2B-Li"): -.250,

    ("H2O-2B-Li", "Methanol"): -.0878,
    ("H2O-2B-Li", "Ethanol"): -.0878,
    ("H2O-2B-Li", "Butanol"): lambda T: 2.94e-4 * T - .102,

    # Figiel 2025 - Predicting Thermodynamic Properties of Ions in Single Solvents
    # and in Mixed Solvents Using a Modified Born Term within the ePC-SAFT Framework
    # 10.1021/acs.iecr.5c00475
    # Water
    # Cation - Solvent
    ("H+", "H2O-2B-Li"): 0.0,
    ("Li+", "H2O-2B-Li"): -.4,
    ("Na+", "H2O-2B-Li"): -.3,
    ("K+", "H2O-2B-Li"): -.1,
    # Anion - Solvent
    ("Cl-", "H2O-2B-Li"): -.3,
    ("Br-", "H2O-2B-Li"): -.3,
    ("I-", "H2O-2B-Li"): -.05,

    # Methanol
    # Cation - Solvent
    ("H+", "Methanol"): -.3,
    ("Li+", "Methanol"): -.9,
    ("Na+", "Methanol"): -.25,
    ("K+", "Methanol"): .32,
    # Anion - Solvent
    ("Cl-", "Methanol"): .5,
    ("Br-", "Methanol"): .15,
    ("I-", "Methanol"): .37,

    # Ethanol
    # Cation - Solvent
    ("H+", "Ethanol"): -.6,
    ("Li+", "Ethanol"): -.8,
    ("Na+", "Ethanol"): .05,
    ("K+", "Ethanol"): .53,
    # Anion - Solvent
    ("Cl-", "Ethanol"): .8,
    ("Br-", "Ethanol"): 0.,
    ("I-", "Ethanol"): .18,

    # Cation - Anion
    ("H+", "Cl-"): -.9,
    ("H+", "Br-"): -.7,
    ("Li+", "Cl-"): .8,
    ("Li+", "Br-"): .5,
    ("Na+", "Cl-"): .8,
    ("Na+", "Br-"): .65,
    ("Na+", "I-"): .45,
    ("K+", "Cl-"): 0.,
    ("K+", "Br-"): -.35,
    ("K+", "I-"): 0.,

    # Water
    # Cation - Solvent
    ("Na+", "H2O-2B-NaCl"): .0045,
    # Anion - Solvent
    ("Cl-", "H2O-2B-NaCl"): -.25,
    ("Na+", "Cl-"): .317,
}

k_ij_dict[("H2O-2B-CC", "MEAH+")] = k_ij_dict[("MEA-2B", "H2O-2B-CC")]
k_ij_dict[("H2O-2B-CC", "MEACOO-")] = k_ij_dict[("MEA-2B", "H2O-2B-CC")]

unique_strings = set()

for key in k_ij_dict:
    unique_strings.update(key)  # key is a tuple, so this adds both elements

# Convert to list if needed
unique_list = list(unique_strings)
for k in unique_list:
    k_ij_dict[(k, k)] = 0.0

for k in k_ij_dict.copy().keys():
    k1, k2 = k
    k_ij_dict[(k2, k1)] = k_ij_dict[(k1, k2)]

# Create the binary interaction parameter dictionary for association forces

k_hb_dict = {
    ("Li+", "TOP"): .3,
    ("Li+", "IL"): 1,
    ("Li+", "H2O-2B-Li"): 1,
    ("Butanol", "H2O-2B-Li"): .026,
}

unique_strings_hb = set()

for key in k_hb_dict:
    unique_strings_hb.update(key)  # key is a tuple, so this adds both elements

# Convert to list if needed
unique_list_hb = list(unique_strings_hb)
for k in unique_list_hb:
    k_hb_dict[(k, k)] = 0.0

for k in k_hb_dict.copy().keys():
    k1, k2 = k
    k_hb_dict[(k2, k1)] = k_hb_dict[(k1, k2)]

l_ij_dict = {
    ("H2O-2B-Li", "Butanol"): -.0044,
}

unique_strings_l_ij = set()

for key in l_ij_dict:
    unique_strings_l_ij.update(key)  # key is a tuple, so this adds both elements

# Convert to list if needed
unique_list_l_ij = list(unique_strings_l_ij)
for k in unique_list_l_ij:
    l_ij_dict[(k, k)] = 0.0

for k in l_ij_dict.copy().keys():
    k1, k2 = k
    l_ij_dict[(k2, k1)] = l_ij_dict[(k1, k2)]


BASE_KEYS = [
    'MW', 'm', 's', 'e',
    'e_assoc', 'vol_a', 'assoc_scheme',
    'dipm', 'dip_num', 'z', 'dielc'
]
OPTIONAL_KEYS = ['d_born', 'f_solv']


def _resolve_species_params(sp, user_params):
    if user_params is not None and sp in user_params:
        return user_params[sp], 'user_params'
    if sp in pcsaft_prop:
        return pcsaft_prop[sp], 'default'
    raise KeyError("Species '{}' not found in user_params or default pcsaft_prop.".format(sp))


def _get_value(entry, prop, T):
    value = entry[prop]
    if isinstance(value, types.FunctionType):
        return value(T)
    return value


def get_prop_dict(species, x, T, user_params=None, user_options=None):
    """
    species: list of species names that match dictionary keys in pcsaft_prop
    T: Temperature (K) (often not used, used in calculations of temperature-dependent binary interaction parameters)
    user_params: optional dict in the form {component: {m, s, e, ...}}
    """

    # DH_model 0 and 1 are aliases in the current C++ implementation; 2 is reserved.
    default_options = {
        "dielc_rule": 1,
        "dielc_diff_mode": 0,
        "born_model": 1,
        "DH_model": 1,
        "debug": False,
    }
    if user_options is None:
        user_options = {}
    unknown_options = set(user_options) - set(default_options)
    if unknown_options:
        raise KeyError("Unknown user_options key(s): {}. Allowed keys: {}.".format(
            sorted(unknown_options), sorted(default_options.keys())))
    options = default_options.copy()
    options.update(user_options)

    prop_dic = {}
    entries = []
    for sp in species:
        entry, _ = _resolve_species_params(sp, user_params)
        entries.append(entry)

    for prop in BASE_KEYS:
        prop_list = []
        for sp, entry in zip(species, entries):
            if prop in entry:
                prop_list.append(_get_value(entry, prop, T))
                continue
            if prop == 'dielc':
                z_i = float(_get_value(entry, 'z', T)) if 'z' in entry else 0.0
                if abs(z_i) > 1e-12:
                    prop_list.append(8.0)
                    continue
            raise KeyError("Missing '{}' for species '{}' in {}.".format(prop, sp, 'user_params' if (user_params is not None and sp in user_params) else 'default pcsaft_prop'))
        if prop == 'assoc_scheme':
            prop_dic[prop] = prop_list
        else:
            prop_dic[prop] = np.array(prop_list)

    for prop in OPTIONAL_KEYS:
        if any(prop in entry for entry in entries):
            prop_list = []
            for sp, entry in zip(species, entries):
                prop_list.append(_get_value(entry, prop, T) if prop in entry else 0.0)
            prop_dic[prop] = np.array(prop_list)

    n = len(species)

    # Create the binary interaction parameter dictionary and matrix for dispersion forces
    k_ij = np.zeros((n, n))
    for i, sp1 in enumerate(species):
        for j, sp2 in enumerate(species):
            try:
                if isinstance(k_ij_dict[(sp1, sp2)], types.FunctionType):
                    k_ij[i, j] = k_ij_dict[(sp1, sp2)](T)
                else:
                    k_ij[i, j] = k_ij_dict[(sp1, sp2)]
            except KeyError:
                k_ij[i, j] = 0.0
    prop_dic['k_ij'] = k_ij

    # Create the binary interaction parameter dictionary and matrix for association forces
    assoc_species = []
    for sp, entry in zip(species, entries):
        if entry.get('assoc_scheme') is not None:
            assoc_species.append(sp)
    k_hb = np.zeros((n, n))
    for i, sp1 in enumerate(assoc_species):
        for j, sp2 in enumerate(assoc_species):
            try:
                k_hb[i, j] = k_hb_dict[(sp1, sp2)]
            except KeyError:
                k_hb[i, j] = 0.0
    prop_dic['k_hb'] = k_hb

    l_ij = np.zeros((n, n))
    for i, sp1 in enumerate(species):
        for j, sp2 in enumerate(species):
            try:
                if isinstance(l_ij_dict[(sp1, sp2)], types.FunctionType):
                    l_ij[i, j] = l_ij_dict[(sp1, sp2)](T)
                else:
                    l_ij[i, j] = l_ij_dict[(sp1, sp2)]
            except KeyError:
                l_ij[i, j] = 0.0
    prop_dic['l_ij'] = l_ij

    if np.all(prop_dic['z'] == 0):
        prop_dic['z'] = np.array([])

    prop_dic['born_model'] = options['born_model']
    prop_dic['DH_model'] = options['DH_model']
    prop_dic['dielc_rule'] = options['dielc_rule']
    prop_dic['dielc_diff_mode'] = options['dielc_diff_mode']
    prop_dic['debug'] = bool(options['debug'])

    return prop_dic


def validate_species_params(species, user_params=None):
    """
    Validate that species exist and required keys are present.

    Returns a dict with:
      - missing_species: list of species not found in user_params or default
      - missing_keys: dict of species -> list of missing required keys
    """
    missing_species = []
    missing_keys = {}

    for sp in species:
        entry = None
        if user_params is not None and sp in user_params:
            entry = user_params[sp]
        elif sp in pcsaft_prop:
            entry = pcsaft_prop[sp]
        else:
            missing_species.append(sp)
            continue

        missing = [k for k in BASE_KEYS if k not in entry]
        if missing:
            missing_keys[sp] = missing

    return {
        "missing_species": missing_species,
        "missing_keys": missing_keys,
    }


def _prepare_dielc_inputs(x, dielc, mw=None, z=None):
    x = np.asarray(x, dtype=float)
    dielc = np.asarray(dielc, dtype=float)
    if x.size != dielc.size:
        raise ValueError("x and dielc must have the same length.")
    if mw is not None:
        mw = np.asarray(mw, dtype=float)
        if mw.size != x.size:
            raise ValueError("mw must have the same length as x.")
    if z is not None:
        z = np.asarray(z, dtype=float)
        if z.size != x.size:
            raise ValueError("z must have the same length as x.")
    return x, dielc, mw, z


def dielc_rule(x, dielc, rule=1, mw=None, z=None):
    # Utility/reference implementation. Runtime electrolyte calculations now use the C++ dielectric engine.
    x, dielc, mw, z = _prepare_dielc_inputs(x, dielc, mw=mw, z=z)

    alpha = 7.01

    if rule == 0:
        return float(np.max(dielc))

    if rule == 1:
        # Mole-fraction mixing.
        return float(np.dot(x, dielc))

    if rule == 2:
        # Mass-fraction mixing.
        if mw is None:
            raise ValueError("dielc_rule=2 requires mw.")
        mw_bar = float(np.dot(x, mw))
        if mw_bar <= 0:
            raise ValueError("Average molecular weight must be positive for dielc_rule=2.")
        return float(np.dot(x * mw, dielc) / mw_bar)

    if rule == 3:
        # Combo rule: solvents mass-fraction mixed, ions mole-fraction mixed.
        if mw is None or z is None:
            raise ValueError("dielc_rule=3 requires mw and z.")
        idx_sol = np.where(np.abs(z) <= 1e-12)[0]
        idx_ion = np.where(np.abs(z) > 1e-12)[0]
        if idx_sol.size == 0:
            raise ValueError("dielc_rule=3 requires at least one solvent species (z=0).")
        mw_sol = float(np.dot(x[idx_sol], mw[idx_sol]))
        if mw_sol <= 0:
            raise ValueError("Solvent molecular-weight denominator must be positive for dielc_rule=3.")
        eps_sol_w = float(np.dot(x[idx_sol] * mw[idx_sol], dielc[idx_sol]) / mw_sol)
        x_sol = float(np.sum(x[idx_sol]))
        eps_ion = float(np.dot(x[idx_ion], dielc[idx_ion])) if idx_ion.size > 0 else 0.0
        return float(x_sol * eps_sol_w + eps_ion)

    if rule == 4:
        # New mixing rule (Figiel 2025).
        if mw is None or z is None:
            raise ValueError("dielc_rule=4 requires mw and z.")
        idx_sol = np.where(np.abs(z) <= 1e-12)[0]
        idx_ion = np.where(np.abs(z) > 1e-12)[0]
        if idx_sol.size == 0:
            raise ValueError("dielc_rule=4 requires at least one solvent species (z=0).")
        mw_sol = float(np.dot(x[idx_sol], mw[idx_sol]))
        if mw_sol <= 0:
            raise ValueError("Solvent molecular-weight denominator must be positive for dielc_rule=4.")
        eps_sf = float(np.dot(x[idx_sol] * mw[idx_sol], dielc[idx_sol]) / mw_sol)
        x_ion = float(np.sum(x[idx_ion])) if idx_ion.size > 0 else 0.0
        return float(eps_sf / (1.0 + alpha * x_ion))

    if rule == 5:
        # Same eps_mix as rule 4; derivative behavior is overridden in dielc_diff_rule.
        if mw is None or z is None:
            raise ValueError("dielc_rule=5 requires mw and z.")
        idx_sol = np.where(np.abs(z) <= 1e-12)[0]
        idx_ion = np.where(np.abs(z) > 1e-12)[0]
        if idx_sol.size == 0:
            raise ValueError("dielc_rule=5 requires at least one solvent species (z=0).")
        mw_sol = float(np.dot(x[idx_sol], mw[idx_sol]))
        if mw_sol <= 0:
            raise ValueError("Solvent molecular-weight denominator must be positive for dielc_rule=5.")
        eps_sf = float(np.dot(x[idx_sol] * mw[idx_sol], dielc[idx_sol]) / mw_sol)
        x_ion = float(np.sum(x[idx_ion])) if idx_ion.size > 0 else 0.0
        return float(eps_sf / (1.0 + alpha * x_ion))

    if rule == 6:
        # Same functional form as rule 4, but eps_sf is treated as concentration-independent.
        if z is None:
            raise ValueError("dielc_rule=6 requires z.")
        idx_sol = np.where(np.abs(z) <= 1e-12)[0]
        idx_ion = np.where(np.abs(z) > 1e-12)[0]
        if idx_sol.size == 0:
            raise ValueError("dielc_rule=6 requires at least one solvent species (z=0).")
        eps_sf_const = float(np.mean(dielc[idx_sol]))
        x_ion = float(np.sum(x[idx_ion])) if idx_ion.size > 0 else 0.0
        return float(eps_sf_const / (1.0 + alpha * x_ion))

    raise ValueError("Unknown dielc_rule: {}. Supported rules: 0, 1, 2, 3, 4, 5, 6.".format(rule))


def dielc_diff_rule(x, dielc, rule=1, mw=None, z=None):
    # Utility/reference implementation. Runtime electrolyte calculations now use the C++ dielectric engine.
    x, dielc, mw, z = _prepare_dielc_inputs(x, dielc, mw=mw, z=z)

    alpha = 7.01

    if rule == 0:
        return np.zeros(len(x), dtype=float)

    if rule == 1:
        # d/dx_i sum_j x_j eps_j = eps_i
        return np.asarray(dielc, dtype=float)

    if rule == 2:
        # (MW_i / MW_bar) * (eps_i - eps_mix)
        if mw is None:
            raise ValueError("dielc_rule=2 requires mw.")
        mw_bar = float(np.dot(x, mw))
        if mw_bar <= 0:
            raise ValueError("Average molecular weight must be positive for dielc_rule=2.")
        eps_mix = float(np.dot(x * mw, dielc) / mw_bar)
        return np.asarray((mw / mw_bar) * (dielc - eps_mix), dtype=float)

    if rule == 3:
        # Solvent components use combo-rule derivative; ions contribute with eps_i.
        if mw is None or z is None:
            raise ValueError("dielc_rule=3 requires mw and z.")
        idx_sol = np.where(np.abs(z) <= 1e-12)[0]
        idx_ion = np.where(np.abs(z) > 1e-12)[0]
        if idx_sol.size == 0:
            raise ValueError("dielc_rule=3 requires at least one solvent species (z=0).")
        mw_sol = float(np.dot(x[idx_sol], mw[idx_sol]))
        if mw_sol <= 0:
            raise ValueError("Solvent molecular-weight denominator must be positive for dielc_rule=3.")
        eps_sol_w = float(np.dot(x[idx_sol] * mw[idx_sol], dielc[idx_sol]) / mw_sol)
        x_sol = float(np.sum(x[idx_sol]))
        deps = np.zeros(len(x), dtype=float)
        deps[idx_sol] = eps_sol_w + x_sol * (mw[idx_sol] / mw_sol) * (dielc[idx_sol] - eps_sol_w)
        deps[idx_ion] = dielc[idx_ion]
        return deps

    if rule == 4:
        # Figiel 2025 dielectric mixing-rule derivative.
        if mw is None or z is None:
            raise ValueError("dielc_rule=4 requires mw and z.")
        idx_sol = np.where(np.abs(z) <= 1e-12)[0]
        idx_ion = np.where(np.abs(z) > 1e-12)[0]
        if idx_sol.size == 0:
            raise ValueError("dielc_rule=4 requires at least one solvent species (z=0).")
        mw_sol = float(np.dot(x[idx_sol], mw[idx_sol]))
        if mw_sol <= 0:
            raise ValueError("Solvent molecular-weight denominator must be positive for dielc_rule=4.")
        eps_sf = float(np.dot(x[idx_sol] * mw[idx_sol], dielc[idx_sol]) / mw_sol)
        x_ion = float(np.sum(x[idx_ion])) if idx_ion.size > 0 else 0.0
        den = 1.0 + alpha * x_ion
        deps = np.zeros(len(x), dtype=float)
        deps[idx_sol] = (1.0 / den) * (mw[idx_sol] / mw_sol) * (dielc[idx_sol] - eps_sf)
        deps[idx_ion] = -alpha * eps_sf / (den * den)
        return deps

    if rule == 5:
        # Same eps_mix as rule 4, but force rule-1 differential.
        return np.asarray(dielc, dtype=float)

    if rule == 6:
        # Solvent-only differential: solvents use their dielectric values, ions are zero.
        if z is None:
            raise ValueError("dielc_rule=6 requires z.")
        idx_sol = np.where(np.abs(z) <= 1e-12)[0]
        idx_ion = np.where(np.abs(z) > 1e-12)[0]
        if idx_sol.size == 0:
            raise ValueError("dielc_rule=6 requires at least one solvent species (z=0).")
        deps = np.zeros(len(x), dtype=float)
        deps[idx_sol] = dielc[idx_sol]
        deps[idx_ion] = 0.0
        return deps

    raise ValueError("Unknown dielc_rule: {}. Supported rules: 0, 1, 2, 3, 4, 5, 6.".format(rule))

def molality_to_molefraction(molality, species=None, solvent=None, basis_mass_kg=1.0):
    """
    Convert salt molality (mol/kg solvent) into mole fractions for a solvent + cation + anion system.

    Parameters
    ----------
    molality : float
        Molality of the salt in mol per kg of solvent.
    species : list[str]
        Ordered list of species names to align the returned mole-fraction array.
        Must include exactly one solvent (neutral), one cation (name ends with '+'),
        and one anion (name ends with '-').
    solvent : str, optional
        Name of the solvent species. If omitted, the lone neutral (non-ionic) species in `species`
        is used automatically. If multiple neutrals are present, raises an error.
    basis_mass_kg : float, optional
        Mass of solvent (kg) used as the molality basis. Default is 1 kg solvent.

    Returns
    -------
    np.ndarray
        Mole fractions aligned with `species` order.
    """
    if species is None:
        raise ValueError("`species` must be provided to align mole fractions.")
    molality = float(molality)
    basis_mass_kg = float(basis_mass_kg)

    # Identify charged species by name suffix; fall back to charge sign if available.
    cations = [sp for sp in species if sp.endswith("+")]
    anions = [sp for sp in species if sp.endswith("-")]

    if len(cations) != 1 or len(anions) != 1:
        # Fallback using charge sign from property table in case names lack +/- suffix
        cations = [sp for sp in species if pcsaft_prop.get(sp, {}).get("z", 0) > 0]
        anions = [sp for sp in species if pcsaft_prop.get(sp, {}).get("z", 0) < 0]
    if len(cations) != 1 or len(anions) != 1:
        raise ValueError("Expected exactly one cation and one anion in `species`.")

    cation, anion = cations[0], anions[0]

    if solvent is None:
        neutrals = [sp for sp in species if (not sp.endswith("+") and not sp.endswith("-") and pcsaft_prop.get(sp, {}).get("z", 0) == 0)]
        if len(neutrals) != 1:
            raise ValueError("Expected exactly one neutral solvent species when `solvent` is not provided.")
        solvent = neutrals[0]
    elif solvent not in species:
        raise ValueError(f"Solvent '{solvent}' not found in provided `species` list.")

    z_cat = pcsaft_prop[cation]["z"]
    z_an = pcsaft_prop[anion]["z"]
    if z_cat <= 0 or z_an >= 0:
        raise ValueError("Charges for cation/anion must be positive/negative respectively.")

    # Stoichiometric coefficients per formula unit (charge-balanced).
    z_cat_abs, z_an_abs = int(round(abs(z_cat))), int(round(abs(z_an)))
    gcd_z = math.gcd(z_cat_abs, z_an_abs)
    v_cat = z_an_abs // gcd_z  # moles of cation per formula unit
    v_an = z_cat_abs // gcd_z  # moles of anion per formula unit

    mw_solvent = pcsaft_prop[solvent]["MW"]
    n_solvent = basis_mass_kg / mw_solvent
    n_cation = molality * basis_mass_kg * v_cat
    n_anion = molality * basis_mass_kg * v_an

    n_totals = {sp: 0.0 for sp in species}
    n_totals[solvent] += n_solvent
    n_totals[cation] += n_cation
    n_totals[anion] += n_anion

    total_moles = sum(n_totals.values())
    if total_moles <= 0:
        raise ValueError("Total moles computed as zero or negative; check inputs.")

    return np.array([n_totals[sp] / total_moles for sp in species])

