import numpy as np
import types
import math
import json
import copy
import warnings
from pathlib import Path

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
        'z': 0., 'dielc': 32., 'd_born': 0., 'f_solv': 1.0,
    },
    'MEA-4C': {
        'MW': 61.08e-3,  # kg/mol
        'm': 4.5208, 's': 2.6574, 'e': 237.6864,
        'e_assoc': 989.8984, 'vol_a': 0.187533, 'assoc_scheme': '4C',
        'dipm': 0., 'dip_num': 1,
        'z': 0., 'dielc': 0., 'd_born': 0., 'f_solv': 1.0,
    },
    'MEA': {
        'MW': 61.08e-3,  # kg/mol
        'm': {'2B': 3.0353, '4C': 4.5208},
        's': {'2B': 3.0435, '4C': 2.6574},
        'e': {'2B': 277.174, '4C': 237.6864},
        'e_assoc': {'2B': 2586.3, '4C': 989.8984},
        'vol_a': {'2B': 0.037470, '4C': 0.187533},
        'assoc_scheme': {'2B': '2B', '4C': '4C'},
        'dipm': 0., 'dip_num': 1,
        'z': 0., 'dielc': {'2B': 32., '4C': 0.}, 'd_born': 0., 'f_solv': 1.0,
    },
    'H2O-2B-MEA': {
        'MW': 18.01528e-3,  # kg/mol
        'm': 1.9599, 's': 2.362, 'e': 279.42,
        'e_assoc': 2059.28, 'vol_a': 0.1750, 'assoc_scheme': '2B',
        'dipm': 0., 'dip_num': 1,
        'z': 0.,
        'dielc': 78.09,
    },
    'H2O-4C-MEA': {
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
    'H2O-Salt-2001': {
        'MW': 18.01528e-3,  # kg/mol
        'm': 1.09528, 's': 2.88980, 'e': 365.956,
        'e_assoc': 2515.6706, 'vol_a': 0.0348679836, 'assoc_scheme': '2B',
        'dipm': 0., 'dip_num': 1,
        'z': 0.,
        # 'dielc': lambda T: -105.2*np.log(T) + 677.480,
        'dielc': 78.09,
        'f_solv': 1.5,
        'd_born': 0,
    },

    'H2O-Salt-': {
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
    # ("MEA-2B", "H2O-2B-MEA"): -0.0420, # Baygi 2015
    ("MEA-2B", "H2O-2B-MEA"): 0.250,  # Baygi 2015
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

k_ij_dict[("H2O-2B-MEA", "MEAH+")] = k_ij_dict[("MEA-2B", "H2O-2B-MEA")]
k_ij_dict[("H2O-2B-MEA", "MEACOO-")] = k_ij_dict[("MEA-2B", "H2O-2B-MEA")]

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

CATALOG_PATH = Path(__file__).resolve().parent / "epc-saft_parameters" / "epc-saft_parameters.json"
_CATALOG_CACHE = None

_RULE_ALIASES = {
    "constant": 0,
    "linear-mixing-mole": 1,
    "linear-mixing-weight": 2,
    "combined": 3,
    "empirical": 4,
    "rule4": 4,
    "rule5": 5,
    "rule6": 6,
}
_DIELC_DIFF_RULE_ALIASES = {
    "same": "same",
    "rule1": "rule1",
    "rule1-override": "rule1",
    "constant_saltfree": "constant_saltfree",
}
_DIFFC_MODE_ALIASES = {"analytic": 0, "numeric": 1}
_BORN_DIFF_MODE_ALIASES = {"analytic": 0, "numeric": 1, "eq133": 2, "no_dielc_dep": 3}
_BORN_RADIUS_MODE_ALIASES = {
    "sigma": 1,
    "sigma_const": 2,
    "sigma_temp": 3,
    "dborn": 4,
    "dborn_delta": 5,
}

_CANONICAL_ELEC_MODEL = {
    "born_contrib": True,
    "ssm_ds": False,
    "dielc_rule": "linear-mixing-mole",
    "dielc_diff_rule": "same",
    "dielc_diff_mode": "analytic",
    "born_radius_model": 1,
    "born_diff_model": "analytic",
    "born_diff_options": {
        "include_sum_term": True,
        "include_dielc_conc_dep": True,
        "include_delta_d_i_conc_dep": True,
    },
    "eps_r_bulk": "mix",
    "bjeruum_treatment": False,
}

_FALLBACK_ELEC_PRESETS = {
    "legacy_default": {
        "parameter_set_key": "default",
        "component_set_key": {"default": "default", "H2O": "salt_2025"},
        "kij_set": "2025",
        "model": copy.deepcopy(_CANONICAL_ELEC_MODEL),
    },
    "2005": {
        "parameter_set_key": "2005",
        "component_set_key": {"default": "2005", "H2O": "salt_2005"},
        "kij_set": "2005",
        "model": {
            "born_contrib": False,
            "ssm_ds": False,
            "dielc_rule": "constant",
            "dielc_diff_rule": "same",
            "dielc_diff_mode": "analytic",
            "born_diff_model": "analytic",
            "born_diff_options": {
                "include_sum_term": False,
                "include_dielc_conc_dep": False,
                "include_delta_d_i_conc_dep": False,
            },
            "eps_r_bulk": "mix",
            "bjeruum_treatment": False,
        },
    },
    "2008": {
        "parameter_set_key": "2008",
        "component_set_key": {"default": "2008", "H2O": "salt_2008"},
        "kij_set": "2008",
        "model": {
            "born_contrib": False,
            "ssm_ds": False,
            "dielc_rule": "constant",
            "dielc_diff_rule": "same",
            "dielc_diff_mode": "analytic",
            "born_diff_model": "analytic",
            "born_diff_options": {
                "include_sum_term": False,
                "include_dielc_conc_dep": False,
                "include_delta_d_i_conc_dep": False,
            },
            "eps_r_bulk": "mix",
            "bjeruum_treatment": False,
        },
    },
    "2014_s1": {
        "parameter_set_key": "2014_s1",
        "component_set_key": {"default": "2014_s1", "H2O": "salt_2014_s1"},
        "kij_set": "2014_s1",
        "model": {
            "born_contrib": False,
            "ssm_ds": False,
            "dielc_rule": "constant",
            "dielc_diff_rule": "same",
            "dielc_diff_mode": "analytic",
            "born_diff_model": "analytic",
            "born_diff_options": {
                "include_sum_term": False,
                "include_dielc_conc_dep": False,
                "include_delta_d_i_conc_dep": False,
            },
            "eps_r_bulk": "mix",
            "bjeruum_treatment": False,
        },
    },
    "2014_s2": {
        "parameter_set_key": "2014_s2",
        "component_set_key": {"default": "2014_s2", "H2O": "salt_2014_s2"},
        "kij_set": "2014_s2",
        "model": {
            "born_contrib": False,
            "ssm_ds": False,
            "dielc_rule": "constant",
            "dielc_diff_rule": "same",
            "dielc_diff_mode": "analytic",
            "born_diff_model": "analytic",
            "born_diff_options": {
                "include_sum_term": False,
                "include_dielc_conc_dep": False,
                "include_delta_d_i_conc_dep": False,
            },
            "eps_r_bulk": "mix",
            "bjeruum_treatment": False,
        },
    },
    "2020": {
        "parameter_set_key": "2020",
        "component_set_key": {"default": "2020", "H2O": "salt_2020"},
        "kij_set": "2020",
        "model": {
            "born_contrib": True,
            "ssm_ds": False,
            "dielc_rule": "linear-mixing-mole",
            "dielc_diff_rule": "same",
            "dielc_diff_mode": "analytic",
            "born_diff_model": "analytic",
            "born_diff_options": {
                "include_sum_term": True,
                "include_dielc_conc_dep": True,
                "include_delta_d_i_conc_dep": True,
            },
            "eps_r_bulk": "mix",
            "bjeruum_treatment": False,
        },
    },
    "2025": {
        "parameter_set_key": "2025",
        "component_set_key": {"default": "2025", "H2O": "salt_2025"},
        "kij_set": "2025",
        "model": {
            "born_contrib": True,
            "ssm_ds": True,
            "dielc_rule": "empirical",
            "dielc_diff_rule": "same",
            "dielc_diff_mode": "analytic",
            "born_diff_model": "analytic",
            "born_diff_options": {
                "include_sum_term": True,
                "include_dielc_conc_dep": True,
                "include_delta_d_i_conc_dep": True,
            },
            "eps_r_bulk": "mix",
            "bjeruum_treatment": False,
        },
    },
}

_FALLBACK_ALIASES = {
    "preset_int": {"1": "2005", "2": "2008", "3": "2014_s1", "4": "2014_s2", "5": "2020", "6": "2025"},
    "preset_str": {"2005": "2005", "2008": "2008", "2014_s1": "2014_s1", "2014_s2": "2014_s2", "2020": "2020", "2025": "2025"},
    "component_aliases": {
        "H2O-2B-Li": "H2O",
        "H2O-2B-NaCl": "H2O",
        "H2O-Salt-2001": "H2O",
        "H2O-Salt-": "H2O",
    },
}


def _load_epcsaft_catalog():
    global _CATALOG_CACHE
    if _CATALOG_CACHE is not None:
        return _CATALOG_CACHE
    if CATALOG_PATH.exists():
        with CATALOG_PATH.open("r", encoding="utf-8") as handle:
            _CATALOG_CACHE = json.load(handle)
    else:
        _CATALOG_CACHE = None
    return _CATALOG_CACHE


def _coerce_bool(value):
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, np.integer)):
        return bool(int(value))
    if isinstance(value, str):
        v = value.strip().lower()
        if v in {"1", "true", "yes", "on"}:
            return True
        if v in {"0", "false", "no", "off"}:
            return False
    raise ValueError("Could not interpret boolean value: {!r}".format(value))


def _deep_update(base, updates):
    out = copy.deepcopy(base)
    for key, value in updates.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _deep_update(out[key], value)
        else:
            out[key] = value
    return out


def _normalize_component_name(sp, catalog):
    aliases = (catalog or {}).get("aliases", {}).get("component_aliases", {})
    aliases = {**_FALLBACK_ALIASES["component_aliases"], **aliases}
    return aliases.get(sp, sp)


def _resolve_descriptor(value, T):
    if isinstance(value, dict):
        kind = value.get("type", "").strip().lower()
        if kind == "water_sigma_2008":
            return 2.7927 + (10.11 * np.exp(-0.01775 * T) - 1.417 * np.exp(-0.01146 * T))
        if kind == "linear_t":
            return float(value["a"]) * T + float(value["b"])
    return value


def _resolve_set_value(set_values, set_key, T):
    if not isinstance(set_values, dict):
        return _resolve_descriptor(set_values, T)
    candidates = [set_key]
    if not str(set_key).startswith("salt_"):
        candidates.append("salt_{}".format(set_key))
    candidates.extend(["default", "legacy_default"])
    for key in candidates:
        if key in set_values:
            return _resolve_descriptor(set_values[key], T)
    return None


def _normalize_elec_model(model):
    out = copy.deepcopy(_CANONICAL_ELEC_MODEL)
    if model is None:
        return out
    if not isinstance(model, dict):
        raise TypeError("elec_model dict expected, got {}.".format(type(model).__name__))

    normalized = {}
    for key, value in model.items():
        k = key
        if key == "SSM+DS":
            k = "ssm_ds"
        if key == "emperical":
            k = "empirical"
        if key == "born_rel_perm":
            k = "eps_r_bulk"
        normalized[k] = value

    if "born_model" in normalized and "born_contrib" not in normalized:
        normalized["born_contrib"] = _coerce_bool(normalized["born_model"])
    if "bjeruum_treatment" in normalized:
        normalized["bjeruum_treatment"] = _coerce_bool(normalized["bjeruum_treatment"])
    if "born_contrib" in normalized:
        normalized["born_contrib"] = _coerce_bool(normalized["born_contrib"])
    if "ssm_ds" in normalized:
        normalized["ssm_ds"] = _coerce_bool(normalized["ssm_ds"])
    if "born_diff_options" in normalized and not isinstance(normalized["born_diff_options"], dict):
        raise TypeError("elec_model['born_diff_options'] must be a dict.")

    out = _deep_update(out, normalized)
    out["born_diff_options"]["include_sum_term"] = _coerce_bool(out["born_diff_options"]["include_sum_term"])
    out["born_diff_options"]["include_dielc_conc_dep"] = _coerce_bool(out["born_diff_options"]["include_dielc_conc_dep"])
    out["born_diff_options"]["include_delta_d_i_conc_dep"] = _coerce_bool(out["born_diff_options"]["include_delta_d_i_conc_dep"])
    return out


def _resolve_elec_preset(elec_model):
    catalog = _load_epcsaft_catalog()
    presets = dict(_FALLBACK_ELEC_PRESETS)
    aliases = copy.deepcopy(_FALLBACK_ALIASES)
    if catalog is not None:
        presets.update(catalog.get("elec_model_presets", {}))
        aliases.update(catalog.get("aliases", {}))

    if elec_model is None:
        preset_key = "legacy_default"
        preset = presets[preset_key]
        model = _normalize_elec_model(preset.get("model"))
        return {
            "catalog": catalog,
            "preset_key": preset_key,
            "preset": preset,
            "model": model,
        }

    if isinstance(elec_model, dict):
        base_key = elec_model.get("preset") or elec_model.get("base") or "legacy_default"
        if base_key not in presets:
            raise KeyError("Unknown elec_model base preset '{}'. Available: {}.".format(base_key, sorted(presets.keys())))
        preset = presets[base_key]
        model = _normalize_elec_model(_deep_update(preset.get("model", {}), elec_model))
        return {
            "catalog": catalog,
            "preset_key": base_key,
            "preset": preset,
            "model": model,
        }

    if isinstance(elec_model, (int, np.integer)):
        preset_key = aliases.get("preset_int", {}).get(str(int(elec_model)))
        if preset_key is None:
            raise KeyError("Unknown elec_model id '{}'. Supported ids: {}.".format(elec_model, sorted(aliases.get("preset_int", {}).keys())))
    elif isinstance(elec_model, str):
        raw = elec_model.strip()
        preset_key = aliases.get("preset_str", {}).get(raw, aliases.get("preset_int", {}).get(raw))
        if preset_key is None and raw in presets:
            preset_key = raw
        if preset_key is None:
            raise KeyError("Unknown elec_model preset '{}'. Available: {}.".format(elec_model, sorted(presets.keys())))
    else:
        raise TypeError("elec_model must be int, str, dict, or None.")

    preset = presets[preset_key]
    model = _normalize_elec_model(preset.get("model"))
    return {
        "catalog": catalog,
        "preset_key": preset_key,
        "preset": preset,
        "model": model,
    }


def _warn_legacy_key(option_key):
    warnings.warn(
        "user_options['{}'] is legacy. Prefer user_options['elec_model'] / user_options['bjeruum_treatment'].".format(option_key),
        DeprecationWarning,
        stacklevel=3,
    )


def _as_rule_number(value, aliases):
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, str):
        v = value.strip().lower()
        if v.isdigit() or (v.startswith("-") and v[1:].isdigit()):
            return int(v)
        if v in aliases:
            return int(aliases[v])
    raise ValueError("Unknown rule option '{}'. Supported aliases: {}.".format(value, sorted(aliases.keys())))


def _resolve_born_eps_mode(value):
    if isinstance(value, (int, np.integer)):
        mode = int(value)
        if mode in (0, 1):
            return mode
        raise ValueError("Unknown born permittivity mode '{}'. Supported values are 0 (mix) and 1 (solvent).".format(value))
    eps_key = str(value).strip().lower()
    if eps_key in {"mix", "bulk", "eps_r_mix"}:
        return 0
    if eps_key in {"solvent", "eps_r_solvent"}:
        return 1
    raise ValueError("Unknown born permittivity selector '{}'. Supported values are 'mix' and 'solvent'.".format(value))


def _resolve_born_radius_mode(value):
    if isinstance(value, (int, np.integer)):
        mode = int(value)
    else:
        v = str(value).strip().lower()
        if v.isdigit() or (v.startswith("-") and v[1:].isdigit()):
            mode = int(v)
        elif v in _BORN_RADIUS_MODE_ALIASES:
            mode = int(_BORN_RADIUS_MODE_ALIASES[v])
        else:
            raise ValueError(
                "Unknown born_radius_model '{}'. Supported values are 1, 2, 3, 4, 5 and aliases {}.".format(
                    value, sorted(_BORN_RADIUS_MODE_ALIASES.keys())
                )
            )
    if mode < 1 or mode > 5:
        raise ValueError("Unknown born_radius_model '{}'. Supported values are 1, 2, 3, 4, 5.".format(value))
    return mode


def _flatten_model_to_runtime(model):
    dielc_rule = _as_rule_number(model["dielc_rule"], _RULE_ALIASES)
    dielc_diff_rule = str(model.get("dielc_diff_rule", "same")).strip().lower()
    dielc_diff_rule = _DIELC_DIFF_RULE_ALIASES.get(dielc_diff_rule, dielc_diff_rule)
    if dielc_rule == 4:
        if dielc_diff_rule == "rule1":
            dielc_rule = 5
        elif dielc_diff_rule == "constant_saltfree":
            dielc_rule = 6

    if model["born_contrib"]:
        if model["ssm_ds"]:
            use_deps = model["born_diff_options"]["include_dielc_conc_dep"]
            use_delta = model["born_diff_options"]["include_delta_d_i_conc_dep"]
            if use_deps and use_delta:
                born_model = 5
            elif use_deps:
                born_model = 3
            elif use_delta:
                born_model = 4
            else:
                born_model = 2
        else:
            born_model = 1
    else:
        born_model = 0

    born_diff_model = model.get("born_diff_model", "analytic")
    if born_model == 1:
        include_dielc = _coerce_bool(model["born_diff_options"].get("include_dielc_conc_dep", True))
        if not include_dielc:
            born_diff_mode = 3
        elif born_diff_model is None:
            born_diff_mode = 2 if not model["born_diff_options"]["include_sum_term"] else 0
        else:
            born_diff_mode = _as_rule_number(born_diff_model, _BORN_DIFF_MODE_ALIASES)
    else:
        born_diff_mode = _as_rule_number(born_diff_model, _BORN_DIFF_MODE_ALIASES)

    dielc_diff_mode = _as_rule_number(model["dielc_diff_mode"], _DIFFC_MODE_ALIASES)
    born_radius_model = _resolve_born_radius_mode(model.get("born_radius_model", 1))
    if born_model == 1 and born_radius_model == 5:
        raise ValueError("born_model=1 supports born_radius_model values 1, 2, 3, 4.")
    if born_model >= 2 and born_radius_model != 5:
        raise ValueError("born_model >= 2 requires born_radius_model=5.")

    dh_model = 2 if _coerce_bool(model.get("bjeruum_treatment", False)) else 1
    born_eps_mode = _resolve_born_eps_mode(model.get("eps_r_bulk", "mix"))

    return {
        "born_model": int(born_model),
        "born_radius_model": int(born_radius_model),
        "born_diff_mode": int(born_diff_mode),
        "born_eps_mode": int(born_eps_mode),
        "dielc_rule": int(dielc_rule),
        "dielc_diff_mode": int(dielc_diff_mode),
        "DH_model": int(dh_model),
        "bjeruum_treatment": bool(dh_model == 2),
    }


def _resolve_runtime_options(user_options):
    allowed = {
        "elec_model",
        "bjeruum_treatment",
        "dielc_rule",
        "dielc_diff_mode",
        "born_model",
        "born_radius_model",
        "born_diff_mode",
        "born_eps_mode",
        "born_rel_perm",
        "eps_r_bulk",
        "DH_model",
        "debug",
    }
    if user_options is None:
        user_options = {}
    unknown_options = set(user_options) - allowed
    if unknown_options:
        raise KeyError("Unknown user_options key(s): {}. Allowed keys: {}.".format(
            sorted(unknown_options), sorted(allowed)))

    resolved = _resolve_elec_preset(user_options.get("elec_model"))
    model = copy.deepcopy(resolved["model"])
    if "bjeruum_treatment" in user_options:
        model["bjeruum_treatment"] = _coerce_bool(user_options["bjeruum_treatment"])

    runtime = _flatten_model_to_runtime(model)

    for key in ("dielc_rule", "dielc_diff_mode", "born_model", "born_diff_mode", "born_eps_mode", "DH_model"):
        if key in user_options:
            _warn_legacy_key(key)
            runtime[key] = int(user_options[key])
    if "born_radius_model" in user_options:
        runtime["born_radius_model"] = _resolve_born_radius_mode(user_options["born_radius_model"])
    if "born_rel_perm" in user_options and "eps_r_bulk" in user_options:
        raise ValueError("Use one Born permittivity selector: either born_rel_perm or eps_r_bulk.")
    if "born_rel_perm" in user_options:
        runtime["born_eps_mode"] = _resolve_born_eps_mode(user_options["born_rel_perm"])
    elif "eps_r_bulk" in user_options:
        runtime["born_eps_mode"] = _resolve_born_eps_mode(user_options["eps_r_bulk"])

    if runtime["born_model"] == 1 and runtime["born_radius_model"] == 5:
        raise ValueError("born_model=1 supports born_radius_model values 1, 2, 3, 4.")
    if runtime["born_model"] >= 2 and runtime["born_radius_model"] != 5:
        raise ValueError("born_model >= 2 requires born_radius_model=5.")

    if "DH_model" in user_options and int(user_options["DH_model"]) == 2:
        runtime["bjeruum_treatment"] = True
    elif "DH_model" in user_options:
        runtime["bjeruum_treatment"] = False

    runtime["debug"] = bool(user_options.get("debug", False))
    return {
        "runtime": runtime,
        "model": model,
        "preset_key": resolved["preset_key"],
        "preset": resolved["preset"],
        "catalog": resolved["catalog"],
    }


def _build_catalog_entry(sp, T, preset_data):
    catalog = preset_data["catalog"]
    if catalog is None:
        return None
    comp_key = _normalize_component_name(sp, catalog)
    comp_data = catalog.get("component_parameters", {}).get(comp_key)
    if comp_data is None:
        return None

    component_set_map = preset_data["preset"].get("component_set_key", {})
    default_set = preset_data["preset"].get("parameter_set_key", "default")
    set_key = component_set_map.get(comp_key, component_set_map.get(sp, component_set_map.get("default", default_set)))

    out = {}
    for prop in (BASE_KEYS + OPTIONAL_KEYS):
        if prop in comp_data:
            value = _resolve_set_value(comp_data[prop], set_key, T)
            if value is not None or prop == "assoc_scheme":
                out[prop] = value
    if out:
        return out
    return None


def _resolve_species_params(sp, user_params, T, preset_data):
    base_entry = _build_catalog_entry(sp, T, preset_data)
    if base_entry is None and sp in pcsaft_prop:
        base_entry = dict(pcsaft_prop[sp])
    if user_params is not None and sp in user_params:
        merged = {}
        if base_entry is not None:
            merged.update(base_entry)
        merged.update(user_params[sp])
        return merged, 'user_params'
    if base_entry is not None:
        return base_entry, 'catalog' if preset_data["catalog"] is not None else 'default'
    raise KeyError("Species '{}' not found in user_params, parameter catalog, or default pcsaft_prop.".format(sp))


def _get_value(entry, prop, T):
    value = entry[prop]
    if isinstance(value, types.FunctionType):
        return value(T)
    return _resolve_descriptor(value, T)


def _resolve_catalog_kij(sp1, sp2, z1, z2, T, preset_data):
    catalog = preset_data["catalog"]
    if catalog is None:
        return None
    kij_set = preset_data["preset"].get("kij_set")
    if kij_set is None:
        return None
    kij_data = catalog.get("kij_parameters", {}).get(kij_set)
    if kij_data is None:
        return None

    c1 = _normalize_component_name(sp1, catalog)
    c2 = _normalize_component_name(sp2, catalog)
    key12 = "{}|{}".format(c1, c2)
    key21 = "{}|{}".format(c2, c1)
    pairs = kij_data.get("pairs", {})
    if key12 in pairs:
        return float(_resolve_descriptor(pairs[key12], T))
    if key21 in pairs:
        return float(_resolve_descriptor(pairs[key21], T))
    if z1 * z2 < 0 and abs(z1) > 1e-12 and abs(z2) > 1e-12 and "ion_pair_default" in kij_data:
        return float(kij_data["ion_pair_default"])
    if kij_data.get("default_zero", False):
        return 0.0
    return None


def get_prop_dict(species, x, T, user_params=None, user_options=None):
    """
    species: list of species names that match dictionary keys in pcsaft_prop
    T: Temperature (K) (often not used, used in calculations of temperature-dependent binary interaction parameters)
    user_params: optional dict in the form {component: {m, s, e, ...}}
    """

    preset_data = _resolve_runtime_options(user_options)
    runtime_options = preset_data["runtime"]

    prop_dic = {}
    entries = []
    for sp in species:
        entry, _ = _resolve_species_params(sp, user_params, T, preset_data)
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
            raise KeyError("Missing '{}' for species '{}' in {}.".format(prop, sp, 'user_params' if (user_params is not None and sp in user_params) else 'defaults/catalog'))
        if prop == 'assoc_scheme':
            prop_dic[prop] = prop_list
        else:
            prop_dic[prop] = np.array(prop_list)

    for prop in OPTIONAL_KEYS:
        if any(prop in entry for entry in entries):
            prop_list = []
            for entry in entries:
                prop_list.append(_get_value(entry, prop, T) if prop in entry else 0.0)
            prop_dic[prop] = np.array(prop_list)

    n = len(species)
    z_vals = np.asarray(prop_dic.get('z', np.zeros(n)), dtype=float)

    # Create the binary interaction parameter dictionary and matrix for dispersion forces
    k_ij = np.zeros((n, n))
    for i, sp1 in enumerate(species):
        for j, sp2 in enumerate(species):
            kij_val = _resolve_catalog_kij(sp1, sp2, z_vals[i], z_vals[j], T, preset_data)
            if kij_val is not None:
                k_ij[i, j] = kij_val
                continue
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

    prop_dic['elec_model'] = copy.deepcopy(preset_data["model"])
    prop_dic['elec_model_preset'] = preset_data["preset_key"]
    prop_dic['bjeruum_treatment'] = bool(runtime_options['bjeruum_treatment'])
    prop_dic['born_model'] = runtime_options['born_model']
    prop_dic['born_radius_model'] = runtime_options['born_radius_model']
    prop_dic['born_diff_mode'] = runtime_options['born_diff_mode']
    prop_dic['born_eps_mode'] = runtime_options['born_eps_mode']
    prop_dic['DH_model'] = runtime_options['DH_model']
    prop_dic['dielc_rule'] = runtime_options['dielc_rule']
    prop_dic['dielc_diff_mode'] = runtime_options['dielc_diff_mode']
    prop_dic['debug'] = bool(runtime_options['debug'])

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

    preset_data = _resolve_runtime_options(None)
    for sp in species:
        entry = None
        if user_params is not None and sp in user_params:
            entry = user_params[sp]
        else:
            entry = _build_catalog_entry(sp, 298.15, preset_data)
            if entry is None and sp in pcsaft_prop:
                entry = pcsaft_prop[sp]
        if entry is None:
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

    def _default_species_entry(sp):
        if sp in pcsaft_prop:
            return pcsaft_prop[sp]
        catalog = _load_epcsaft_catalog()
        if catalog is None:
            raise KeyError(sp)
        comp_key = _normalize_component_name(sp, catalog)
        comp_data = catalog.get("component_parameters", {}).get(comp_key)
        if comp_data is None:
            raise KeyError(sp)
        mw = _resolve_set_value(comp_data.get("MW", {}), "default", 298.15)
        z = _resolve_set_value(comp_data.get("z", {}), "default", 298.15)
        return {"MW": mw, "z": z}

    # Identify charged species by name suffix; fall back to charge sign if available.
    cations = [sp for sp in species if sp.endswith("+")]
    anions = [sp for sp in species if sp.endswith("-")]

    if len(cations) != 1 or len(anions) != 1:
        # Fallback using charge sign from property table in case names lack +/- suffix
        cations = [sp for sp in species if _default_species_entry(sp).get("z", 0) > 0]
        anions = [sp for sp in species if _default_species_entry(sp).get("z", 0) < 0]
    if len(cations) != 1 or len(anions) != 1:
        raise ValueError("Expected exactly one cation and one anion in `species`.")

    cation, anion = cations[0], anions[0]

    if solvent is None:
        neutrals = [sp for sp in species if (not sp.endswith("+") and not sp.endswith("-") and _default_species_entry(sp).get("z", 0) == 0)]
        if len(neutrals) != 1:
            raise ValueError("Expected exactly one neutral solvent species when `solvent` is not provided.")
        solvent = neutrals[0]
    elif solvent not in species:
        raise ValueError(f"Solvent '{solvent}' not found in provided `species` list.")

    z_cat = _default_species_entry(cation)["z"]
    z_an = _default_species_entry(anion)["z"]
    if z_cat <= 0 or z_an >= 0:
        raise ValueError("Charges for cation/anion must be positive/negative respectively.")

    # Stoichiometric coefficients per formula unit (charge-balanced).
    z_cat_abs, z_an_abs = int(round(abs(z_cat))), int(round(abs(z_an)))
    gcd_z = math.gcd(z_cat_abs, z_an_abs)
    v_cat = z_an_abs // gcd_z  # moles of cation per formula unit
    v_an = z_cat_abs // gcd_z  # moles of anion per formula unit

    mw_solvent = _default_species_entry(solvent)["MW"]
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

