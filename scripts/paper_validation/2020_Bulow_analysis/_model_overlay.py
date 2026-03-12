from __future__ import annotations

from functools import lru_cache

import numpy as np

from data.epcsaft_properties import get_prop_dict
from pcsaft import pcsaft_den, pcsaft_gsolv, pcsaft_lnfugcoef_terms


R_GAS = 8.31446261815324
T_REF = 298.15
P_REF = 1.0e5
EPS = 1.0e-8

SOLVENT_SPECIES = {
    "water": "Water",
    "methanol": "Methanol",
    "ethanol": "Ethanol",
}

VARIANT_DATASET = {
    "advanced": "bulow_2020",
    "revised": "held_2014",
}
CONTRIBUTION_KEYS = {
    "hc": "lnfugcoef_hc",
    "disp": "lnfugcoef_disp",
    "polar": "lnfugcoef_polar",
    "assoc": "lnfugcoef_assoc",
    "dh": "lnfugcoef_ion",
    "born": "lnfugcoef_born",
}


def _species_for_ion(ion: str, solvent: str) -> list[str]:
    solvent_species = SOLVENT_SPECIES[solvent]
    if ion in {"Li+", "Na+", "K+"}:
        return [ion, "Cl-", solvent_species]
    if ion == "F-":
        return ["Na+", "F-", solvent_species]
    if ion in {"Cl-", "Br-", "I-"}:
        return ["Na+", ion, solvent_species]
    raise KeyError(f"Unsupported ion '{ion}'.")


@lru_cache(maxsize=None)
def gsolv_ion(variant: str, ion: str, solvent: str) -> float:
    dataset_name = VARIANT_DATASET[variant]
    species = _species_for_ion(ion, solvent)
    x = np.asarray([EPS, EPS, 1.0 - 2.0 * EPS], dtype=float)
    params = get_prop_dict(dataset_name, species, x, T_REF, user_options={})
    rho = pcsaft_den(T_REF, P_REF, x, params, phase="liq")
    values = pcsaft_gsolv(T_REF, rho, x, params, species=species)
    return float(values[ion]) / 1000.0


@lru_cache(maxsize=None)
def contribution_breakdown(variant: str, ion: str, solvent: str) -> dict[str, float]:
    dataset_name = VARIANT_DATASET[variant]
    species = _species_for_ion(ion, solvent)
    x = np.asarray([EPS, EPS, 1.0 - 2.0 * EPS], dtype=float)
    params = get_prop_dict(dataset_name, species, x, T_REF, user_options={})
    rho = pcsaft_den(T_REF, P_REF, x, params, phase="liq")
    terms = pcsaft_lnfugcoef_terms(T_REF, rho, x, params)
    idx = species.index(ion)

    out = {
        key: float(R_GAS * T_REF * terms[term_key][idx] / 1000.0)
        for key, term_key in CONTRIBUTION_KEYS.items()
    }
    out["total"] = float(R_GAS * T_REF * terms["lnfugcoef_total"][idx] / 1000.0)
    return out


def transfer_total(variant: str, ion: str, organic_solvent: str) -> float:
    return gsolv_ion(variant, ion, organic_solvent) - gsolv_ion(variant, ion, "water")


def transfer_breakdown(variant: str, ion: str, organic_solvent: str) -> dict[str, float]:
    organic = contribution_breakdown(variant, ion, organic_solvent)
    water = contribution_breakdown(variant, ion, "water")
    return {key: organic[key] - water[key] for key in organic}
