from __future__ import annotations

from functools import lru_cache

import numpy as np

from scripts._env import require_pcsaft_install

require_pcsaft_install()

from pcsaft.parameters import get_prop_dict
from pcsaft import pcsaft_den, pcsaft_gsolv, pcsaft_lnfugcoef_terms, pcsaft_p


R_GAS = 8.31446261815324
T_REF = 298.15
P_REF = 1.0e5
EPS = 1.0e-8
EPS_INF = 1.0e-12

SOLVENT_SPECIES = {
    "water": "Water",
    "methanol": "Methanol",
    "ethanol": "Ethanol",
}

VARIANT_DATASET = {
    "advanced": "2020_Bulow",
    "revised": "2014_Held",
}
CONTRIBUTION_KEYS = {
    "hc": "lnfugcoef_hc",
    "disp": "lnfugcoef_disp",
    "polar": "lnfugcoef_polar",
    "assoc": "lnfugcoef_assoc",
    "dh": "lnfugcoef_ion",
    "born": "lnfugcoef_born",
}
MU_CONTRIBUTION_KEYS = {
    "hc": "mu_hc",
    "disp": "mu_disp",
    "polar": "mu_polar",
    "assoc": "mu_assoc",
    "dh": "mu_ion",
    "born": "mu_born",
}


def _user_options_for_overrides(d_born_mode: int | None = None) -> dict:
    if d_born_mode is None:
        return {}
    return {"elec_model": {"born_model": {"d_Born_mode": int(d_born_mode)}}}


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
def gsolv_ion(variant: str, ion: str, solvent: str, d_born_mode: int | None = None) -> float:
    dataset_name = VARIANT_DATASET[variant]
    species = _species_for_ion(ion, solvent)
    x = np.asarray([EPS, EPS, 1.0 - 2.0 * EPS], dtype=float)
    params = get_prop_dict(dataset_name, species, x, T_REF, user_options=_user_options_for_overrides(d_born_mode))
    rho = pcsaft_den(T_REF, P_REF, x, params, phase="liq")
    values = pcsaft_gsolv(T_REF, rho, x, params, species=species)
    return float(values[ion]) / 1000.0


@lru_cache(maxsize=None)
def _infinite_dilution_terms(variant: str, ion: str, solvent: str, d_born_mode: int | None = None):
    dataset_name = VARIANT_DATASET[variant]
    species = _species_for_ion(ion, solvent)
    x = np.asarray([EPS, EPS, 1.0 - 2.0 * EPS], dtype=float)
    params = get_prop_dict(dataset_name, species, x, T_REF, user_options=_user_options_for_overrides(d_born_mode))
    rho = pcsaft_den(T_REF, P_REF, x, params, phase="liq")

    z = np.asarray(params.get("z", []), dtype=float)
    idx_ion = np.where(np.abs(z) > 1.0e-12)[0]
    idx_solv = np.where(np.abs(z) <= 1.0e-12)[0]
    x_ref = x.copy()
    x_ref[idx_ion] = 0.0
    solv_sum = float(np.sum(x_ref[idx_solv]))
    if solv_sum > 0.0:
        x_ref[idx_solv] = x_ref[idx_solv] / solv_sum
    else:
        x_ref[idx_solv] = 1.0 / len(idx_solv)

    p_ref = pcsaft_p(T_REF, rho, x_ref, params)
    x_inf = x_ref.copy()
    ion_idx = species.index(ion)
    x_inf[ion_idx] = EPS_INF
    x_inf /= np.sum(x_inf)
    phase = "vap" if rho < 900.0 else "liq"
    rho_inf = pcsaft_den(T_REF, p_ref, x_inf, params, phase=phase)
    terms = pcsaft_lnfugcoef_terms(T_REF, rho_inf, x_inf, params)
    return terms, ion_idx


@lru_cache(maxsize=None)
def contribution_breakdown(
    variant: str,
    ion: str,
    solvent: str,
    basis: str = "lnfug",
    d_born_mode: int | None = None,
) -> dict[str, float]:
    terms, idx = _infinite_dilution_terms(variant, ion, solvent, d_born_mode=d_born_mode)
    if basis == "lnfug":
        key_map = CONTRIBUTION_KEYS
    elif basis == "mu":
        key_map = MU_CONTRIBUTION_KEYS
    else:
        raise ValueError(f"Unsupported contribution basis '{basis}'.")

    out = {
        key: float(R_GAS * T_REF * terms[term_key][idx] / 1000.0)
        for key, term_key in key_map.items()
    }
    total_key = "lnfugcoef_total" if basis == "lnfug" else "mu_total"
    out["total"] = float(R_GAS * T_REF * terms[total_key][idx] / 1000.0)
    return out


def transfer_total(variant: str, ion: str, organic_solvent: str, d_born_mode: int | None = None) -> float:
    return gsolv_ion(variant, ion, organic_solvent, d_born_mode=d_born_mode) - gsolv_ion(
        variant, ion, "water", d_born_mode=d_born_mode
    )


def transfer_breakdown(
    variant: str,
    ion: str,
    organic_solvent: str,
    basis: str = "lnfug",
    d_born_mode: int | None = None,
) -> dict[str, float]:
    organic = contribution_breakdown(variant, ion, organic_solvent, basis=basis, d_born_mode=d_born_mode)
    water = contribution_breakdown(variant, ion, "water", basis=basis, d_born_mode=d_born_mode)
    return {key: organic[key] - water[key] for key in organic}
