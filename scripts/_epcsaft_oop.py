from __future__ import annotations

import numpy as np

from epcsaft import InputError
from epcsaft import ePCSAFTMixture
from epcsaft import ePCSAFTState


def _as_array(values) -> np.ndarray:
    return np.asarray(values, dtype=float)


def as_mixture(model, species=None) -> ePCSAFTMixture:
    if isinstance(model, ePCSAFTMixture):
        return model
    return ePCSAFTMixture.from_params(model, species=species)


def state_from_params(
    t: float,
    x,
    params_or_mixture,
    *,
    species=None,
    phase: str = "liq",
    P: float | None = None,
    rho: float | None = None,
) -> ePCSAFTState:
    mixture = as_mixture(params_or_mixture, species=species)
    return mixture.state(T=float(t), x=_as_array(x), P=P, rho=rho, phase=phase)


def state_from_dataset(
    dataset_name: str,
    species,
    x,
    t: float,
    *,
    user_options: dict | None = None,
    phase: str = "liq",
    P: float | None = None,
    rho: float | None = None,
) -> ePCSAFTState:
    mixture = ePCSAFTMixture.from_dataset(dataset_name, species, _as_array(x), float(t), user_options=user_options)
    return mixture.state(T=float(t), x=_as_array(x), P=P, rho=rho, phase=phase)


def epcsaft_density(t, p, x, params_or_mixture, phase="liq", species=None):
    return float(state_from_params(t, x, params_or_mixture, species=species, P=p, phase=phase).density())


def epcsaft_pressure(t, rho, x, params_or_mixture, phase="liq", species=None):
    return float(state_from_params(t, x, params_or_mixture, species=species, rho=rho, phase=phase).pressure())


def epcsaft_compressibility_factor(t, rho, x, params_or_mixture, phase="liq", species=None):
    return float(state_from_params(t, x, params_or_mixture, species=species, rho=rho, phase=phase).compressibility_factor())


def epcsaft_residual_helmholtz(t, rho, x, params_or_mixture, phase="liq", species=None):
    return float(state_from_params(t, x, params_or_mixture, species=species, rho=rho, phase=phase).residual_helmholtz())


def epcsaft_temperature_derivative_residual_helmholtz(t, rho, x, params_or_mixture, phase="liq", species=None):
    return float(state_from_params(t, x, params_or_mixture, species=species, rho=rho, phase=phase).temperature_derivative_residual_helmholtz())


def epcsaft_composition_derivative_residual_helmholtz(t, rho, x, params_or_mixture, phase="liq", species=None):
    return state_from_params(t, x, params_or_mixture, species=species, rho=rho, phase=phase).composition_derivative_residual_helmholtz()


def epcsaft_residual_enthalpy(t, rho, x, params_or_mixture, phase="liq", species=None):
    return float(state_from_params(t, x, params_or_mixture, species=species, rho=rho, phase=phase).residual_enthalpy())


def epcsaft_residual_entropy(t, rho, x, params_or_mixture, phase="liq", species=None):
    return float(state_from_params(t, x, params_or_mixture, species=species, rho=rho, phase=phase).residual_entropy())


def epcsaft_residual_gibbs(t, rho, x, params_or_mixture, phase="liq", species=None):
    return float(state_from_params(t, x, params_or_mixture, species=species, rho=rho, phase=phase).residual_gibbs())


def epcsaft_residual_chemical_potential(t, rho, x, params_or_mixture, phase="liq", species=None):
    return state_from_params(t, x, params_or_mixture, species=species, rho=rho, phase=phase).residual_chemical_potential()


def epcsaft_fugacity_coefficient(t, rho, x, params_or_mixture, phase="liq", species=None, natural_log=False):
    return state_from_params(t, x, params_or_mixture, species=species, rho=rho, phase=phase).fugacity_coefficient(natural_log=natural_log)


def epcsaft_fugacity_coefficient_terms(t, rho, x, params_or_mixture, phase="liq", species=None):
    state = state_from_params(t, x, params_or_mixture, species=species, rho=rho, phase=phase)
    fug = state.fugacity_coefficient(return_contribution_terms=True)
    mu = state.residual_chemical_potential(return_contribution_terms=True)
    dadx = state.composition_derivative_residual_helmholtz()
    terms = {
        "mu_hc": np.asarray(mu["terms"]["hc"], dtype=float),
        "mu_disp": np.asarray(mu["terms"]["disp"], dtype=float),
        "mu_assoc": np.asarray(mu["terms"]["assoc"], dtype=float),
        "mu_ion": np.asarray(mu["terms"]["ion"], dtype=float),
        "mu_born": np.asarray(mu["terms"]["born"], dtype=float),
        "mu_total": np.asarray(mu["total"], dtype=float),
        "lnfugcoef_hc": np.asarray(fug["terms"]["hc"], dtype=float),
        "lnfugcoef_disp": np.asarray(fug["terms"]["disp"], dtype=float),
        "lnfugcoef_assoc": np.asarray(fug["terms"]["assoc"], dtype=float),
        "lnfugcoef_ion": np.asarray(fug["terms"]["ion"], dtype=float),
        "lnfugcoef_born": np.asarray(fug["terms"]["born"], dtype=float),
        "lnfugcoef_total": np.asarray(fug["terms_total_natural_log"], dtype=float),
        "dadx_hc": np.asarray(dadx["terms"]["hc"], dtype=float),
        "dadx_disp": np.asarray(dadx["terms"]["disp"], dtype=float),
        "dadx_assoc": np.asarray(dadx["terms"]["assoc"], dtype=float),
        "dadx_ion": np.asarray(dadx["terms"]["ion"], dtype=float),
        "dadx_born": np.asarray(dadx["terms"]["born"], dtype=float),
        "a_hc": float(dadx["ares_terms"]["hc"]),
        "a_disp": float(dadx["ares_terms"]["disp"]),
        "a_assoc": float(dadx["ares_terms"]["assoc"]),
        "a_ion": float(dadx["ares_terms"]["ion"]),
        "a_born": float(dadx["ares_terms"]["born"]),
        "sum_x_dadx_hc": float(dadx["sum_x_terms"]["hc"]),
        "sum_x_dadx_disp": float(dadx["sum_x_terms"]["disp"]),
        "sum_x_dadx_assoc": float(dadx["sum_x_terms"]["assoc"]),
        "sum_x_dadx_ion": float(dadx["sum_x_terms"]["ion"]),
        "sum_x_dadx_born": float(dadx["sum_x_terms"]["born"]),
        "z_raw_hc": float(dadx["z_raw_terms"]["hc"]),
        "z_raw_disp": float(dadx["z_raw_terms"]["disp"]),
        "z_raw_assoc": float(dadx["z_raw_terms"]["assoc"]),
        "z_raw_ion": float(dadx["z_raw_terms"]["ion"]),
        "z_raw_born": float(dadx["z_raw_terms"]["born"]),
        "z_hc": float(dadx["z_terms"]["hc"]),
        "z_disp": float(dadx["z_terms"]["disp"]),
        "z_assoc": float(dadx["z_terms"]["assoc"]),
        "z_ion": float(dadx["z_terms"]["ion"]),
        "z_born": float(dadx["z_terms"]["born"]),
        "z_total": float(dadx["z_total"]),
    }
    zero_vector = np.zeros_like(np.asarray(terms["lnfugcoef_total"], dtype=float))
    terms["mu_polar"] = zero_vector.copy()
    terms["lnfugcoef_polar"] = zero_vector.copy()
    terms["dadx_polar"] = zero_vector.copy()
    terms["a_polar"] = 0.0
    terms["sum_x_dadx_polar"] = 0.0
    terms["z_raw_polar"] = 0.0
    terms["z_polar"] = 0.0
    return terms


def epcsaft_relative_permittivity(x, params_or_mixture, species=None, t: float = 298.15, phase: str = "liq"):
    return state_from_params(t, x, params_or_mixture, species=species, P=1.0e5, phase=phase).relative_permittivity()


def epcsaft_activity_coefficient(
    t,
    rho,
    x,
    params_or_mixture,
    species=None,
    solvent=None,
    phase="liq",
    mean_ionic_form=False,
    basis="mole",
):
    return state_from_params(t, x, params_or_mixture, species=species, rho=rho, phase=phase).activity_coefficient(
        species=species,
        solvent=solvent,
        mean_ionic_form=mean_ionic_form,
        basis=basis,
    )


def epcsaft_solvation_free_energy(t, rho, x, params_or_mixture, species=None, phase="liq"):
    return state_from_params(t, x, params_or_mixture, species=species, rho=rho, phase=phase).solvation_free_energy(species=species)


def epcsaft_osmotic_coefficient(t, rho, x, params_or_mixture, species=None, phase="liq"):
    return state_from_params(t, x, params_or_mixture, species=species, rho=rho, phase=phase).osmotic_coefficient()
