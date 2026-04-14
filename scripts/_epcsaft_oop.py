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
    return state_from_params(t, x, params_or_mixture, species=species, rho=rho, phase=phase).fugacity_coefficient_terms()


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
