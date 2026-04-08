from __future__ import annotations

from dataclasses import asdict
from typing import Any

import numpy as np

from pcsaft import InputError
from pcsaft import PCSAFTMixture
from pcsaft import PCSAFTState


def _as_array(values) -> np.ndarray:
    return np.asarray(values, dtype=float)


def as_mixture(model, species=None) -> PCSAFTMixture:
    if isinstance(model, PCSAFTMixture):
        return model
    return PCSAFTMixture.from_params(model, species=species)


def state_from_params(
    t: float,
    x,
    params_or_mixture,
    *,
    species=None,
    phase: str = "liq",
    P: float | None = None,
    rho: float | None = None,
) -> PCSAFTState:
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
) -> PCSAFTState:
    mixture = PCSAFTMixture.from_dataset(dataset_name, species, _as_array(x), float(t), user_options=user_options)
    return mixture.state(T=float(t), x=_as_array(x), P=P, rho=rho, phase=phase)


def pcsaft_den(t, p, x, params_or_mixture, phase="liq", species=None):
    return float(state_from_params(t, x, params_or_mixture, species=species, P=p, phase=phase).density())


def pcsaft_p(t, rho, x, params_or_mixture, phase="liq", species=None):
    return float(state_from_params(t, x, params_or_mixture, species=species, rho=rho, phase=phase).pressure())


def pcsaft_Z(t, rho, x, params_or_mixture, phase="liq", species=None):
    return float(state_from_params(t, x, params_or_mixture, species=species, rho=rho, phase=phase).Z())


def pcsaft_ares(t, rho, x, params_or_mixture, phase="liq", species=None):
    return float(state_from_params(t, x, params_or_mixture, species=species, rho=rho, phase=phase).a_res())


def pcsaft_dadt(t, rho, x, params_or_mixture, phase="liq", species=None):
    return float(state_from_params(t, x, params_or_mixture, species=species, rho=rho, phase=phase).dadt())


def pcsaft_hres(t, rho, x, params_or_mixture, phase="liq", species=None):
    return float(state_from_params(t, x, params_or_mixture, species=species, rho=rho, phase=phase).h_res())


def pcsaft_sres(t, rho, x, params_or_mixture, phase="liq", species=None):
    return float(state_from_params(t, x, params_or_mixture, species=species, rho=rho, phase=phase).s_res())


def pcsaft_gres(t, rho, x, params_or_mixture, phase="liq", species=None):
    return float(state_from_params(t, x, params_or_mixture, species=species, rho=rho, phase=phase).g_res())


def pcsaft_fugcoef(t, rho, x, params_or_mixture, phase="liq", species=None):
    return np.exp(state_from_params(t, x, params_or_mixture, species=species, rho=rho, phase=phase).lnfugcoef())


def pcsaft_lnfugcoef_terms(t, rho, x, params_or_mixture, phase="liq", species=None):
    return state_from_params(t, x, params_or_mixture, species=species, rho=rho, phase=phase).lnfugcoef_terms()


def pcsaft_dielc_eval(x, params_or_mixture, species=None, t: float = 298.15, phase: str = "liq"):
    return state_from_params(t, x, params_or_mixture, species=species, P=1.0e5, phase=phase).dielectric_eval()


def pcsaft_miac_m(t, rho, x, params_or_mixture, species=None, phase="liq"):
    return state_from_params(t, x, params_or_mixture, species=species, rho=rho, phase=phase).miac_m(species=species)


def pcsaft_miac(t, rho, x, params_or_mixture, species=None, phase="liq"):
    return state_from_params(t, x, params_or_mixture, species=species, rho=rho, phase=phase).miac(species=species)


def pcsaft_gsolv(t, rho, x, params_or_mixture, species=None, phase="liq"):
    return state_from_params(t, x, params_or_mixture, species=species, rho=rho, phase=phase).gsolv(species=species)


def flashTQ(t, q, x, params_or_mixture, p_guess=None, species=None):
    state = state_from_params(t, x, params_or_mixture, species=species, P=p_guess if p_guess is not None else 1.0e5, phase="liq")
    result = state.flashTQ(q, p_guess=p_guess)
    return result.value, result.phases[0].x, result.phases[1].x


def flashPQ(p, q, x, params_or_mixture, t_guess=None, species=None):
    state = state_from_params(t_guess if t_guess is not None else 298.15, x, params_or_mixture, species=species, P=p, phase="liq")
    result = state.flashPQ(p, q, t_guess=t_guess)
    return result.value, result.phases[0].x, result.phases[1].x


def pcsaft_Hvap(t, x, params_or_mixture, p_guess=None, species=None):
    state = state_from_params(t, x, params_or_mixture, species=species, P=p_guess if p_guess is not None else 1.0e5, phase="liq")
    result = state.Hvap(p_guess=p_guess)
    return result.value, result.pressure


def pcsaft_cp(t, rho, x, params_or_mixture, aly_lee_params, phase="liq", species=None):
    return float(state_from_params(t, x, params_or_mixture, species=species, rho=rho, phase=phase).cp(aly_lee_params))


def _phase_result_to_dict(phase) -> dict[str, Any]:
    return {
        "beta": float(phase.beta),
        "x": np.asarray(phase.x, dtype=float),
        "rho": float(phase.rho),
        "lnfugcoef": np.asarray(phase.lnfugcoef, dtype=float),
        "lnfug": np.asarray(phase.lnfug, dtype=float),
    }

