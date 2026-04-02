"""Public API surface for the object-oriented pcsaft runtime."""

from __future__ import annotations

import warnings
from typing import Iterable

import numpy as np

from .oop import (
    FlashResult,
    MultiphaseLLEResult,
    PCSAFTMixture,
    PCSAFTState,
    PhaseResult,
    VaporizationResult,
    mixture_from_legacy_params,
)
from .pcsaft import (
    InputError,
    SolutionError,
    aly_lee,
    dielc_water,
    flashPQ as _legacy_flashPQ,
    flashTQ as _legacy_flashTQ,
    pcsaft_Hvap as _legacy_Hvap,
    pcsaft_dielc_eval as _legacy_dielc_eval,
    pcsaft_multiphase_lle as _legacy_multiphase_lle,
)

__all__ = [
    "FlashResult",
    "InputError",
    "MultiphaseLLEResult",
    "PCSAFTMixture",
    "PCSAFTState",
    "PhaseResult",
    "SolutionError",
    "VaporizationResult",
    "aly_lee",
    "dielc_water",
    "flashPQ",
    "flashTQ",
    "pcsaft_Hvap",
    "pcsaft_Z",
    "pcsaft_ares",
    "pcsaft_cp",
    "pcsaft_dadt",
    "pcsaft_den",
    "pcsaft_dielc_eval",
    "pcsaft_fugcoef",
    "pcsaft_gres",
    "pcsaft_gsolv",
    "pcsaft_hres",
    "pcsaft_lnfugcoef",
    "pcsaft_lnfugcoef_terms",
    "pcsaft_miac",
    "pcsaft_miac_m",
    "pcsaft_mures",
    "pcsaft_multiphase_lle",
    "pcsaft_osmoticC",
    "pcsaft_p",
    "pcsaft_sres",
]


def _warn_legacy(name: str) -> None:
    warnings.warn(
        f"{name} is deprecated; use PCSAFTMixture/PCSAFTState instead.",
        DeprecationWarning,
        stacklevel=2,
    )


def _state_from_legacy(
    t: float,
    x,
    params: dict,
    *,
    P: float | None = None,
    rho: float | None = None,
    phase: str = "liq",
    species: Iterable[str] | None = None,
) -> PCSAFTState:
    mixture = mixture_from_legacy_params(params, species=species)
    return mixture.state(T=t, x=np.asarray(x, dtype=float), P=P, rho=rho, phase=phase)


def pcsaft_p(t, rho, x, params):
    _warn_legacy("pcsaft_p")
    state = _state_from_legacy(t, x, params, rho=rho)
    return state.pressure()


def pcsaft_den(t, p, x, params, phase="liq"):
    _warn_legacy("pcsaft_den")
    state = _state_from_legacy(t, x, params, P=p, phase=phase)
    return state.density()


def pcsaft_Z(t, rho, x, params):
    _warn_legacy("pcsaft_Z")
    state = _state_from_legacy(t, x, params, rho=rho)
    return state.Z()


def pcsaft_ares(t, rho, x, params):
    _warn_legacy("pcsaft_ares")
    state = _state_from_legacy(t, x, params, rho=rho)
    return state.ares()


def pcsaft_dadt(t, rho, x, params):
    _warn_legacy("pcsaft_dadt")
    state = _state_from_legacy(t, x, params, rho=rho)
    return state.dadt()


def pcsaft_hres(t, rho, x, params):
    _warn_legacy("pcsaft_hres")
    state = _state_from_legacy(t, x, params, rho=rho)
    return state.hres()


def pcsaft_sres(t, rho, x, params):
    _warn_legacy("pcsaft_sres")
    state = _state_from_legacy(t, x, params, rho=rho)
    return state.sres()


def pcsaft_gres(t, rho, x, params):
    _warn_legacy("pcsaft_gres")
    state = _state_from_legacy(t, x, params, rho=rho)
    return state.gres()


def pcsaft_lnfugcoef(t, rho, x, params):
    _warn_legacy("pcsaft_lnfugcoef")
    state = _state_from_legacy(t, x, params, rho=rho)
    return state.lnfugcoef()


def pcsaft_fugcoef(t, rho, x, params):
    _warn_legacy("pcsaft_fugcoef")
    state = _state_from_legacy(t, x, params, rho=rho)
    return state.fugcoef()


def pcsaft_lnfugcoef_terms(t, rho, x, params):
    _warn_legacy("pcsaft_lnfugcoef_terms")
    state = _state_from_legacy(t, x, params, rho=rho)
    return state.lnfugcoef_terms()


def pcsaft_cp(t, rho, aly_lee_params, x, params):
    _warn_legacy("pcsaft_cp")
    state = _state_from_legacy(t, x, params, rho=rho)
    return state.cp(aly_lee_params)


def pcsaft_dielc_eval(x, params):
    _warn_legacy("pcsaft_dielc_eval")
    return _legacy_dielc_eval(x, params)


def pcsaft_osmoticC(t, rho, x, params):
    _warn_legacy("pcsaft_osmoticC")
    state = _state_from_legacy(t, x, params, rho=rho)
    return state.osmoticC()


def pcsaft_miac_m(t, rho, x, params, species=None):
    _warn_legacy("pcsaft_miac_m")
    state = _state_from_legacy(t, x, params, rho=rho, species=species)
    return state.miac_m(species=species)


def pcsaft_miac(t, rho, x, params, species=None):
    _warn_legacy("pcsaft_miac")
    state = _state_from_legacy(t, x, params, rho=rho, species=species)
    return state.miac(species=species)


def pcsaft_gsolv(t, rho, x, params, species=None):
    _warn_legacy("pcsaft_gsolv")
    state = _state_from_legacy(t, x, params, rho=rho, species=species)
    return state.gsolv(species=species)


def flashTQ(t, q, x, params, p_guess=None):
    _warn_legacy("flashTQ")
    if p_guess is not None:
        return _legacy_flashTQ(t, q, x, params, p_guess=p_guess)
    return _legacy_flashTQ(t, q, x, params)


def flashPQ(p, q, x, params, t_guess=None):
    _warn_legacy("flashPQ")
    if t_guess is not None:
        return _legacy_flashPQ(p, q, x, params, t_guess=t_guess)
    return _legacy_flashPQ(p, q, x, params)


def pcsaft_Hvap(t, x, params, p_guess=None):
    _warn_legacy("pcsaft_Hvap")
    if p_guess is not None:
        return _legacy_Hvap(t, x, params, p_guess=p_guess)
    return _legacy_Hvap(t, x, params)


def pcsaft_multiphase_lle(t, p, z_feed, params, species, options=None):
    _warn_legacy("pcsaft_multiphase_lle")
    return _legacy_multiphase_lle(t, p, z_feed, params, species, options=options)


def pcsaft_mures(t, rho, x, params):
    state = _state_from_legacy(t, x, params, rho=rho)
    return state.mures()
