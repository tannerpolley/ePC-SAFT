"""Object-oriented runtime facade for pcsaft.

This module keeps the compiled PC-SAFT/ePC-SAFT kernels intact and layers a
model/state API on top of them. The object layer caches resolved parameter
payloads per state so repeated property calls do not need to rebuild runtime
parameter dictionaries.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import copy
import math
from typing import Iterable, Iterator

import numpy as np

from .parameters import get_prop_dict
from .pcsaft import (
    InputError,
    SolutionError,
    aly_lee,
    dielc_water,
    flashPQ,
    flashTQ,
    pcsaft_Hvap,
    pcsaft_Z,
    pcsaft_ares,
    pcsaft_cp,
    pcsaft_dadt,
    pcsaft_den,
    pcsaft_dielc_eval,
    pcsaft_fugcoef,
    pcsaft_gres,
    pcsaft_gsolv,
    pcsaft_hres,
    pcsaft_lnfugcoef,
    pcsaft_lnfugcoef_terms,
    pcsaft_miac,
    pcsaft_miac_m,
    pcsaft_multiphase_lle,
    pcsaft_osmoticC,
    pcsaft_p,
    pcsaft_sres,
)


R_GAS = 8.31446261815324


def _as_species_tuple(species: Iterable[str]) -> tuple[str, ...]:
    species_tuple = tuple(str(item) for item in species)
    if not species_tuple:
        raise InputError("species must contain at least one component.")
    return species_tuple


def _as_float_array(values, *, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float).flatten()
    if arr.ndim != 1:
        raise InputError(f"{name} must be one-dimensional.")
    arr = np.asarray(arr, dtype=float)
    arr.setflags(write=False)
    return arr


def _copy_params(params: dict) -> dict:
    return copy.deepcopy(params)


def _validate_state_inputs(T: float, x, P, rho) -> tuple[float, np.ndarray, float | None, float | None]:
    x_arr = _as_float_array(x, name="x")
    payload = {"temperature": float(T)}
    if P is not None and rho is not None:
        raise InputError("Exactly one of P or rho must be provided.")
    if P is None and rho is None:
        raise InputError("Exactly one of P or rho must be provided.")
    if P is not None:
        payload["pressure"] = float(P)
    if rho is not None:
        payload["density"] = float(rho)
    # Reuse the legacy validation messages for temperature/pressure/density and x sum checks.
    try:
        from .pcsaft import check_input

        check_input(x_arr, payload)
    except Exception:
        raise
    return float(T), x_arr, (float(P) if P is not None else None), (float(rho) if rho is not None else None)


def _infer_species_from_params(params: dict, ncomp: int | None = None) -> tuple[str, ...]:
    for key in ("species", "component_names", "labels"):
        if key in params:
            return _as_species_tuple(params[key])
    if ncomp is None:
        for key in ("m", "s", "e", "z", "dielc", "MW"):
            if key in params:
                ncomp = len(np.asarray(params[key], dtype=float).flatten())
                break
    if ncomp is None:
        raise InputError("Could not infer species labels from params; pass species explicitly.")
    return tuple(str(i) for i in range(int(ncomp)))


def _result_array(values) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    arr.setflags(write=False)
    return arr


@dataclass(frozen=True, slots=True, eq=False)
class FlashResult:
    """Structured result for flash calculations."""

    value: float
    xl: np.ndarray
    xv: np.ndarray
    kind: str

    def to_legacy(self):
        return (float(self.value), np.asarray(self.xl, dtype=float), np.asarray(self.xv, dtype=float))

    def __iter__(self) -> Iterator[object]:
        yield float(self.value)
        yield np.asarray(self.xl, dtype=float)
        yield np.asarray(self.xv, dtype=float)


@dataclass(frozen=True, slots=True, eq=False)
class VaporizationResult:
    """Structured result for enthalpy-of-vaporization calculations."""

    hvap: float
    pressure: float

    def to_legacy(self):
        return [float(self.hvap), float(self.pressure)]

    def __iter__(self) -> Iterator[object]:
        yield float(self.hvap)
        yield float(self.pressure)


@dataclass(frozen=True, slots=True, eq=False)
class PhaseResult:
    """One phase in a multiphase LLE result."""

    beta: float
    x: np.ndarray
    rho: float
    lnfugcoef: np.ndarray
    lnfug: np.ndarray

    def to_legacy(self) -> dict:
        return {
            "beta": float(self.beta),
            "x": np.asarray(self.x, dtype=float),
            "rho": float(self.rho),
            "lnfugcoef": np.asarray(self.lnfugcoef, dtype=float),
            "lnfug": np.asarray(self.lnfug, dtype=float),
        }


@dataclass(frozen=True, slots=True, eq=False)
class MultiphaseLLEResult:
    """Structured result for the multiphase electrolyte LLE solver."""

    n_phases: int
    phases: tuple[PhaseResult, ...]
    tpdf_min: float
    tpdf_seed_x: np.ndarray
    converged: bool
    status: int
    message: str
    residual_norm: float
    e_matrix: np.ndarray
    ion_pair_rows: tuple[dict, ...]
    charged_species: tuple[str, ...]
    charged_species_indices: np.ndarray

    def as_dict(self) -> dict:
        return {
            "n_phases": int(self.n_phases),
            "phases": [phase.to_legacy() for phase in self.phases],
            "tpdf_min": float(self.tpdf_min),
            "tpdf_seed_x": np.asarray(self.tpdf_seed_x, dtype=float),
            "converged": bool(self.converged),
            "status": int(self.status),
            "message": str(self.message),
            "residual_norm": float(self.residual_norm),
            "e_matrix": np.asarray(self.e_matrix, dtype=float),
            "ion_pair_rows": list(self.ion_pair_rows),
            "charged_species": list(self.charged_species),
            "charged_species_indices": np.asarray(self.charged_species_indices, dtype=int),
        }

    def to_legacy(self) -> dict:
        return self.as_dict()


@dataclass(frozen=True, slots=True, eq=False)
class PCSAFTMixture:
    """Model definition for a PC-SAFT/ePC-SAFT system."""

    species: tuple[str, ...]
    source: str
    dataset_name: str | None = None
    user_options: dict = field(default_factory=dict, repr=False, compare=False)
    params: dict | None = field(default=None, repr=False, compare=False)

    @classmethod
    def from_dataset(
        cls,
        dataset: str,
        species: Iterable[str],
        *,
        user_options: dict | None = None,
    ) -> "PCSAFTMixture":
        return cls(
            species=_as_species_tuple(species),
            source="dataset",
            dataset_name=str(dataset),
            user_options=_copy_params(user_options or {}),
        )

    @classmethod
    def from_params(
        cls,
        species: Iterable[str],
        params: dict,
    ) -> "PCSAFTMixture":
        return cls(
            species=_as_species_tuple(species),
            source="params",
            params=_copy_params(params),
        )

    @classmethod
    def from_raw(
        cls,
        species: Iterable[str],
        params: dict,
    ) -> "PCSAFTMixture":
        return cls.from_params(species, params)

    def _resolve_params(self, *, T: float, x) -> dict:
        if self.source == "dataset":
            if self.dataset_name is None:
                raise InputError("dataset_name is required for dataset-backed mixtures.")
            return get_prop_dict(
                self.dataset_name,
                self.species,
                np.asarray(x, dtype=float),
                float(T),
                user_options=_copy_params(self.user_options),
            )
        if self.params is None:
            raise InputError("params are required for raw mixtures.")
        return _copy_params(self.params)

    def state(
        self,
        *,
        T: float,
        x,
        P: float | None = None,
        rho: float | None = None,
        phase: str = "liq",
    ) -> "PCSAFTState":
        return PCSAFTState.from_mixture(self, T=T, x=x, P=P, rho=rho, phase=phase)

    @property
    def n_components(self) -> int:
        return len(self.species)


@dataclass(frozen=True, slots=True, eq=False)
class PCSAFTState:
    """Bound thermodynamic state for a PC-SAFT/ePC-SAFT mixture."""

    mixture: PCSAFTMixture
    T: float
    x: np.ndarray
    phase: str = "liq"
    P: float | None = None
    rho: float | None = None
    _params: dict = field(repr=False, compare=False, default_factory=dict)
    _pressure_cache: float | None = field(repr=False, compare=False, default=None)
    _rho_cache: float | None = field(repr=False, compare=False, default=None)

    def __post_init__(self) -> None:
        if not isinstance(self.mixture, PCSAFTMixture):
            raise InputError("mixture must be a PCSAFTMixture instance.")
        T, x_arr, P, rho = _validate_state_inputs(self.T, self.x, self.P, self.rho)
        if len(x_arr) != self.mixture.n_components:
            raise InputError(
                f"x length ({len(x_arr)}) must match mixture size ({self.mixture.n_components})."
            )
        resolved = self._params or self.mixture._resolve_params(T=T, x=x_arr)
        object.__setattr__(self, "T", T)
        object.__setattr__(self, "x", x_arr)
        object.__setattr__(self, "P", P)
        object.__setattr__(self, "rho", rho)
        object.__setattr__(self, "_params", _copy_params(resolved))
        object.__setattr__(self, "_pressure_cache", P)
        object.__setattr__(self, "_rho_cache", rho)

    @classmethod
    def from_mixture(
        cls,
        mixture: PCSAFTMixture,
        *,
        T: float,
        x,
        P: float | None = None,
        rho: float | None = None,
        phase: str = "liq",
    ) -> "PCSAFTState":
        return cls(mixture=mixture, T=T, x=_as_float_array(x, name="x"), phase=str(phase), P=P, rho=rho)

    def with_conditions(
        self,
        *,
        T: float | None = None,
        x=None,
        P: float | None = None,
        rho: float | None = None,
        phase: str | None = None,
    ) -> "PCSAFTState":
        return self.mixture.state(
            T=self.T if T is None else T,
            x=self.x if x is None else x,
            P=self.P if P is None else P,
            rho=self.rho if rho is None else rho,
            phase=self.phase if phase is None else phase,
        )

    def _params_for_call(self) -> dict:
        return _copy_params(self._params)

    def _density_for_call(self) -> float:
        rho = self._rho_cache
        if rho is not None:
            return float(rho)
        pressure = self._pressure_cache if self._pressure_cache is not None else self.pressure()
        rho = float(pcsaft_den(self.T, pressure, self.x, self._params_for_call(), phase=self.phase))
        object.__setattr__(self, "_rho_cache", rho)
        return rho

    def _pressure_for_call(self) -> float:
        pressure = self._pressure_cache
        if pressure is not None:
            return float(pressure)
        rho = self._rho_cache if self._rho_cache is not None else self.density()
        pressure = float(pcsaft_p(self.T, rho, self.x, self._params_for_call()))
        object.__setattr__(self, "_pressure_cache", pressure)
        return pressure

    def density(self) -> float:
        return self._density_for_call()

    def pressure(self) -> float:
        return self._pressure_for_call()

    def Z(self) -> float:
        rho = self.density()
        return float(pcsaft_Z(self.T, rho, self.x, self._params_for_call()))

    def ares(self) -> float:
        rho = self.density()
        return float(pcsaft_ares(self.T, rho, self.x, self._params_for_call()))

    def dadt(self) -> float:
        rho = self.density()
        return float(pcsaft_dadt(self.T, rho, self.x, self._params_for_call()))

    def hres(self) -> float:
        rho = self.density()
        return float(pcsaft_hres(self.T, rho, self.x, self._params_for_call()))

    def sres(self) -> float:
        rho = self.density()
        return float(pcsaft_sres(self.T, rho, self.x, self._params_for_call()))

    def gres(self) -> float:
        rho = self.density()
        return float(pcsaft_gres(self.T, rho, self.x, self._params_for_call()))

    def lnfugcoef(self) -> np.ndarray:
        rho = self.density()
        return np.asarray(pcsaft_lnfugcoef(self.T, rho, self.x, self._params_for_call()), dtype=float)

    def fugcoef(self) -> np.ndarray:
        rho = self.density()
        return np.asarray(pcsaft_fugcoef(self.T, rho, self.x, self._params_for_call()), dtype=float)

    def mures(self) -> np.ndarray:
        return R_GAS * self.T * self.lnfugcoef()

    def lnfugcoef_terms(self) -> dict:
        rho = self.density()
        return pcsaft_lnfugcoef_terms(self.T, rho, self.x, self._params_for_call())

    def dielc_eval(self):
        return pcsaft_dielc_eval(self.x, self._params_for_call())

    def osmoticC(self) -> float:
        rho = self.density()
        return np.asarray(pcsaft_osmoticC(self.T, rho, self.x, self._params_for_call()), dtype=float)

    def miac_m(self, species: Iterable[str] | None = None) -> dict:
        rho = self.density()
        use_species = tuple(species) if species is not None else self.mixture.species
        return pcsaft_miac_m(self.T, rho, self.x, self._params_for_call(), species=list(use_species))

    def miac(self, species: Iterable[str] | None = None) -> dict:
        rho = self.density()
        use_species = tuple(species) if species is not None else self.mixture.species
        return pcsaft_miac(self.T, rho, self.x, self._params_for_call(), species=list(use_species))

    def gsolv(self, species: Iterable[str] | None = None):
        rho = self.density()
        use_species = tuple(species) if species is not None else self.mixture.species
        return pcsaft_gsolv(self.T, rho, self.x, self._params_for_call(), species=list(use_species))

    def cp(self, aly_lee_params) -> float:
        rho = self.density()
        return float(pcsaft_cp(self.T, rho, aly_lee_params, self.x, self._params_for_call()))

    def flashTQ(self, q: float, p_guess: float | None = None) -> FlashResult:
        result = flashTQ(self.T, q, self.x, self._params_for_call(), p_guess=p_guess)
        value, xl, xv = result
        return FlashResult(float(value), _result_array(xl), _result_array(xv), kind="pressure")

    def flashPQ(self, q: float, t_guess: float | None = None) -> FlashResult:
        pressure = self.pressure()
        result = flashPQ(pressure, q, self.x, self._params_for_call(), t_guess=t_guess)
        value, xl, xv = result
        return FlashResult(float(value), _result_array(xl), _result_array(xv), kind="temperature")

    def Hvap(self, p_guess: float | None = None) -> VaporizationResult:
        hvap = pcsaft_Hvap(self.T, self.x, self._params_for_call(), p_guess=p_guess)
        return VaporizationResult(float(hvap[0]), float(hvap[1]))

    def multiphase_lle(self, z_feed=None, options: dict | None = None) -> MultiphaseLLEResult:
        feed = self.x if z_feed is None else _as_float_array(z_feed, name="z_feed")
        result = pcsaft_multiphase_lle(
            self.T,
            self.pressure(),
            feed,
            self._params_for_call(),
            list(self.mixture.species),
            options=options or {},
        )
        phases = tuple(
            PhaseResult(
                beta=float(phase["beta"]),
                x=_result_array(phase["x"]),
                rho=float(phase["rho"]),
                lnfugcoef=_result_array(phase["lnfugcoef"]),
                lnfug=_result_array(phase["lnfug"]),
            )
            for phase in result["phases"]
        )
        return MultiphaseLLEResult(
            n_phases=int(result["n_phases"]),
            phases=phases,
            tpdf_min=float(result["tpdf_min"]),
            tpdf_seed_x=_result_array(result["tpdf_seed_x"]),
            converged=bool(result["converged"]),
            status=int(result["status"]),
            message=str(result["message"]),
            residual_norm=float(result["residual_norm"]),
            e_matrix=_result_array(result["e_matrix"]),
            ion_pair_rows=tuple(result["ion_pair_rows"]),
            charged_species=tuple(result["charged_species"]),
            charged_species_indices=np.asarray(result["charged_species_indices"], dtype=int),
        )


def mixture_from_legacy_params(params: dict, species: Iterable[str] | None = None) -> PCSAFTMixture:
    """Build a raw-parameter mixture from the legacy flat API payload."""

    if species is None:
        species = _infer_species_from_params(params)
    return PCSAFTMixture.from_params(species, params)
