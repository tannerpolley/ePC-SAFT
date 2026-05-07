"""Low-overhead repeated property helpers."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from ._types import InputError
from .epcsaft import ePCSAFTMixture


def evaluate_fugacity_coefficients(
    mixture: ePCSAFTMixture,
    *,
    T: float,
    x: Any,
    P: float | None = None,
    rho: float | None = None,
    phase: str = "liq",
    rho_seed: float | None = None,
    rho_guess: float | None = None,
    natural_log: bool = True,
) -> dict[str, Any]:
    """Evaluate fugacity coefficients plus resolved density for one state."""

    if not isinstance(mixture, ePCSAFTMixture):
        raise InputError("mixture must be an ePCSAFTMixture instance.")
    if rho_seed is not None and rho_guess is not None:
        raise InputError("Provide only one of rho_seed or rho_guess.")
    seed = rho_guess if rho_guess is not None else rho_seed
    state = mixture.state(T=T, x=x, P=P, rho=rho, phase=phase, rho_guess=seed)
    ln_phi = np.asarray(state.fugacity_coefficient(natural_log=True), dtype=float)
    phi = np.exp(ln_phi)
    return {
        "ln_fugacity_coefficient": ln_phi,
        "fugacity_coefficient": phi,
        "selected_fugacity_output": ln_phi if natural_log else phi,
        "natural_log": bool(natural_log),
        "density": float(state.molar_density()),
        "pressure": float(state.pressure()),
        "temperature": float(state.T),
        "phase": _phase_label(phase),
    }


def evaluate_fugacity_coefficients_batch(
    mixture: ePCSAFTMixture,
    *,
    rows: Sequence[Mapping[str, Any]],
    natural_log: bool = True,
    continuation: str = "auto",
) -> list[dict[str, Any]]:
    """Evaluate fugacity coefficients for many rows, reusing nearby density seeds."""

    mode = str(continuation).strip().lower()
    if mode not in {"auto", "none"}:
        raise InputError("continuation must be 'auto' or 'none'.")
    out: list[dict[str, Any]] = []
    last_density: float | None = None
    for row in rows:
        seed = row.get("rho_seed", row.get("rho_guess"))
        if seed is None and mode == "auto" and "P" in row and last_density is not None:
            seed = last_density
        payload = evaluate_fugacity_coefficients(
            mixture,
            T=float(row["T"]),
            x=row["x"],
            P=row.get("P"),
            rho=row.get("rho"),
            phase=str(row.get("phase", "liq")),
            rho_seed=seed,
            natural_log=natural_log,
        )
        out.append(payload)
        last_density = float(payload["density"])
    return out


def _phase_label(value: str) -> str:
    token = str(value).strip().lower()
    if token in {"liq", "liquid", "aq", "org"}:
        return "liq"
    if token in {"vap", "vapor", "gas"}:
        return "vap"
    return token
