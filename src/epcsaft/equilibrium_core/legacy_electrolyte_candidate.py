"""Legacy-style electrolyte LLE candidate generation for diagnostics only."""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from .electrolyte_basis import build_electrolyte_basis


def legacy_electrolyte_lle_candidate(
    mixture: Any,
    *,
    T: float,
    P: float,
    feed: Any,
    options: Any,
) -> dict[str, Any]:
    """Return a best-effort two-phase electrolyte candidate without accepting it."""
    from scipy.optimize import least_squares

    z_feed = np.asarray(feed, dtype=float).flatten()
    charges = np.asarray(mixture.parameters["z"], dtype=float).flatten()
    species = list(getattr(mixture, "species", [str(i) for i in range(z_feed.size)]))
    basis = build_electrolyte_basis(species, charges, z_feed)
    beta_bounds = (1.0e-8, 1.0 - 1.0e-8)
    starts = _candidate_starts(z_feed, charges, basis.to_dict(), float(options.min_composition))
    beta_starts = (0.02, 0.05, 0.10, 0.25, 0.50)
    candidates: list[dict[str, Any]] = []

    for start in starts:
        for beta0 in beta_starts:
            y0 = _pack_variables(start, beta0, beta_bounds)
            try:
                sol = least_squares(
                    _residual,
                    y0,
                    args=(
                        mixture,
                        float(T),
                        float(P),
                        z_feed,
                        charges,
                        basis.to_dict(),
                        beta_bounds,
                        float(options.min_composition),
                        False,
                    ),
                    method="trf",
                    ftol=max(float(options.tolerance), 1.0e-10),
                    xtol=max(float(options.tolerance), 1.0e-10),
                    gtol=max(float(options.tolerance), 1.0e-10),
                    max_nfev=int(options.legacy_candidate_max_iterations),
                )
            except Exception:
                continue
            candidate = _candidate_from_variables(
                mixture,
                float(T),
                float(P),
                z_feed,
                charges,
                basis.to_dict(),
                sol.x,
                beta_bounds,
                float(options.min_composition),
            )
            if candidate is not None:
                candidate["status"] = int(sol.status)
                candidate["message"] = str(sol.message)
                candidates.append(candidate)

    best = _choose_best_candidate(candidates, float(options.legacy_candidate_split_tolerance))
    if best is None:
        best = _heuristic_candidate(z_feed, charges, starts, float(options.min_composition))
        if best is None:
            return {
                "legacy_candidate_found": False,
                "legacy_candidate_message": "legacy candidate fallback did not produce a feasible two-phase candidate",
            }

    if best["phase_distance"] <= float(options.legacy_candidate_split_tolerance):
        anchored = _anchored_retry(
            mixture,
            float(T),
            float(P),
            z_feed,
            charges,
            basis.to_dict(),
            starts,
            beta_bounds,
            float(options.min_composition),
            options,
        )
        if anchored is not None:
            best = anchored

    loose_ok = (
        best["residual_norm"] <= float(options.legacy_candidate_residual_tolerance)
        and best["material_balance_error"] <= 1.0e-8
        and best["charge_balance_error"] <= 1.0e-6
        and best["phase_distance"] > float(options.legacy_candidate_split_tolerance)
    )
    strict_ok = (
        best["residual_norm"] <= float(options.tolerance)
        and best["material_balance_error"] <= max(float(options.tolerance), 1.0e-10)
        and best["charge_balance_error"] <= 1.0e-8
        and best["gibbs_delta"] < 0.0
        and best["phase_distance"] > max(1.0e-4, float(options.legacy_candidate_split_tolerance))
    )
    return {
        "legacy_candidate_found": True,
        "legacy_candidate_accepted_by_loose_gate": bool(loose_ok),
        "legacy_candidate_rejected_by_strict_gate": bool(not strict_ok),
        "legacy_candidate_phase_distance": float(best["phase_distance"]),
        "legacy_candidate_residual_norm": float(best["residual_norm"]),
        "legacy_candidate_material_balance_error": float(best["material_balance_error"]),
        "legacy_candidate_charge_balance_error": float(best["charge_balance_error"]),
        "legacy_candidate_gibbs_delta": float(best["gibbs_delta"]),
        "legacy_candidate_aq_composition": np.asarray(best["aq"], dtype=float).tolist(),
        "legacy_candidate_org_composition": np.asarray(best["org"], dtype=float).tolist(),
        "legacy_candidate_phase_fraction": float(best["beta_org"]),
        "legacy_candidate_split_tolerance": float(options.legacy_candidate_split_tolerance),
        "legacy_candidate_residual_tolerance": float(options.legacy_candidate_residual_tolerance),
        "legacy_candidate_status": int(best.get("status", 0)),
        "legacy_candidate_message": str(best.get("message", "")),
    }


def _safe_log(x: np.ndarray) -> np.ndarray:
    return np.log(np.clip(np.asarray(x, dtype=float), 1.0e-300, None))


def _softmax(u: np.ndarray) -> np.ndarray:
    shifted = np.asarray(u, dtype=float) - float(np.max(u))
    exp_u = np.exp(np.clip(shifted, -700.0, 700.0))
    return exp_u / float(np.sum(exp_u))


def _logit(value: float) -> float:
    v = min(max(float(value), 1.0e-12), 1.0 - 1.0e-12)
    return math.log(v / (1.0 - v))


def _sigmoid(value: float) -> float:
    if value >= 0.0:
        exp_v = math.exp(-value)
        return 1.0 / (1.0 + exp_v)
    exp_v = math.exp(value)
    return exp_v / (1.0 + exp_v)


def _pack_variables(org: np.ndarray, beta: float, beta_bounds: tuple[float, float]) -> np.ndarray:
    org = np.clip(np.asarray(org, dtype=float), 1.0e-300, None)
    org = org / float(np.sum(org))
    lo, hi = beta_bounds
    frac = (float(beta) - lo) / (hi - lo)
    return np.asarray([*np.log(org[:-1] / org[-1]), _logit(frac)], dtype=float)


def _unpack_variables(y: np.ndarray, beta_bounds: tuple[float, float]) -> tuple[np.ndarray, float]:
    logits = np.concatenate([np.asarray(y[:-1], dtype=float), np.zeros(1, dtype=float)])
    org = _softmax(logits)
    lo, hi = beta_bounds
    beta = lo + (hi - lo) * _sigmoid(float(y[-1]))
    return org, float(beta)


def _phase_state(mixture: Any, T: float, P: float, composition: np.ndarray) -> dict[str, Any]:
    state = mixture.state(T=T, P=P, x=composition, phase="liq")
    ln_phi = np.asarray(state.fugacity_coefficient(), dtype=float)
    return {"ln_phi": ln_phi, "g": float(np.sum(composition * (_safe_log(composition) + ln_phi)))}


def _equilibrium_residual(
    aq: np.ndarray,
    org: np.ndarray,
    aq_state: dict[str, Any],
    org_state: dict[str, Any],
    basis: dict[str, Any],
    charges: np.ndarray,
) -> np.ndarray:
    aq_lnf = _safe_log(aq) + np.asarray(aq_state["ln_phi"], dtype=float)
    org_lnf = _safe_log(org) + np.asarray(org_state["ln_phi"], dtype=float)
    residuals = [float(org_lnf[int(index)] - aq_lnf[int(index)]) for index in basis["neutral_indices"]]
    for pair in basis["salt_pairs"]:
        c = int(pair["cation"])
        a = int(pair["anion"])
        residuals.append(
            float(pair.get("cation_stoich", 1)) * float(org_lnf[c] - aq_lnf[c])
            + float(pair.get("anion_stoich", 1)) * float(org_lnf[a] - aq_lnf[a])
        )
    residuals.append(float(np.dot(charges, org)))
    return np.asarray(residuals, dtype=float)


def _residual(
    y: np.ndarray,
    mixture: Any,
    T: float,
    P: float,
    feed: np.ndarray,
    charges: np.ndarray,
    basis: dict[str, Any],
    beta_bounds: tuple[float, float],
    min_composition: float,
    anchored: bool,
) -> np.ndarray:
    penalty_size = len(basis["neutral_indices"]) + len(basis["salt_pairs"]) + 1 + int(anchored)
    try:
        org, beta = _unpack_variables(y, beta_bounds)
        aq = (feed - beta * org) / (1.0 - beta)
        if np.any(~np.isfinite(aq)) or np.any(aq <= min_composition):
            return np.full(penalty_size, 1.0e6, dtype=float)
        aq = aq / float(np.sum(aq))
        aq_state = _phase_state(mixture, T, P, aq)
        org_state = _phase_state(mixture, T, P, org)
        residual = _equilibrium_residual(aq, org, aq_state, org_state, basis, charges)
        if anchored:
            organic_indices = [int(i) for i in basis["neutral_indices"][1:]] or [int(basis["neutral_indices"][0])]
            anchor_value = sum(org[i] for i in organic_indices)
            residual = np.concatenate([residual, np.asarray([anchor_value - 0.75], dtype=float)])
        return residual
    except Exception:
        return np.full(penalty_size, 1.0e6, dtype=float)


def _candidate_from_variables(
    mixture: Any,
    T: float,
    P: float,
    feed: np.ndarray,
    charges: np.ndarray,
    basis: dict[str, Any],
    variables: np.ndarray,
    beta_bounds: tuple[float, float],
    min_composition: float,
) -> dict[str, Any] | None:
    try:
        org, beta = _unpack_variables(variables, beta_bounds)
        aq = (feed - beta * org) / (1.0 - beta)
        if np.any(~np.isfinite(aq)) or np.any(aq <= min_composition):
            return None
        aq = aq / float(np.sum(aq))
        aq_state = _phase_state(mixture, T, P, aq)
        org_state = _phase_state(mixture, T, P, org)
        residual = _equilibrium_residual(aq, org, aq_state, org_state, basis, charges)
    except Exception:
        return None
    material = (1.0 - beta) * aq + beta * org - feed
    charge_error = max(
        abs(float(np.dot(charges, feed))), abs(float(np.dot(charges, aq))), abs(float(np.dot(charges, org)))
    )
    gibbs_feed = _phase_state(mixture, T, P, feed)["g"]
    gibbs_split = (1.0 - beta) * float(aq_state["g"]) + beta * float(org_state["g"])
    return {
        "aq": aq,
        "org": org,
        "beta_org": float(beta),
        "residual_norm": float(np.max(np.abs(residual))) if residual.size else 0.0,
        "material_balance_error": float(np.max(np.abs(material))),
        "charge_balance_error": float(charge_error),
        "gibbs_delta": float(gibbs_split - gibbs_feed),
        "phase_distance": float(np.max(np.abs(aq - org))),
    }


def _choose_best_candidate(candidates: list[dict[str, Any]], split_tolerance: float) -> dict[str, Any] | None:
    if not candidates:
        return None
    split = [item for item in candidates if item["phase_distance"] > split_tolerance]
    pool = split if split else candidates
    return min(pool, key=lambda item: (item["residual_norm"], -item["phase_distance"]))


def _heuristic_candidate(
    feed: np.ndarray, charges: np.ndarray, starts: list[np.ndarray], min_composition: float
) -> dict[str, Any] | None:
    candidates: list[dict[str, Any]] = []
    for org_seed in starts:
        org = np.clip(np.asarray(org_seed, dtype=float), min_composition, None)
        org = org / float(np.sum(org))
        beta = _max_material_balanced_beta(feed, org, min_composition)
        if beta <= 0.0:
            continue
        aq = (feed - beta * org) / (1.0 - beta)
        if np.any(~np.isfinite(aq)) or np.any(aq <= min_composition):
            continue
        aq = aq / float(np.sum(aq))
        material = (1.0 - beta) * aq + beta * org - feed
        candidates.append(
            {
                "aq": aq,
                "org": org,
                "beta_org": float(beta),
                "residual_norm": 1.0e300,
                "material_balance_error": float(np.max(np.abs(material))),
                "charge_balance_error": max(
                    abs(float(np.dot(charges, feed))),
                    abs(float(np.dot(charges, aq))),
                    abs(float(np.dot(charges, org))),
                ),
                "gibbs_delta": 0.0,
                "phase_distance": float(np.max(np.abs(aq - org))),
                "status": -1,
                "message": "thermodynamic fallback failed; preserving material-balanced heuristic split",
            }
        )
    if not candidates:
        return None
    return max(candidates, key=lambda item: item["phase_distance"])


def _max_material_balanced_beta(feed: np.ndarray, org: np.ndarray, min_composition: float) -> float:
    limits: list[float] = []
    for feed_i, org_i in zip(np.asarray(feed, dtype=float), np.asarray(org, dtype=float), strict=True):
        if org_i > feed_i:
            limits.append((float(feed_i) - min_composition) / max(float(org_i) - min_composition, min_composition))
    upper = min(limits) if limits else 0.5
    return float(max(0.0, min(0.25, 0.8 * upper)))


def _candidate_starts(
    feed: np.ndarray, charges: np.ndarray, basis: dict[str, Any], min_composition: float
) -> list[np.ndarray]:
    starts = [np.asarray(feed, dtype=float)]
    neutral_indices = [int(i) for i in basis["neutral_indices"]]
    charged_indices = [i for i, charge in enumerate(charges) if abs(float(charge)) > 1.0e-12]
    organic_indices = neutral_indices[1:] if len(neutral_indices) > 1 else neutral_indices
    salt_total = max(float(np.sum(feed[charged_indices])) if charged_indices else 0.0, 2.0 * min_composition)
    for water_share, organic_share, salt_share in ((0.20, 0.75, 0.05), (0.05, 0.90, 0.05), (0.65, 0.30, 0.05)):
        candidate = np.full(feed.size, min_composition, dtype=float)
        if neutral_indices:
            candidate[neutral_indices[0]] = water_share
        for index in organic_indices:
            candidate[index] = organic_share / max(len(organic_indices), 1)
        _copy_scaled_ions(candidate, feed, charges, max(salt_share, salt_total), min_composition)
        starts.append(candidate / float(np.sum(candidate)))
    for index in neutral_indices:
        candidate = np.full(feed.size, min_composition, dtype=float)
        candidate[index] = 0.90
        for other in neutral_indices:
            if other != index:
                candidate[other] = 0.05 / max(len(neutral_indices) - 1, 1)
        _copy_scaled_ions(candidate, feed, charges, salt_total, min_composition)
        starts.append(candidate / float(np.sum(candidate)))
    return starts


def _copy_scaled_ions(
    candidate: np.ndarray, feed: np.ndarray, charges: np.ndarray, ion_total: float, min_composition: float
) -> None:
    charged_indices = [i for i, charge in enumerate(charges) if abs(float(charge)) > 1.0e-12]
    if not charged_indices:
        return
    feed_ions = np.asarray([feed[i] for i in charged_indices], dtype=float)
    total = float(np.sum(feed_ions))
    if total <= 0.0:
        return
    for local, index in enumerate(charged_indices):
        candidate[index] = max(min_composition, ion_total * float(feed_ions[local]) / total)


def _anchored_retry(
    mixture: Any,
    T: float,
    P: float,
    feed: np.ndarray,
    charges: np.ndarray,
    basis: dict[str, Any],
    starts: list[np.ndarray],
    beta_bounds: tuple[float, float],
    min_composition: float,
    options: Any,
) -> dict[str, Any] | None:
    from scipy.optimize import least_squares

    candidates: list[dict[str, Any]] = []
    for start in starts:
        for beta0 in (0.02, 0.05, 0.10, 0.25):
            y0 = _pack_variables(start, beta0, beta_bounds)
            try:
                sol = least_squares(
                    _residual,
                    y0,
                    args=(mixture, T, P, feed, charges, basis, beta_bounds, min_composition, True),
                    method="trf",
                    ftol=max(float(options.tolerance), 1.0e-10),
                    xtol=max(float(options.tolerance), 1.0e-10),
                    gtol=max(float(options.tolerance), 1.0e-10),
                    max_nfev=int(options.legacy_candidate_max_iterations),
                )
            except Exception:
                continue
            candidate = _candidate_from_variables(
                mixture, T, P, feed, charges, basis, sol.x, beta_bounds, min_composition
            )
            if candidate is not None:
                candidate["status"] = int(sol.status)
                candidate["message"] = str(sol.message)
                candidates.append(candidate)
    split = [item for item in candidates if item["phase_distance"] > float(options.legacy_candidate_split_tolerance)]
    if not split:
        return None
    return min(split, key=lambda item: (item["residual_norm"], -item["phase_distance"]))
