from __future__ import annotations

import argparse
import importlib
import json
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts._env import require_pcsaft_install

require_pcsaft_install()

try:
    from cyipopt import minimize_ipopt
except Exception as exc:  # pragma: no cover
    minimize_ipopt = None
    _CYIPOPT_IMPORT_ERROR = exc
else:
    _CYIPOPT_IMPORT_ERROR = None

pcs_core = importlib.import_module("pcsaft.pcsaft")


def _require_cyipopt() -> None:
    if minimize_ipopt is None:
        raise ImportError(
            "cyipopt is required for solve_two_phase_lle_cyipopt(). "
            "Install it with `conda install -n PC-SAFT -c conda-forge ipopt cyipopt`."
        ) from _CYIPOPT_IMPORT_ERROR


def _softmax(logits: np.ndarray) -> np.ndarray:
    shifted = np.asarray(logits, dtype=float) - np.max(logits)
    expo = np.exp(shifted)
    return expo / np.sum(expo)


def _sigmoid(value: float) -> float:
    if value >= 0.0:
        exp_neg = np.exp(-value)
        return 1.0 / (1.0 + exp_neg)
    exp_pos = np.exp(value)
    return exp_pos / (1.0 + exp_pos)


def _logit(value: float) -> float:
    clipped = min(max(float(value), 1.0e-15), 1.0 - 1.0e-15)
    return float(np.log(clipped / (1.0 - clipped)))


def _pack_start(x_start: np.ndarray, beta0: float, beta_bounds: tuple[float, float]) -> np.ndarray:
    x_start = np.asarray(x_start, dtype=float)
    y0 = np.zeros(x_start.size, dtype=float)
    y0[:-1] = np.log(np.maximum(x_start[:-1], 1.0e-15) / max(float(x_start[-1]), 1.0e-15))
    beta_lo, beta_hi = float(beta_bounds[0]), float(beta_bounds[1])
    frac = (float(beta0) - beta_lo) / (beta_hi - beta_lo)
    y0[-1] = _logit(frac)
    return y0


def _unpack_y(y: np.ndarray, z_feed: np.ndarray, beta_bounds: tuple[float, float], x_floor: float = 1.0e-15) -> dict:
    y = np.asarray(y, dtype=float).flatten()
    ncomp = int(len(z_feed))
    if y.size != ncomp:
        return {"valid": False, "x1": np.full(ncomp, np.nan), "x2": np.full(ncomp, np.nan), "beta": np.nan}
    x1 = _softmax(np.concatenate([y[:-1], np.zeros(1, dtype=float)]))
    beta_lo, beta_hi = float(beta_bounds[0]), float(beta_bounds[1])
    beta = beta_lo + (beta_hi - beta_lo) * _sigmoid(float(y[-1]))
    denom = 1.0 - beta
    if denom <= 0.0:
        return {"valid": False, "x1": x1, "x2": np.full(ncomp, np.nan), "beta": float(beta)}
    x2 = (np.asarray(z_feed, dtype=float) - beta * x1) / denom
    valid = bool(
        np.all(np.isfinite(x2))
        and np.all(x2 > x_floor)
        and abs(float(np.sum(x2)) - 1.0) <= 1.0e-8
        and 0.0 < beta < 1.0
    )
    return {"valid": valid, "x1": x1, "x2": x2, "beta": float(beta)}


def _build_starts(seed_x: np.ndarray, z_feed: np.ndarray, random_seed: int) -> list[np.ndarray]:
    rng = np.random.default_rng(int(random_seed))
    starts = [
        np.asarray(seed_x, dtype=float),
        np.asarray(z_feed, dtype=float),
        0.8 * np.asarray(seed_x, dtype=float) + 0.2 * np.asarray(z_feed, dtype=float),
        0.2 * np.asarray(seed_x, dtype=float) + 0.8 * np.asarray(z_feed, dtype=float),
    ]
    for _ in range(8):
        blend = float(rng.uniform(0.1, 0.9))
        jitter = rng.random(len(z_feed)) + 1.0e-10
        starts.append(blend * np.asarray(seed_x, dtype=float) + (1.0 - blend) * jitter / np.sum(jitter))
    return [np.maximum(start, 1.0e-12) / np.sum(np.maximum(start, 1.0e-12)) for start in starts]


def _residual_vector(y: np.ndarray, ctx: dict, anchor: dict | None = None) -> np.ndarray:
    unpacked = _unpack_y(y, ctx["z_feed"], ctx["beta_bounds"], x_floor=ctx["x_floor"])
    nres = ctx["ncomp"] + (1 if anchor is not None else 0)
    if not unpacked["valid"]:
        out = np.full(nres, ctx["invalid_scale"], dtype=float)
        x2 = np.asarray(unpacked["x2"], dtype=float)
        if np.any(np.isfinite(x2)):
            out *= 1.0 + float(np.sum(np.maximum(ctx["x_floor"] - x2, 0.0)))
        return out
    x1 = unpacked["x1"]
    x2 = unpacked["x2"]
    try:
        st1 = pcs_core._phase_state_liq(ctx["t"], ctx["p"], x1, ctx["params"])
        st2 = pcs_core._phase_state_liq(ctx["t"], ctx["p"], x2, ctx["params"])
    except Exception:
        return np.full(nres, 2.0 * ctx["invalid_scale"], dtype=float)
    res_neutral = np.asarray(st1["lnfug"], dtype=float)[ctx["neutral_idx"]] - np.asarray(st2["lnfug"], dtype=float)[ctx["neutral_idx"]]
    delta_ch = np.asarray(st1["lnfug"], dtype=float)[ctx["charged_idx"]] - np.asarray(st2["lnfug"], dtype=float)[ctx["charged_idx"]]
    res_ionic = ctx["E"].dot(delta_ch)
    res_charge = np.array([float(ctx["charge_weight"] * np.dot(ctx["z"], x1))], dtype=float)
    parts = [res_neutral, res_ionic, res_charge]
    if anchor is not None:
        parts.append(np.array([float(anchor.get("weight", 1.0)) * (x1[int(anchor["component"])] - float(anchor["target"]))], dtype=float))
    return np.concatenate(parts)


def _objective(y: np.ndarray, ctx: dict, anchor: dict | None = None) -> float:
    residual = _residual_vector(y, ctx, anchor=anchor)
    unpacked = _unpack_y(y, ctx["z_feed"], ctx["beta_bounds"], x_floor=ctx["x_floor"])
    smooth_penalty = 0.0
    if np.any(np.isfinite(unpacked["x2"])):
        smooth_penalty += float(np.sum(np.maximum(ctx["x_floor"] - unpacked["x2"], 0.0) ** 2))
    if not unpacked["valid"]:
        smooth_penalty += 1.0
    return float(0.5 * np.dot(residual, residual) + ctx["invalid_objective"] * smooth_penalty)


def _x2_positivity(y: np.ndarray, ctx: dict) -> np.ndarray:
    unpacked = _unpack_y(y, ctx["z_feed"], ctx["beta_bounds"], x_floor=ctx["x_floor"])
    x2 = np.asarray(unpacked["x2"], dtype=float)
    if np.any(~np.isfinite(x2)):
        return np.full(ctx["ncomp"], -1.0, dtype=float)
    return x2 - ctx["x_floor"]


def _charge_balance(y: np.ndarray, ctx: dict) -> np.ndarray:
    unpacked = _unpack_y(y, ctx["z_feed"], ctx["beta_bounds"], x_floor=ctx["x_floor"])
    if not unpacked["valid"]:
        return np.array([1.0], dtype=float)
    return np.array([float(np.dot(ctx["z"], unpacked["x1"]))], dtype=float)


def _ipopt_options(options: dict) -> dict:
    return {
        "max_iter": int(options.get("max_nfev", 200)),
        "print_level": int(options.get("ipopt_print_level", 0)),
        "mu_strategy": options.get("ipopt_mu_strategy", "adaptive"),
        "hessian_approximation": options.get("ipopt_hessian_approximation", "limited-memory"),
        "nlp_scaling_method": options.get("ipopt_scaling", "gradient-based"),
        "eps": float(options.get("fd_eps", 1.0e-7)),
    }


def _solve_start(y0: np.ndarray, ctx: dict, options: dict, anchor: dict | None = None):
    constraints = (
        {"type": "ineq", "fun": _x2_positivity, "args": (ctx,)},
        {"type": "eq", "fun": _charge_balance, "args": (ctx,)},
    )
    return minimize_ipopt(
        _objective,
        y0,
        args=(ctx, anchor),
        jac=None,
        hess=None,
        constraints=constraints,
        tol=float(options.get("solver_tol", 1.0e-9)),
        options=_ipopt_options(options),
    )


def _candidate_from_result(sol, ctx: dict, start_index: int, anchored_retry: bool) -> dict | None:
    y = np.asarray(getattr(sol, "x", []), dtype=float).flatten()
    if y.size != ctx["ncomp"] or np.any(~np.isfinite(y)):
        return None
    unpacked = _unpack_y(y, ctx["z_feed"], ctx["beta_bounds"], x_floor=ctx["x_floor"])
    if not unpacked["valid"]:
        return None
    x1 = unpacked["x1"]
    x2 = unpacked["x2"]
    try:
        st1 = pcs_core._phase_state_liq(ctx["t"], ctx["p"], x1, ctx["params"])
        st2 = pcs_core._phase_state_liq(ctx["t"], ctx["p"], x2, ctx["params"])
    except Exception:
        return None
    residual = pcs_core._phase_equilibrium_residual(
        x1, x2, ctx["z"], ctx["E"], ctx["neutral_idx"], ctx["charged_idx"], st1["lnfug"], st2["lnfug"]
    )
    return {
        "sol": sol,
        "x1": x1,
        "x2": x2,
        "beta": float(unpacked["beta"]),
        "st1": st1,
        "st2": st2,
        "residual": np.asarray(residual, dtype=float),
        "residual_norm": float(np.linalg.norm(residual)),
        "split_norm": float(np.max(np.abs(x1 - x2))),
        "objective_value": float(getattr(sol, "fun", 0.5 * np.dot(residual, residual))),
        "solver_info": {
            "backend": "cyipopt",
            "executed": True,
            "status": int(getattr(sol, "status", -999) if getattr(sol, "status", None) is not None else -999),
            "message": str(getattr(sol, "message", "")),
            "success": bool(getattr(sol, "success", False)),
            "objective_value": float(getattr(sol, "fun", 0.5 * np.dot(residual, residual))),
            "nit": int(getattr(sol, "nit", 0) if getattr(sol, "nit", None) is not None else 0),
            "start_index": int(start_index),
            "anchored_retry": bool(anchored_retry),
        },
    }


def _stable_single_phase_result(feed_state: dict, z_feed: np.ndarray, tpdf: dict, E: np.ndarray, e_info: dict) -> dict:
    return {
        "n_phases": 1,
        "phases": [{
            "beta": 1.0,
            "x": np.asarray(z_feed, dtype=float),
            "rho": float(feed_state["rho"]),
            "lnfugcoef": np.asarray(feed_state["lnfugcoef"], dtype=float),
            "lnfug": np.asarray(feed_state["lnfug"], dtype=float),
        }],
        "tpdf_min": float(tpdf["tpdf_min"]),
        "tpdf_seed_x": np.asarray(tpdf["seed_x"], dtype=float),
        "converged": True,
        "status": 1,
        "message": "Feed is stable (single liquid phase).",
        "residual_norm": 0.0,
        "e_matrix": np.asarray(E, dtype=float),
        "ion_pair_rows": e_info["ion_pair_rows"],
        "charged_species": e_info["charged_species"],
        "charged_species_indices": e_info["charged_idx"],
        "solver_info": {"backend": "cyipopt", "executed": False, "reason": "single_phase_stable"},
    }


def solve_two_phase_lle_cyipopt(t, p, z_feed, params, species, options=None) -> dict:
    _require_cyipopt()
    if options is None:
        options = {}
    z_feed = np.asarray(z_feed, dtype=float).flatten()
    if len(species) != z_feed.size:
        raise pcs_core.InputError("`species` length ({}) must match z_feed length ({}).".format(len(species), z_feed.size))
    pcs_core.check_input(z_feed, {"temperature": t, "pressure": p})
    _, params = pcs_core.ensure_numpy_input(z_feed, params)
    params = pcs_core.check_association(params)
    z = np.asarray(params.get("z", []), dtype=float).flatten()
    if z.size != z_feed.size:
        raise pcs_core.InputError("params['z'] must be present and aligned with z_feed for solve_two_phase_lle_cyipopt.")
    if np.allclose(z, 0.0):
        raise pcs_core.InputError("solve_two_phase_lle_cyipopt requires ionic species (non-zero z).")
    feed_charge = float(np.dot(z, z_feed))
    if abs(feed_charge) > 1.0e-8:
        raise pcs_core.InputError("Feed must be electroneutral. Sum(z_i*z_feed_i) = {}".format(feed_charge))

    neutral_idx = np.where(np.abs(z) <= 1.0e-12)[0].astype(int)
    charged_idx = pcs_core._split_ion_groups(z, z_feed)["charged_idx"].astype(int)
    if charged_idx.size < 2:
        raise pcs_core.InputError("At least two charged species are required.")
    E, e_info = pcs_core._build_e_matrix(z, z_feed, species=species)
    tpdf = pcs_core._find_tpdf_seed(t, p, z_feed, params, neutral_idx, charged_idx, E, options)
    tpdf_tol = float(options.get("tpdf_tol", -1.0e-8))
    seed_x = options.get("seed_x", None)
    force_seed_solve = bool(options.get("force_seed_solve", False))
    seed_x = np.asarray(seed_x, dtype=float).flatten() if seed_x is not None else None
    if seed_x is not None:
        if seed_x.size != z_feed.size:
            raise pcs_core.InputError("options['seed_x'] length ({}) must match z_feed length ({}).".format(seed_x.size, z_feed.size))
        if np.any(seed_x <= 0.0):
            raise pcs_core.InputError("options['seed_x'] must contain only positive entries.")
        seed_x = seed_x / np.sum(seed_x)
    if tpdf["tpdf_min"] >= tpdf_tol and not (force_seed_solve and seed_x is not None):
        return _stable_single_phase_result(tpdf["feed_state"], z_feed, tpdf, E, e_info)

    beta_bounds = options.get("beta_bounds", (1.0e-8, 1.0 - 1.0e-8))
    if len(beta_bounds) != 2:
        raise ValueError("beta_bounds must be a two-element tuple.")
    beta_lo, beta_hi = float(beta_bounds[0]), float(beta_bounds[1])
    if not (0.0 < beta_lo < beta_hi < 1.0):
        raise ValueError("beta_bounds must satisfy 0 < low < high < 1.")
    charge_tol = float(options.get("charge_tol", 1.0e-6))
    mass_balance_tol = float(options.get("mass_balance_tol", 1.0e-8))
    split_tol = float(options.get("split_tol", 1.0e-3))
    solver_tol = float(options.get("solver_tol", 1.0e-9))
    solver_accept_norm = float(options.get("solver_accept_norm", 2.0e-1))
    ctx = {
        "t": float(t),
        "p": float(p),
        "z_feed": np.asarray(z_feed, dtype=float),
        "params": params,
        "z": np.asarray(z, dtype=float),
        "E": np.asarray(E, dtype=float),
        "neutral_idx": neutral_idx,
        "charged_idx": charged_idx,
        "charge_weight": float(options.get("charge_weight", 1000.0)),
        "beta_bounds": (beta_lo, beta_hi),
        "ncomp": int(z_feed.size),
        "x_floor": 1.0e-15,
        "invalid_scale": 1.0e6,
        "invalid_objective": 1.0e12,
    }

    trial_seed = seed_x if seed_x is not None else np.asarray(tpdf["seed_x"], dtype=float)
    starts = _build_starts(trial_seed, z_feed, int(options.get("random_seed", 12345)))
    beta0 = 0.5 * (beta_lo + beta_hi)
    candidates = []
    for start_index, x_start in enumerate(starts):
        try:
            cand = _candidate_from_result(_solve_start(_pack_start(x_start, beta0, (beta_lo, beta_hi)), ctx, options), ctx, start_index, False)
        except Exception:
            cand = None
        if cand is not None:
            candidates.append(cand)
    if len(candidates) == 0:
        return {
            "converged": False, "status": -1, "message": "cyipopt did not produce a candidate solution.", "residual_norm": np.inf,
            "phases": [], "n_phases": 2, "tpdf_min": float(tpdf["tpdf_min"]), "tpdf_seed_x": np.asarray(tpdf["seed_x"], dtype=float),
            "e_matrix": np.asarray(E, dtype=float), "ion_pair_rows": e_info["ion_pair_rows"], "charged_species": e_info["charged_species"],
            "charged_species_indices": e_info["charged_idx"], "solver_info": {"backend": "cyipopt", "executed": True, "status": -1, "message": "No valid candidate."},
        }

    best_res = min(c["residual_norm"] for c in candidates)
    split_candidates = [c for c in candidates if c["split_norm"] > 1.0e-3 and c["residual_norm"] <= max(1.0e-3, 50.0 * best_res)]
    best_cand = min(split_candidates if split_candidates else candidates, key=lambda c: (c["residual_norm"], c["objective_value"]))
    if best_cand["split_norm"] <= 1.0e-4:
        split_delta = np.abs(np.asarray(trial_seed, dtype=float) - np.asarray(z_feed, dtype=float))
        if split_delta.size > 0 and float(np.max(split_delta)) > 1.0e-8:
            anchor = {"component": int(np.argmax(split_delta)), "target": float(trial_seed[int(np.argmax(split_delta))]), "weight": 1.0}
            anchored = []
            for start_index, x_start in enumerate(starts):
                try:
                    cand = _candidate_from_result(_solve_start(_pack_start(x_start, beta0, (beta_lo, beta_hi)), ctx, options, anchor=anchor), ctx, start_index, True)
                except Exception:
                    cand = None
                if cand is not None and cand["split_norm"] > 1.0e-4:
                    anchored.append(cand)
            if anchored:
                best_cand = min(anchored, key=lambda c: (c["residual_norm"], c["objective_value"]))

    x1 = best_cand["x1"]
    x2 = best_cand["x2"]
    beta = float(best_cand["beta"])
    residual_norm = float(best_cand["residual_norm"])
    tol_ok = residual_norm <= max(solver_accept_norm, 1.0e3 * solver_tol)
    charge_ok = abs(float(np.dot(z, x1))) <= charge_tol and abs(float(np.dot(z, x2))) <= charge_tol
    mb_ok = np.max(np.abs(z_feed - (beta * x1 + (1.0 - beta) * x2))) <= mass_balance_tol
    frac_ok = (beta > beta_lo) and (beta < beta_hi)
    split_ok = best_cand["split_norm"] > split_tol
    return {
        "n_phases": 2,
        "phases": [
            {"beta": float(beta), "x": np.asarray(x1, dtype=float), "rho": float(best_cand["st1"]["rho"]), "lnfugcoef": np.asarray(best_cand["st1"]["lnfugcoef"], dtype=float), "lnfug": np.asarray(best_cand["st1"]["lnfug"], dtype=float)},
            {"beta": float(1.0 - beta), "x": np.asarray(x2, dtype=float), "rho": float(best_cand["st2"]["rho"]), "lnfugcoef": np.asarray(best_cand["st2"]["lnfugcoef"], dtype=float), "lnfug": np.asarray(best_cand["st2"]["lnfug"], dtype=float)},
        ],
        "tpdf_min": float(tpdf["tpdf_min"]),
        "tpdf_seed_x": np.asarray(tpdf["seed_x"], dtype=float),
        "converged": bool(charge_ok and mb_ok and frac_ok and split_ok and tol_ok),
        "status": int(best_cand["solver_info"]["status"]),
        "message": best_cand["solver_info"]["message"],
        "residual_norm": residual_norm,
        "e_matrix": np.asarray(E, dtype=float),
        "ion_pair_rows": e_info["ion_pair_rows"],
        "charged_species": e_info["charged_species"],
        "charged_species_indices": e_info["charged_idx"],
        "solver_info": {**best_cand["solver_info"], "split_norm": float(best_cand["split_norm"]), "charge_ok": bool(charge_ok), "mass_balance_ok": bool(mb_ok), "tol_ok": bool(tol_ok), "split_ok": bool(split_ok), "fraction_ok": bool(frac_ok)},
    }


def benchmark_ascani_case2(fd_eps: float | None = None) -> list[dict]:
    from scripts.multiphase_model_analysis import ascani_case2_dataset_comparison as case2

    rows = []
    for config in case2._default_model_configs():
        species, z_feed, _ = case2._case2_feed()
        params = case2._build_params_for_config(config, species, z_feed)
        current = case2._solve_lle_with_retries(case2.T_REF, case2.P_REF, z_feed, params, species)
        cy_options = dict(current.get("_solve_options", {}))
        if fd_eps is not None:
            cy_options["fd_eps"] = float(fd_eps)
        trial = solve_two_phase_lle_cyipopt(case2.T_REF, case2.P_REF, z_feed, params, species, options=cy_options)
        rows.append({"config_key": config["key"], "current": {"converged": bool(current.get("converged", False)), "status": int(current.get("status", -1)), "residual_norm": float(current.get("residual_norm", np.inf))}, "cyipopt": {"converged": bool(trial.get("converged", False)), "status": int(trial.get("status", -1)), "residual_norm": float(trial.get("residual_norm", np.inf)), "solver_info": trial.get("solver_info", {})}})
    return rows


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Experimental cyipopt two-phase electrolyte LLE solver.")
    parser.add_argument("--bench-ascani-case2", action="store_true", help="Run the current solver and the cyipopt prototype on the Ascani case-2 configs.")
    parser.add_argument("--fd-eps", type=float, default=None, help="Override the cyipopt finite-difference step.")
    args = parser.parse_args(argv)
    if args.bench_ascani_case2:
        print(json.dumps(benchmark_ascani_case2(fd_eps=args.fd_eps), indent=2))
    else:
        parser.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
