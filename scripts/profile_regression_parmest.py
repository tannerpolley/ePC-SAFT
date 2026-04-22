"""Benchmark local Pyomo/parmest against native pure-neutral regression backends."""

from __future__ import annotations

import csv
import json
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

LOCAL_PYOMO_ROOT = Path(r"C:\Users\Tanner\Documents\git\pyomo")
if LOCAL_PYOMO_ROOT.exists() and str(LOCAL_PYOMO_ROOT) not in sys.path:
    sys.path.insert(0, str(LOCAL_PYOMO_ROOT))

import pandas as pd
import pyomo.environ as pyo
import pyomo.contrib.parmest.parmest as parmest
from pyomo.contrib.parmest.experiment import Experiment

from epcsaft import ePCSAFTMixture
from epcsaft import fit_pure_neutral
from epcsaft.regression import _debug_native_pure_neutral_objective
from epcsaft.regression import _fit_pure_neutral_ipopt_explicit_internal
from epcsaft.regression import _fit_pure_neutral_least_squares_internal
from tests.test_regression import _load_workbook_reference_rows
from tests.test_regression import _neutral_fixed_parameters
from tests.test_regression import _real_saturation_records


REPORT_DIR = REPO_ROOT / "build" / "runtime_profile"
REPORT_CSV = REPORT_DIR / "parmest_regression_profile.csv"
REPORT_MD = REPORT_DIR / "parmest_regression_profile.md"

PI = float(np.pi)
N_AV = 6.02214076e23
KB = 1.380649e-23
R_GAS = KB * N_AV
K_REGRESSION_GRADIENT_FLOOR = 1.0e-300

K_DISPERSION_A0 = [0.9105631445, 0.6361281449, 2.6861347891, -26.547362491, 97.759208784, -159.59154087, 91.297774084]
K_DISPERSION_A1 = [-0.3084016918, 0.1860531159, -2.5030047259, 21.419793629, -65.255885330, 83.318680481, -33.746922930]
K_DISPERSION_A2 = [-0.0906148351, 0.4527842806, 0.5962700728, -1.7241829131, -4.1302112531, 13.776631870, -8.6728470368]
K_DISPERSION_B0 = [0.7240946941, 2.2382791861, -4.0025849485, -21.003576815, 26.855641363, 206.55133841, -355.60235612]
K_DISPERSION_B1 = [-0.5755498075, 0.6995095521, 3.8925673390, -17.215471648, 192.67226447, -161.82646165, -165.20769346]
K_DISPERSION_B2 = [0.0976883116, -0.2557574982, -9.1558561530, 20.642075974, -38.804430052, 93.626774077, -29.666905585]


def _should_run_perf() -> bool:
    return os.environ.get("ePCSAFT_RUN_PERF", "").strip().lower() in {"1", "true", "yes", "on"}


def _case_list() -> tuple[str, ...]:
    raw = os.environ.get("ePCSAFT_PARMEST_CASES", "").strip()
    if not raw:
        return ("Methane", "Ethane", "Propane")
    cases = tuple(part.strip() for part in raw.split(",") if part.strip())
    if not cases:
        raise ValueError("ePCSAFT_PARMEST_CASES did not contain any usable component names.")
    return cases


def _rho_molar_from_record(component: str, record: dict[str, Any]) -> float:
    mw_value = float(_neutral_fixed_parameters(component)["MW"])
    return float(record["rho_sat_liq_kg_m3"]) / mw_value


def _candidate_kwargs(component: str) -> dict[str, Any]:
    refs = _load_workbook_reference_rows()
    ref = refs[component]
    return {
        "records": _real_saturation_records(component),
        "component": component,
        "assoc_scheme": "",
        "fixed_parameters": _neutral_fixed_parameters(component),
        "initial_guess": {
            "m": ref["m"] * 1.08,
            "s": ref["s"] * 0.96,
            "e": ref["e"] * 1.05,
        },
        "bounds": {
            "m": (0.5, 3.5),
            "s": (2.0, 5.0),
            "e": (50.0, 400.0),
        },
    }


def _pure_neutral_expr(model: pyo.ConcreteModel, rho, temperature: float):
    mseg = model.m
    sigma = model.s
    epsilon = model.e
    den = rho * (N_AV / 1.0e30)
    d = sigma * (1.0 - 0.12 * pyo.exp(-3.0 * epsilon / temperature))
    d2 = d * d
    d3 = d2 * d
    prefactor = (PI / 6.0) * den * mseg
    zeta0 = prefactor
    zeta1 = prefactor * d
    zeta2 = prefactor * d2
    zeta3 = prefactor * d3
    eta = zeta3
    pair_diameter = d / 2.0
    pair_diameter2 = pair_diameter * pair_diameter
    zeta2_sq = zeta2 * zeta2
    zeta2_cu = zeta2_sq * zeta2
    zeta3_sq = zeta3 * zeta3
    ghs = (
        1.0 / (1.0 - zeta3)
        + pair_diameter * 3.0 * zeta2 / (1.0 - zeta3) ** 2
        + pair_diameter2 * 2.0 * zeta2_sq / (1.0 - zeta3) ** 3
    )
    ares_hs = (
        3.0 * zeta1 * zeta2 / (1.0 - zeta3)
        + zeta2_cu / (zeta3 * (1.0 - zeta3) ** 2)
        + ((zeta2_cu / zeta3_sq) - zeta0) * pyo.log(1.0 - zeta3)
    ) / zeta0
    ares_hc = mseg * ares_hs - (mseg - 1.0) * pyo.log(ghs)

    c1 = (mseg - 1.0) / mseg
    c2 = (mseg - 2.0) / mseg
    I1 = 0.0
    I2 = 0.0
    d_eta_I1_deta = 0.0
    d_eta_I2_deta = 0.0
    for i in range(7):
        a_i = K_DISPERSION_A0[i] + c1 * K_DISPERSION_A1[i] + c1 * c2 * K_DISPERSION_A2[i]
        b_i = K_DISPERSION_B0[i] + c1 * K_DISPERSION_B1[i] + c1 * c2 * K_DISPERSION_B2[i]
        I1 += a_i * eta ** i
        I2 += b_i * eta ** i
        d_eta_I1_deta += a_i * float(i + 1) * eta ** i
        d_eta_I2_deta += b_i * float(i + 1) * eta ** i

    C1 = 1.0 / (
        1.0
        + mseg * (8.0 * eta - 2.0 * eta * eta) / (1.0 - eta) ** 4
        + (1.0 - mseg) * (20.0 * eta - 27.0 * eta * eta + 12.0 * eta ** 3 - 2.0 * eta ** 4)
        / ((1.0 - eta) * (2.0 - eta)) ** 2
    )
    C2 = -C1 * C1 * (
        mseg * (-4.0 * eta * eta + 20.0 * eta + 8.0) / (1.0 - eta) ** 5
        + (1.0 - mseg) * (2.0 * eta ** 3 + 12.0 * eta * eta - 48.0 * eta + 40.0)
        / ((1.0 - eta) * (2.0 - eta)) ** 3
    )

    e_over_t = epsilon / temperature
    m2es3 = mseg * mseg * e_over_t * sigma ** 3
    m2e2s3 = mseg * mseg * e_over_t * e_over_t * sigma ** 3
    ares_disp = -2.0 * PI * den * I1 * m2es3 - PI * den * mseg * C1 * I2 * m2e2s3
    ares_total = ares_hc + ares_disp

    dghs_drho = (
        zeta3 / (1.0 - zeta3) ** 2
        + pair_diameter * (3.0 * zeta2 / (1.0 - zeta3) ** 2 + 6.0 * zeta2 * zeta3 / (1.0 - zeta3) ** 3)
        + pair_diameter2 * (4.0 * zeta2_sq / (1.0 - zeta3) ** 3 + 6.0 * zeta2_sq * zeta3 / (1.0 - zeta3) ** 4)
    )
    dadrho_hs = (
        zeta3 / (1.0 - zeta3)
        + 3.0 * zeta1 * zeta2 / zeta0 / (1.0 - zeta3) ** 2
        + (3.0 * zeta2_cu - zeta3 * zeta2_cu) / zeta0 / (1.0 - zeta3) ** 3
    )
    zraw_hc = mseg * dadrho_hs - (mseg - 1.0) * dghs_drho / ghs
    zraw_disp = -2.0 * PI * den * d_eta_I1_deta * m2es3 - PI * den * mseg * (C1 * d_eta_I2_deta + C2 * eta * I2) * m2e2s3
    zraw_total = zraw_hc + zraw_disp
    Z = 1.0 + zraw_total
    pressure = Z * R_GAS * temperature * rho
    lnfug = ares_total + zraw_total - pyo.log(Z)
    return {"pressure": pressure, "lnfug": lnfug}


@dataclass
class DensityExperimentData:
    component: str
    temperature: float
    pressure: float
    rho_exp: float
    initial_guess: dict[str, float]
    bounds: dict[str, tuple[float, float]]
    density_scale: float


@dataclass
class VLEExperimentData:
    component: str
    temperature: float
    pressure: float
    rho_liq_guess: float
    rho_vap_guess: float
    initial_guess: dict[str, float]
    bounds: dict[str, tuple[float, float]]
    pure_vle_scale: float


def _square_seed_densities(component: str, temperature: float, pressure: float, theta: dict[str, float]) -> tuple[float, float]:
    params = dict(_neutral_fixed_parameters(component))
    params["m"] = np.asarray([theta["m"]], dtype=float)
    params["s"] = np.asarray([theta["s"]], dtype=float)
    params["e"] = np.asarray([theta["e"]], dtype=float)
    params["MW"] = np.asarray([float(params["MW"])], dtype=float)
    params["e_assoc"] = np.asarray([float(params["e_assoc"])], dtype=float)
    params["vol_a"] = np.asarray([float(params["vol_a"])], dtype=float)
    params["z"] = np.asarray([float(params["z"])], dtype=float)
    params["dielc"] = np.asarray([float(params["dielc"])], dtype=float)
    params["d_born"] = np.asarray([float(params["d_born"])], dtype=float)
    params["f_solv"] = np.asarray([float(params["f_solv"])], dtype=float)
    mixture = ePCSAFTMixture.from_params(params, species=[component])
    x = np.asarray([1.0], dtype=float)
    rho_liq = float(mixture.state(temperature, x, P=pressure, phase="liq").density())
    rho_vap = float(mixture.state(temperature, x, P=pressure, phase="vap").density())
    return rho_liq, rho_vap


class PureNeutralDensityExperiment(Experiment):
    def __init__(self, data: DensityExperimentData):
        self.data = data
        self.model = None

    def create_model(self):
        data = self.data
        model = pyo.ConcreteModel()
        model.m = pyo.Var(initialize=data.initial_guess["m"], bounds=data.bounds["m"], within=pyo.PositiveReals)
        model.s = pyo.Var(initialize=data.initial_guess["s"], bounds=data.bounds["s"], within=pyo.PositiveReals)
        model.e = pyo.Var(initialize=data.initial_guess["e"], bounds=data.bounds["e"], within=pyo.PositiveReals)
        rho_lower = max(1.0e-9, data.rho_exp * 0.25)
        rho_upper = max(rho_lower * 1.01, data.rho_exp * 4.0)
        model.rho = pyo.Var(initialize=data.rho_exp, bounds=(rho_lower, rho_upper), within=pyo.PositiveReals)
        model.scaled_output = pyo.Var(initialize=data.density_scale)
        expr = _pure_neutral_expr(model, model.rho, data.temperature)
        denom = max(abs(data.rho_exp), K_REGRESSION_GRADIENT_FLOOR)
        model.pressure_closure = pyo.Constraint(expr=expr["pressure"] == data.pressure)
        model.output_link = pyo.Constraint(expr=model.scaled_output == data.density_scale * model.rho / denom)
        self.model = model
        return model

    def finalize_model(self):
        return self.model

    def label_model(self):
        data = self.data
        model = self.model
        denom = max(abs(data.rho_exp), K_REGRESSION_GRADIENT_FLOOR)
        target_value = data.density_scale * data.rho_exp / denom
        model.experiment_outputs = pyo.Suffix(direction=pyo.Suffix.LOCAL)
        model.experiment_outputs.update([(model.scaled_output, target_value)])
        model.unknown_parameters = pyo.Suffix(direction=pyo.Suffix.LOCAL)
        model.unknown_parameters.update((var, pyo.ComponentUID(var)) for var in (model.m, model.s, model.e))
        return model

    def get_labeled_model(self):
        self.create_model()
        self.finalize_model()
        return self.label_model()


class PureNeutralVLEExperiment(Experiment):
    def __init__(self, data: VLEExperimentData):
        self.data = data
        self.model = None

    def create_model(self):
        data = self.data
        model = pyo.ConcreteModel()
        model.m = pyo.Var(initialize=data.initial_guess["m"], bounds=data.bounds["m"], within=pyo.PositiveReals)
        model.s = pyo.Var(initialize=data.initial_guess["s"], bounds=data.bounds["s"], within=pyo.PositiveReals)
        model.e = pyo.Var(initialize=data.initial_guess["e"], bounds=data.bounds["e"], within=pyo.PositiveReals)
        rho_liq_guess = max(data.rho_liq_guess, 1.0e-6)
        rho_vap_guess = max(data.rho_vap_guess, 1.0e-9)
        rho_liq_lower = max(1.0e-9, rho_liq_guess * 0.25)
        rho_liq_upper = max(rho_liq_lower * 1.01, rho_liq_guess * 4.0)
        rho_vap_lower = max(1.0e-12, rho_vap_guess * 0.05)
        rho_vap_upper = max(rho_vap_lower * 1.01, min(rho_vap_guess * 25.0, rho_liq_guess * 0.8))
        rho_liq_lower = max(rho_liq_lower, rho_vap_upper * 1.2)
        rho_liq_upper = max(rho_liq_upper, rho_liq_lower * 1.01)
        model.rho_liq = pyo.Var(initialize=rho_liq_guess, bounds=(rho_liq_lower, rho_liq_upper), within=pyo.PositiveReals)
        model.rho_vap = pyo.Var(initialize=rho_vap_guess, bounds=(rho_vap_lower, rho_vap_upper), within=pyo.PositiveReals)
        model.scaled_gap = pyo.Var(initialize=0.0)
        liq_expr = _pure_neutral_expr(model, model.rho_liq, data.temperature)
        vap_expr = _pure_neutral_expr(model, model.rho_vap, data.temperature)
        model.liq_pressure_closure = pyo.Constraint(expr=liq_expr["pressure"] == data.pressure)
        model.vap_pressure_closure = pyo.Constraint(expr=vap_expr["pressure"] == data.pressure)
        model.output_link = pyo.Constraint(expr=model.scaled_gap == data.pure_vle_scale * (liq_expr["lnfug"] - vap_expr["lnfug"]))
        self.model = model
        return model

    def finalize_model(self):
        return self.model

    def label_model(self):
        model = self.model
        model.experiment_outputs = pyo.Suffix(direction=pyo.Suffix.LOCAL)
        model.experiment_outputs.update([(model.scaled_gap, 0.0)])
        model.unknown_parameters = pyo.Suffix(direction=pyo.Suffix.LOCAL)
        model.unknown_parameters.update((var, pyo.ComponentUID(var)) for var in (model.m, model.s, model.e))
        return model

    def get_labeled_model(self):
        self.create_model()
        self.finalize_model()
        return self.label_model()


def _parmest_experiments(component: str, initial_guess: dict[str, float], bounds: dict[str, tuple[float, float]]) -> list[Experiment]:
    experiments: list[Experiment] = []
    records = _real_saturation_records(component)
    density_scale = 1.0
    pure_vle_scale = 1.0
    for record in records:
        rho_exp = _rho_molar_from_record(component, record)
        rho_liq_guess, rho_vap_guess = _square_seed_densities(
            component,
            float(record["T"]),
            float(record["P"]),
            initial_guess,
        )
        experiments.append(
            PureNeutralDensityExperiment(
                DensityExperimentData(
                    component=component,
                    temperature=float(record["T"]),
                    pressure=float(record["P"]),
                    rho_exp=rho_exp,
                    initial_guess=initial_guess,
                    bounds=bounds,
                    density_scale=density_scale,
                )
            )
        )
        experiments.append(
            PureNeutralVLEExperiment(
                VLEExperimentData(
                    component=component,
                    temperature=float(record["T"]),
                    pressure=float(record["P"]),
                    rho_liq_guess=rho_liq_guess,
                    rho_vap_guess=rho_vap_guess,
                    initial_guess=initial_guess,
                    bounds=bounds,
                    pure_vle_scale=pure_vle_scale,
                )
            )
        )
    return experiments


def _theta_dataframe(initial_guess: dict[str, float]) -> pd.DataFrame:
    return pd.DataFrame([{"m": initial_guess["m"], "s": initial_guess["s"], "e": initial_guess["e"]}])


def _benchmark_local_parmest(component: str) -> dict[str, Any]:
    kwargs = _candidate_kwargs(component)
    experiments = _parmest_experiments(component, kwargs["initial_guess"], kwargs["bounds"])
    pest = parmest.Estimator(
        experiments,
        obj_function="SSE",
        solver_options={
            "bound_push": 1.0e-8,
            "bound_frac": 1.0e-8,
            "tol": 1.0e-6,
            "acceptable_tol": 1.0e-5,
            "acceptable_iter": 3,
            "print_level": 0,
            "max_iter": 2000,
        },
    )
    theta0 = _theta_dataframe(kwargs["initial_guess"])
    t0 = time.perf_counter()
    init_obj_value = float("nan")
    init_message = ""
    try:
        init_obj = pest.objective_at_theta(theta0, initialize_parmest_model=True)
        init_obj_value = float(init_obj.iloc[0]["obj"])
    except Exception as exc:  # noqa: BLE001
        detail = str(exc).strip() or exc.__class__.__name__
        init_message = f"objective_at_theta init failed: {detail}"
    try:
        obj, theta = pest.theta_est()
        success = True
        message = "local pyomo parmest"
    except Exception as exc:  # noqa: BLE001
        elapsed = time.perf_counter() - t0
        message = str(exc).strip() or exc.__class__.__name__
        if init_message:
            message = f"{init_message}; theta_est failed: {message}"
        return {
            "case": component,
            "backend": "pyomo_parmest_local",
            "returned_backend": "pyomo_parmest_local",
            "workflow_selected": "pyomo_parmest_local",
            "wall_s": float(elapsed),
            "nfev": int(0),
            "success": False,
            "status": 1,
            "message": message,
            "m": float("nan"),
            "s": float("nan"),
            "e": float("nan"),
            "density_rms": float("nan"),
            "pure_vle_rms": float("nan"),
            "starts_tried": int(1),
            "fallback_triggered": False,
            "initial_cost": init_obj_value,
        }
    elapsed = time.perf_counter() - t0
    theta_map = {str(k): float(v) for k, v in theta.items()}
    debug = _debug_native_pure_neutral_objective(
        kwargs["records"],
        component,
        assoc_scheme=kwargs["assoc_scheme"],
        fixed_parameters=kwargs["fixed_parameters"],
        initial_guess=kwargs["initial_guess"],
        bounds=kwargs["bounds"],
        x=theta_map,
    )
    if init_message:
        message = f"{message}; {init_message}"
    return {
        "case": component,
        "backend": "pyomo_parmest_local",
        "returned_backend": "pyomo_parmest_local",
        "workflow_selected": "pyomo_parmest_local",
        "wall_s": float(elapsed),
        "nfev": int(0),
        "success": success,
        "status": 0,
        "message": message,
        "m": theta_map["m"],
        "s": theta_map["s"],
        "e": theta_map["e"],
        "density_rms": _rms(debug["density_raw_residuals"]),
        "pure_vle_rms": _rms(debug["pure_vle_raw_residuals"]),
        "starts_tried": int(1),
        "fallback_triggered": False,
        "initial_cost": init_obj_value,
    }


def _benchmark_native(component: str, backend: str) -> dict[str, Any]:
    kwargs = _candidate_kwargs(component)
    if backend == "public_default":
        solve = fit_pure_neutral
    elif backend == "least_squares_native":
        solve = _fit_pure_neutral_least_squares_internal
    elif backend == "ipopt_explicit_native":
        solve = _fit_pure_neutral_ipopt_explicit_internal
    else:
        raise ValueError(f"Unsupported backend: {backend}")
    t0 = time.perf_counter()
    result = solve(**kwargs)
    elapsed = time.perf_counter() - t0
    return {
        "case": component,
        "backend": backend,
        "returned_backend": str(result.backend),
        "workflow_selected": str(result.backend),
        "wall_s": float(elapsed),
        "nfev": int(result.nfev),
        "success": bool(result.success),
        "status": int(result.status),
        "message": str(result.message),
        "m": float(result.fitted_values["m"]),
        "s": float(result.fitted_values["s"]),
        "e": float(result.fitted_values["e"]),
        "density_rms": float(result.metrics_by_term["density"]),
        "pure_vle_rms": float(result.metrics_by_term["pure_vle_fugacity_balance"]),
        "starts_tried": int(0),
        "fallback_triggered": False,
        "initial_cost": float("nan"),
    }


def _format_float(value: Any) -> str:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return str(value)
    if numeric != numeric:
        return "nan"
    return f"{numeric:.6g}"


def _rms(values: Any) -> float:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return 0.0
    return float(np.sqrt(np.mean(arr * arr)))


def _write_reports(rows: list[dict[str, Any]]) -> None:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "case",
        "backend",
        "returned_backend",
        "workflow_selected",
        "wall_s",
        "nfev",
        "success",
        "status",
        "message",
        "m",
        "s",
        "e",
        "density_rms",
        "pure_vle_rms",
        "starts_tried",
        "fallback_triggered",
        "initial_cost",
    ]
    with REPORT_CSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})

    lines = [
        "# Local Pyomo/Parmest Regression Profile",
        "",
        "| Case | Backend | Wall (s) | NFEV | Success | Density RMS | Pure VLE RMS | Initial Cost | m | s | e |",
        "| --- | --- | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row["case"]),
                    str(row["backend"]),
                    _format_float(row["wall_s"]),
                    str(row["nfev"]),
                    "yes" if row["success"] else "no",
                    _format_float(row["density_rms"]),
                    _format_float(row["pure_vle_rms"]),
                    _format_float(row["initial_cost"]),
                    _format_float(row["m"]),
                    _format_float(row["s"]),
                    _format_float(row["e"]),
                ]
            )
            + " |"
        )
    REPORT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_parmest_runtime_profile() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for component in _case_list():
        rows.append(_benchmark_native(component, "public_default"))
        rows.append(_benchmark_native(component, "least_squares_native"))
        rows.append(_benchmark_native(component, "ipopt_explicit_native"))
        rows.append(_benchmark_local_parmest(component))
    _write_reports(rows)
    return rows


def main() -> int:
    if not _should_run_perf():
        print("Set ePCSAFT_RUN_PERF=1 to run the local pyomo/parmest regression profile.")
        return 0
    rows = run_parmest_runtime_profile()
    print(json.dumps(rows, indent=2))
    print(f"Wrote {REPORT_CSV}")
    print(f"Wrote {REPORT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
