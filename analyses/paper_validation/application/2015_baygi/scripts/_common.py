from __future__ import annotations

import csv
import math
import os
import platform
import sys


from pathlib import Path
import sys as _bootstrap_sys
from pathlib import Path as _BootstrapPath

for _candidate in _BootstrapPath(__file__).resolve().parents:
    if (_candidate / "scripts" / "plot_outputs.py").is_file():
        if str(_candidate) not in _bootstrap_sys.path:
            _bootstrap_sys.path.insert(0, str(_candidate))
        break
else:
    raise ModuleNotFoundError("Could not locate repo root containing scripts/plot_outputs.py")
from scripts.plot_outputs import REPO_ROOT
from typing import Iterable
from typing import Mapping

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

platform.machine = lambda: os.environ.get("PROCESSOR_ARCHITECTURE", "AMD64")

ROOT = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts._env import require_epcsaft_install
from scripts.plot_outputs import paper_validation_output_path
from scripts.plot_outputs import save_plot_figure

require_epcsaft_install()

from epcsaft import SolutionError
from epcsaft import ePCSAFTMixture

MEA_MW = 0.06108
T_MIN = 303.15
T_MAX = 443.15
TEMPERATURE_GRID = np.linspace(T_MIN, T_MAX, 21)

TABLE1_VAPOR_PRESSURE = {"A": 92.624, "B": -10367.0, "C": -9.4699, "D": 1.9e-18, "E": 6.0}
TABLE1_LIQUID_DENSITY = {"A": 1.0011, "B": 0.22523, "C": 678.2, "D": 0.21515, "E": 0.0}

TABLE2_MEA_PARAMETERS = {
    "2B": {"m": 3.0353, "s": 3.0435, "e": 277.1740, "e_assoc": 2586.3, "vol_a": 0.037470},
    "3B": {"m": 4.5354, "s": 2.6019, "e": 204.0438, "e_assoc": 2383.4744, "vol_a": 0.118488},
    "4C": {"m": 4.5208, "s": 2.6574, "e": 237.6864, "e_assoc": 989.8984, "vol_a": 0.187533},
}
TABLE2_REPORTED_AAD = {
    "2B": {"P_Pa": 0.62, "rho_mol_m3": 0.12},
    "3B": {"P_Pa": 1.75, "rho_mol_m3": 0.43},
    "4C": {"P_Pa": 0.24, "rho_mol_m3": 0.23},
}

SCHEME_STYLES = {
    "2B": {"color": "#2ca02c", "linestyle": "-"},
    "3B": {"color": "#1f77b4", "linestyle": "--"},
    "4C": {"color": "#ff7f0e", "linestyle": "-."},
}


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 10,
            "axes.linewidth": 1.0,
            "axes.grid": False,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": False,
            "ytick.right": False,
            "legend.frameon": False,
            "mathtext.default": "regular",
        }
    )


def save_figure(fig: plt.Figure, path: Path) -> Path:
    path = paper_validation_output_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    save_plot_figure(fig, path, dpi=300, svg_companion=True)
    return path


def output_path(path: Path) -> Path:
    return paper_validation_output_path(path)


def baygi_output_root() -> Path:
    return ROOT.parent / "figures" / "regressed_parameters" / "output"


def regressed_parameters_path() -> Path:
    return baygi_output_root() / "regressed_parameters.csv"


def write_csv_rows(path: Path, fieldnames: Iterable[str], rows: Iterable[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames))
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def _figure_data_result_path(figure: str, filename: str) -> Path:
    return ROOT.parent / "figures" / figure / "output" / "data" / filename


def _diagnostics_cache_candidates(preferred_path: Path) -> list[Path]:
    candidates = [
        preferred_path,
        _figure_data_result_path("figure_2", "figure_2_diagnostics.csv"),
        _figure_data_result_path("figure_3", "figure_3_diagnostics.csv"),
    ]
    unique: list[Path] = []
    for candidate in candidates:
        if candidate not in unique:
            unique.append(candidate)
    return unique


def _regressed_diagnostics_cache_candidates(preferred_path: Path) -> list[Path]:
    candidates = [
        preferred_path,
        _figure_data_result_path("figure_2_regressed", "figure_2_regressed_diagnostics.csv"),
        _figure_data_result_path("figure_3_regressed", "figure_3_regressed_diagnostics.csv"),
    ]
    unique: list[Path] = []
    for candidate in candidates:
        if candidate not in unique:
            unique.append(candidate)
    return unique


def _rows_from_diagnostics(path: Path) -> tuple[list[dict], list[dict], list[dict]]:
    diagnostics = read_csv_rows(path)
    dippr = [row for row in diagnostics if row.get("series") == "DIPPR"]
    non_reference = [row for row in diagnostics if row.get("series") != "DIPPR"]
    model_rows = [row for row in non_reference if row.get("status") == "solved"]
    if not dippr or not model_rows:
        raise RuntimeError(f"Cached Baygi diagnostics are incomplete: {path}")
    return dippr, model_rows, non_reference


def saturation_rows_for_plot(
    diagnostics_path: Path,
    *,
    parameter_sets: Mapping[str, Mapping[str, float]] | None = None,
) -> tuple[list[dict], list[dict], list[dict]]:
    if os.environ.get("EPCSAFT_RECOMPUTE_BAYGI") != "1":
        candidates = (
            _regressed_diagnostics_cache_candidates(diagnostics_path)
            if parameter_sets is not None
            else _diagnostics_cache_candidates(diagnostics_path)
        )
        for candidate in candidates:
            if candidate.exists():
                return _rows_from_diagnostics(candidate)
    dippr = dippr_reference_rows()
    model_rows, diagnostics = model_saturation_rows(parameter_sets=parameter_sets)
    return dippr, model_rows, diagnostics


def baygi_mea_psat_pa(T: float) -> float:
    c = TABLE1_VAPOR_PRESSURE
    return math.exp(c["A"] + c["B"] / T + c["C"] * math.log(T) + c["D"] * T ** c["E"])


def baygi_mea_density_molar(T: float) -> float:
    c = TABLE1_LIQUID_DENSITY
    rho_mol_l = c["A"] / (c["B"] ** (1.0 + (1.0 - T / c["C"]) ** c["D"]))
    return 1000.0 * rho_mol_l


def dippr_reference_rows(temperatures: Iterable[float] = TEMPERATURE_GRID) -> list[dict[str, float | str]]:
    rows = []
    for T in temperatures:
        rows.append(
            {
                "series": "DIPPR",
                "T_K": float(T),
                "P_Pa": baygi_mea_psat_pa(float(T)),
                "rho_mol_m3": baygi_mea_density_molar(float(T)),
                "status": "reference",
                "contribution_terms": "reference",
            }
        )
    return rows


def build_mea_params(assoc_scheme: str, values: Mapping[str, float] | None = None) -> dict:
    values = dict(TABLE2_MEA_PARAMETERS[assoc_scheme] if values is None else values)
    params = {
        "m": np.asarray([values["m"]], dtype=float),
        "s": np.asarray([values["s"]], dtype=float),
        "e": np.asarray([values["e"]], dtype=float),
        "e_assoc": np.asarray([values["e_assoc"]], dtype=float),
        "vol_a": np.asarray([values["vol_a"]], dtype=float),
        "assoc_scheme": [assoc_scheme],
        "z": np.asarray([], dtype=float),
        "dielc": np.asarray([8.0], dtype=float),
        "d_born": np.asarray([0.0], dtype=float),
        "f_solv": np.asarray([1.0], dtype=float),
        "MW": np.asarray([MEA_MW], dtype=float),
        "k_ij": np.zeros((1, 1), dtype=float),
        "l_ij": np.zeros((1, 1), dtype=float),
        "k_hb": np.zeros((1, 1), dtype=float),
    }
    return params


def contribution_terms(params: dict) -> str:
    if np.asarray(params["z"], dtype=float).size != 0:
        raise RuntimeError("Baygi pure MEA validation must not include ion terms.")
    for key in ("k_ij", "l_ij", "k_hb"):
        if not np.allclose(np.asarray(params[key], dtype=float), 0.0):
            raise RuntimeError("Baygi pure MEA validation must not include binary interactions.")
    return "hc;disp;assoc"


def _fugacity_difference(mixture: ePCSAFTMixture, T: float, log_p: float) -> float:
    P = math.exp(float(log_p))
    x = np.asarray([1.0], dtype=float)
    mixture.clear_runtime_caches()
    liquid = mixture.state(T, x, P=P, phase="liq")
    mixture.clear_runtime_caches()
    vapor = mixture.state(T, x, P=P, phase="vap")
    return float(liquid.fugacity_coefficient(natural_log=True)[0] - vapor.fugacity_coefficient(natural_log=True)[0])


def _bracket_saturation(mixture: ePCSAFTMixture, T: float, p_hint: float) -> tuple[float, float] | None:
    guesses = np.clip(p_hint * np.exp(np.linspace(-2.0, 2.0, 33)), 1.0e-3, 1.0e8)
    samples: list[tuple[float, float]] = []
    for P in sorted(set(float(value) for value in guesses)):
        try:
            samples.append((math.log(P), _fugacity_difference(mixture, T, math.log(P))))
        except (SolutionError, ValueError, OverflowError, RuntimeError):
            continue
    for (x0, y0), (x1, y1) in zip(samples, samples[1:]):
        if not np.isfinite(y0) or not np.isfinite(y1):
            continue
        if y0 == 0.0:
            return x0, x0
        if y0 * y1 < 0.0:
            return x0, x1
    return None


def _bisect_log_pressure_root(
    mixture: ePCSAFTMixture,
    T: float,
    lower: float,
    upper: float,
    *,
    tolerance: float = 1.0e-10,
    max_iterations: int = 100,
) -> float:
    if lower == upper:
        return float(lower)
    f_lower = _fugacity_difference(mixture, T, lower)
    f_upper = _fugacity_difference(mixture, T, upper)
    if f_lower == 0.0:
        return float(lower)
    if f_upper == 0.0:
        return float(upper)
    if f_lower * f_upper > 0.0:
        raise RuntimeError("Saturation bracket does not contain a sign change.")

    lo = float(lower)
    hi = float(upper)
    flo = float(f_lower)
    for _ in range(max_iterations):
        mid = 0.5 * (lo + hi)
        f_mid = _fugacity_difference(mixture, T, mid)
        if abs(f_mid) <= tolerance or abs(hi - lo) <= tolerance:
            return float(mid)
        if flo * f_mid <= 0.0:
            hi = mid
        else:
            lo = mid
            flo = float(f_mid)
    return float(0.5 * (lo + hi))


def metric_rows(rows: list[dict], field: str) -> list[dict]:
    reference = {round(float(row["T_K"]), 6): float(row[field]) for row in rows if row["series"] == "DIPPR"}
    out: list[dict] = []
    for assoc_scheme in ("2B", "3B", "4C"):
        residuals = []
        for row in rows:
            if row["series"] != f"MEA {assoc_scheme}" or row["status"] != "solved":
                continue
            expected = reference[round(float(row["T_K"]), 6)]
            residuals.append(abs(float(row[field]) / expected - 1.0) * 100.0)
        reported = TABLE2_REPORTED_AAD[assoc_scheme][field]
        aad = float("nan") if not residuals else sum(residuals) / len(residuals)
        out.append(
            {
                "series": f"MEA {assoc_scheme}",
                "field": field,
                "n_solved": len(residuals),
                "aard_percent": aad,
                "reported_table2_aad_percent": reported,
                "within_reported_aad": bool(residuals and aad <= reported + 0.25),
            }
        )
    return out


def baygi_mea_fit_records(temperatures: Iterable[float] = TEMPERATURE_GRID) -> list[dict[str, float]]:
    return [
        {
            "T": float(T),
            "P": baygi_mea_psat_pa(float(T)),
            "rho": baygi_mea_density_molar(float(T)),
        }
        for T in temperatures
    ]


def _baygi_regression_bounds(values: Mapping[str, float]) -> dict[str, tuple[float, float]]:
    return {
        "m": (0.5 * float(values["m"]), 1.5 * float(values["m"])),
        "s": (0.8 * float(values["s"]), 1.2 * float(values["s"])),
        "e": (0.5 * float(values["e"]), 1.5 * float(values["e"])),
        "e_assoc": (0.5 * float(values["e_assoc"]), 1.5 * float(values["e_assoc"])),
        "vol_a": (0.5 * float(values["vol_a"]), 1.5 * float(values["vol_a"])),
    }


def regressed_parameter_rows(path: Path) -> list[dict[str, str]]:
    if os.environ.get("EPCSAFT_RECOMPUTE_BAYGI_REGRESSION") != "1" and path.exists():
        return read_csv_rows(path)

    from epcsaft.regression import _fit_pure_neutral_associating_native

    records = baygi_mea_fit_records()
    max_nfev = int(os.environ.get("EPCSAFT_BAYGI_REGRESSION_MAX_NFEV", "1"))
    rows: list[dict[str, str | float | int | bool]] = []
    for assoc_scheme, table_values in TABLE2_MEA_PARAMETERS.items():
        result = _fit_pure_neutral_associating_native(
            records,
            "MEA",
            assoc_scheme=assoc_scheme,
            fixed_parameters={"MW": MEA_MW},
            initial_guess=table_values,
            bounds=_baygi_regression_bounds(table_values),
            multistart=0,
            max_nfev=max_nfev,
        )
        row: dict[str, str | float | int | bool] = {
            "assoc_scheme": assoc_scheme,
            "backend": result.backend,
            "success": result.success,
            "status": result.status,
            "message": result.message,
            "nfev": result.nfev,
            "cost": result.cost,
            "residual_norm": result.residual_norm,
            "initial_residual_norm": result.metrics_by_term.get("initial_residual_norm", ""),
        }
        row.update(result.fitted_values)
        rows.append(row)

    write_csv_rows(
        path,
        (
            "assoc_scheme",
            "m",
            "s",
            "e",
            "e_assoc",
            "vol_a",
            "backend",
            "success",
            "status",
            "message",
            "nfev",
            "cost",
            "residual_norm",
            "initial_residual_norm",
        ),
        rows,
    )
    return read_csv_rows(path)


def parameter_sets_from_rows(rows: Iterable[Mapping[str, str]]) -> dict[str, dict[str, float]]:
    out: dict[str, dict[str, float]] = {}
    for row in rows:
        assoc_scheme = str(row["assoc_scheme"])
        out[assoc_scheme] = {
            "m": float(row["m"]),
            "s": float(row["s"]),
            "e": float(row["e"]),
            "e_assoc": float(row["e_assoc"]),
            "vol_a": float(row["vol_a"]),
        }
    return out


def model_saturation_rows(
    temperatures: Iterable[float] = TEMPERATURE_GRID,
    *,
    parameter_sets: Mapping[str, Mapping[str, float]] | None = None,
) -> tuple[list[dict], list[dict]]:
    rows: list[dict] = []
    diagnostics: list[dict] = []
    parameter_sets = TABLE2_MEA_PARAMETERS if parameter_sets is None else parameter_sets
    for assoc_scheme in ("2B", "3B", "4C"):
        params = build_mea_params(assoc_scheme, parameter_sets[assoc_scheme])
        terms = contribution_terms(params)
        for T in temperatures:
            mixture = ePCSAFTMixture.from_params(params, species=["MEA"])
            p_hint = baygi_mea_psat_pa(float(T))
            bracket = _bracket_saturation(mixture, float(T), p_hint)
            if bracket is None:
                diagnostics.append(
                    {
                        "series": f"MEA {assoc_scheme}",
                        "T_K": float(T),
                        "P_Pa": "",
                        "rho_mol_m3": "",
                        "status": "no_saturation_bracket",
                        "contribution_terms": terms,
                    }
                )
                continue
            try:
                log_p = _bisect_log_pressure_root(mixture, float(T), *bracket)
                P = math.exp(float(log_p))
                mixture.clear_runtime_caches()
                liquid = mixture.state(float(T), np.asarray([1.0], dtype=float), P=P, phase="liq")
                rho = float(liquid.molar_density())
            except (SolutionError, ValueError, OverflowError, RuntimeError) as exc:
                diagnostics.append(
                    {
                        "series": f"MEA {assoc_scheme}",
                        "T_K": float(T),
                        "P_Pa": "",
                        "rho_mol_m3": "",
                        "status": f"solve_failed:{type(exc).__name__}",
                        "contribution_terms": terms,
                    }
                )
                continue
            row = {
                "series": f"MEA {assoc_scheme}",
                "T_K": float(T),
                "P_Pa": P,
                "rho_mol_m3": rho,
                "status": "solved",
                "contribution_terms": terms,
            }
            rows.append(row)
            diagnostics.append(row)
    return rows, diagnostics
