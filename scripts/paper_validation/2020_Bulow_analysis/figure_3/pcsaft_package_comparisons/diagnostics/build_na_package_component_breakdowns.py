from __future__ import annotations

import csv
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
PACKAGE_DIR = SCRIPT_DIR.parent
FIGURE3_DIR = PACKAGE_DIR.parent
ANALYSIS_ROOT = FIGURE3_DIR.parent
REPO_ROOT = ANALYSIS_ROOT.parents[2]
FIGURE3_DIAG_DIR = FIGURE3_DIR / "diagnostics"

if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(FIGURE3_DIAG_DIR) not in sys.path:
    sys.path.insert(0, str(FIGURE3_DIAG_DIR))
if str(PACKAGE_DIR) not in sys.path:
    sys.path.insert(0, str(PACKAGE_DIR))

import _model_overlay as overlay
from data.epcsaft_properties import get_prop_dict
import figure3_detailed_bookkeeping as repo_diag
import feos_extractor
from pcsaft import pcsaft_lnfugcoef_terms

try:
    import feos
    import si_units as si
except ImportError as exc:  # pragma: no cover
    raise RuntimeError("feos and si_units must be importable in the PC-SAFT environment.") from exc


CLAPEYRON_ROOT = Path(r"C:\Users\Tanner\Documents\git\Clapeyron.jl")
CLAPEYRON_JSON = SCRIPT_DIR / "clapeyron_na_component_breakdown.json"
CLAPEYRON_SCRIPT = SCRIPT_DIR / "extract_clapeyron_na_component_breakdown.jl"

OUTPUT_COMMON = SCRIPT_DIR / "na_package_state_summary.csv"
OUTPUTS = {
    "hc": SCRIPT_DIR / "na_hc_package_breakdown.csv",
    "disp": SCRIPT_DIR / "na_disp_package_breakdown.csv",
    "assoc": SCRIPT_DIR / "na_assoc_package_breakdown.csv",
    "dh": SCRIPT_DIR / "na_dh_package_breakdown.csv",
    "born": SCRIPT_DIR / "na_born_package_breakdown.csv",
}
REPORT_PATH = SCRIPT_DIR / "na_package_component_breakdown_notes.md"

ION = "Na+"
SPECIES = ["Na+", "Cl-", "Water"]
PCSAFT_TARGET = 0
PCSAFT_COUNTER = 1
PCSAFT_WATER = 2
FEOS_WATER = 0
FEOS_TARGET = 1
FEOS_COUNTER = 2
FD_H = 1.0e-8

FEOS_COMPONENTS = feos_extractor.ION_SETUPS[ION].components
FEOS_PARAMETERS = feos.Parameters.from_json(
    list(FEOS_COMPONENTS),
    str(feos_extractor.FEOS_PARAMETER_PATH),
    str(feos_extractor.FEOS_BINARY_PATH),
)
FEOS_EOS = feos.EquationOfState.epcsaft(FEOS_PARAMETERS, epcsaft_variant="advanced")
FRAME = repo_diag.common.load_indexed_csv(repo_diag.DATA_PATH)


def _to_kj_per_mol(quantity) -> float:
    return float(quantity / (si.KILO * si.JOULE / si.MOL))


def _to_pa(quantity) -> float:
    return float(quantity / si.PASCAL)


def _label_slug(label: str) -> str:
    return label.lower().replace(" ", "_").replace("-", "_")


def _float_or_nan(value) -> float:
    if value is None:
        return float("nan")
    return float(value)


def _forward_simplex_fd(x_base: np.ndarray, idx: int, bal_idx: int, evaluator, h: float = FD_H) -> float:
    step = float(h)
    available = float(x_base[bal_idx])
    if available <= 0.0:
        return float("nan")
    if step >= 0.5 * available:
        step = 0.5 * available
    if step <= 0.0:
        return float("nan")
    x_new = x_base.copy()
    x_new[idx] += step
    x_new[bal_idx] -= step
    if np.any(x_new < 0.0):
        raise ValueError(f"Invalid simplex perturbation for idx={idx}, bal_idx={bal_idx}: {x_new}")
    return float((evaluator(x_new) - evaluator(x_base)) / step)


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _run_clapeyron_json() -> dict:
    cmd = [
        "julia",
        f"--project={CLAPEYRON_ROOT}",
        str(CLAPEYRON_SCRIPT),
        str(CLAPEYRON_JSON),
        str(CLAPEYRON_ROOT),
    ]
    completed = subprocess.run(
        cmd,
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    if completed.stdout.strip():
        print(completed.stdout.strip(), flush=True)
    if completed.stderr.strip():
        print(completed.stderr.strip(), flush=True)
    return json.loads(CLAPEYRON_JSON.read_text(encoding="utf-8"))


def _common_fields(
    package: str,
    config: str,
    component_order: str,
    target_component: str,
    counter_component: str,
    water_component: str,
    x_target: float,
    x_counter: float,
    x_water: float,
    density_mol_m3: float,
    compressibility_factor: float,
    sigma_target: float,
    sigma_counter: float,
    sigma_water: float,
    epsilon_target: float,
    epsilon_counter: float,
    epsilon_water: float,
    born_radius_target: float,
    born_radius_counter: float,
    water_target_k: float,
    water_counter_k: float,
    target_counter_k: float,
    dielectric_target: float,
    dielectric_counter: float,
    dielectric_water: float,
    **extras: object,
) -> dict[str, object]:
    row = {
        "package": package,
        "config": config,
        "component_order": component_order,
        "target_component": target_component,
        "counter_component": counter_component,
        "water_component": water_component,
        "x_target": x_target,
        "x_counter": x_counter,
        "x_water": x_water,
        "density_mol_m3": density_mol_m3,
        "compressibility_factor": compressibility_factor,
        "sigma_target": sigma_target,
        "sigma_counter": sigma_counter,
        "sigma_water": sigma_water,
        "epsilon_target": epsilon_target,
        "epsilon_counter": epsilon_counter,
        "epsilon_water": epsilon_water,
        "born_radius_target": born_radius_target,
        "born_radius_counter": born_radius_counter,
        "water_target_k": water_target_k,
        "water_counter_k": water_counter_k,
        "target_counter_k": target_counter_k,
        "dielectric_target": dielectric_target,
        "dielectric_counter": dielectric_counter,
        "dielectric_water": dielectric_water,
    }
    row.update(extras)
    return row


def _pcsaft_base():
    state = repo_diag._state_for_ion(ION)
    geom = repo_diag._mixture_geometry(state)
    params = state["params"]
    terms = state["terms"]
    x = np.asarray(state["x_inf"], dtype=float)
    rho = float(state["rho_inf"])

    def term_state(x_new: np.ndarray):
        params_new = get_prop_dict("bulow_2020", SPECIES, x_new, overlay.T_REF)
        return pcsaft_lnfugcoef_terms(overlay.T_REF, rho, x_new, params_new)

    common = _common_fields(
        package="pcsaft",
        config="bulow_2020",
        component_order="|".join(SPECIES),
        target_component="Na+",
        counter_component="Cl-",
        water_component="Water",
        x_target=float(x[PCSAFT_TARGET]),
        x_counter=float(x[PCSAFT_COUNTER]),
        x_water=float(x[PCSAFT_WATER]),
        density_mol_m3=rho,
        compressibility_factor=float(terms["z_total"]),
        sigma_target=float(params["s"][PCSAFT_TARGET]),
        sigma_counter=float(params["s"][PCSAFT_COUNTER]),
        sigma_water=float(params["s"][PCSAFT_WATER]),
        epsilon_target=float(params["e"][PCSAFT_TARGET]),
        epsilon_counter=float(params["e"][PCSAFT_COUNTER]),
        epsilon_water=float(params["e"][PCSAFT_WATER]),
        born_radius_target=float(params["d_born"][PCSAFT_TARGET]),
        born_radius_counter=float(params["d_born"][PCSAFT_COUNTER]),
        water_target_k=float(params["k_ij"][PCSAFT_WATER][PCSAFT_TARGET]),
        water_counter_k=float(params["k_ij"][PCSAFT_WATER][PCSAFT_COUNTER]),
        target_counter_k=float(params["k_ij"][PCSAFT_TARGET][PCSAFT_COUNTER]),
        dielectric_target=float(params.get("dielc", np.zeros(3))[PCSAFT_TARGET]),
        dielectric_counter=float(params.get("dielc", np.zeros(3))[PCSAFT_COUNTER]),
        dielectric_water=float(params.get("dielc", np.zeros(3))[PCSAFT_WATER]),
    )
    return state, geom, terms, x, common, term_state


def _pcsaft_exact_extra(contribution: str, state, geom) -> dict[str, object]:
    if contribution == "hc":
        row = next(row for row in repo_diag._hc_breakdown(state, geom, FRAME) if row["ion"] == ION)
        out = {f"pcsaft_exact_{key}": value for key, value in row.items() if key != "ion"}
        out["pcsaft_geom_m_avg"] = float(geom["m_avg"])
        out["pcsaft_geom_eta"] = float(geom["eta"])
        out["pcsaft_geom_ares_hs_kj_mol"] = repo_diag._kjmol(float(geom["ares_hs"]))
        out["pcsaft_geom_zhs_kj_mol"] = repo_diag._kjmol(float(geom["zhs"]))
        return out
    if contribution == "disp":
        row = next(row for row in repo_diag._disp_breakdown(state, geom, FRAME) if row["ion"] == ION)
        out = {f"pcsaft_exact_{key}": value for key, value in row.items() if key != "ion"}
        out["pcsaft_geom_m_avg"] = float(geom["m_avg"])
        out["pcsaft_geom_eta"] = float(geom["eta"])
        out["pcsaft_geom_i1"] = float(geom["i1"])
        out["pcsaft_geom_i2"] = float(geom["i2"])
        out["pcsaft_geom_c1"] = float(geom["c1"])
        out["pcsaft_geom_c2"] = float(geom["c2"])
        out["pcsaft_geom_m2es3"] = float(geom["m2es3"])
        out["pcsaft_geom_m2e2s3"] = float(geom["m2e2s3"])
        return out
    if contribution == "assoc":
        rows = repo_diag._assoc_breakdown(state, geom, FRAME)
        if not rows:
            return {}
        row = next(row for row in rows if row["ion"] == ION)
        return {f"pcsaft_exact_{key}": value for key, value in row.items() if key != "ion"}
    return {}


def _pcsaft_row(contribution: str, state, geom, terms, x, common, term_state) -> dict[str, object]:
    suffix = repo_diag.CONTRIBUTION_MAP[contribution]["suffix"]
    mu_arr = np.asarray(terms[f"mu_{suffix}"], dtype=float)
    lnfug_arr = np.asarray(terms[f"lnfugcoef_{suffix}"], dtype=float)
    a_term = float(terms[f"a_{suffix}"])
    z_term = float(terms[f"z_{suffix}"])
    row = dict(common)
    row.update(
        {
            "contribution": contribution,
            "a_kj_mol": repo_diag._kjmol(a_term),
            "z_kj_mol": repo_diag._kjmol(z_term),
            "a_plus_z_kj_mol": repo_diag._kjmol(a_term + z_term),
            "target_mu_kj_mol": repo_diag._kjmol(float(mu_arr[PCSAFT_TARGET])),
            "counter_mu_kj_mol": repo_diag._kjmol(float(mu_arr[PCSAFT_COUNTER])),
            "water_mu_kj_mol": repo_diag._kjmol(float(mu_arr[PCSAFT_WATER])),
            "target_lnfug_kj_mol": repo_diag._kjmol(float(lnfug_arr[PCSAFT_TARGET])),
            "target_z_correction_kj_mol": repo_diag._kjmol(float(lnfug_arr[PCSAFT_TARGET] - mu_arr[PCSAFT_TARGET])),
            "weighted_mu_kj_mol": repo_diag._kjmol(float(np.dot(x, mu_arr))),
            "weighted_mu_minus_a_plus_z_kj_mol": repo_diag._kjmol(float(np.dot(x, mu_arr) - (a_term + z_term))),
            "target_mu_minus_a_z_kj_mol": repo_diag._kjmol(float(mu_arr[PCSAFT_TARGET] - a_term - z_term)),
            "counter_mu_minus_a_z_kj_mol": repo_diag._kjmol(float(mu_arr[PCSAFT_COUNTER] - a_term - z_term)),
            "water_mu_minus_a_z_kj_mol": repo_diag._kjmol(float(mu_arr[PCSAFT_WATER] - a_term - z_term)),
            "simplex_fd_target_da_dx_kj_mol": repo_diag._kjmol(
                _forward_simplex_fd(x, PCSAFT_TARGET, PCSAFT_WATER, lambda x_new: float(term_state(x_new)[f"a_{suffix}"]))
            ),
            "simplex_fd_counter_da_dx_kj_mol": repo_diag._kjmol(
                _forward_simplex_fd(x, PCSAFT_COUNTER, PCSAFT_WATER, lambda x_new: float(term_state(x_new)[f"a_{suffix}"]))
            ),
            "simplex_fd_water_da_dx_kj_mol": repo_diag._kjmol(
                _forward_simplex_fd(x, PCSAFT_WATER, PCSAFT_TARGET, lambda x_new: float(term_state(x_new)[f"a_{suffix}"]))
            ),
            "pcsaft_exact_target_dadx_kj_mol": repo_diag._kjmol(float(np.asarray(terms[f"dadx_{suffix}"], dtype=float)[PCSAFT_TARGET])),
            "pcsaft_exact_sum_xj_dadx_kj_mol": repo_diag._kjmol(float(terms[f"sum_x_dadx_{suffix}"])),
            "pcsaft_exact_resid_a_plus_z_minus_sumxj_dadx_kj_mol": repo_diag._kjmol(float(a_term + z_term - float(terms[f"sum_x_dadx_{suffix}"]))),
        }
    )
    row.update(_pcsaft_exact_extra(contribution, state, geom))
    return row


def _feos_parameter_lookup(records: list[dict[str, object]], name: str) -> dict[str, object]:
    return next(record for record in records if record.get("identifier", {}).get("name") == name)


def _feos_base():
    state_inf, p_ref_pa = feos_extractor._build_reference_state(FEOS_EOS, FEOS_TARGET)
    x = np.asarray(state_inf.molefracs, dtype=float)
    density = float(state_inf.density / (si.MOL / si.METER**3))
    volume_m3 = float(state_inf.volume / (si.METER**3))
    pressure_pa = _to_pa(state_inf.pressure())
    z_total = pressure_pa * volume_m3 / (overlay.R_GAS * overlay.T_REF)
    z_residual = z_total - 1.0
    records = json.loads(feos_extractor.FEOS_PARAMETER_PATH.read_text(encoding="utf-8"))
    binary_entries = json.loads(feos_extractor.FEOS_BINARY_PATH.read_text(encoding="utf-8"))
    water = _feos_parameter_lookup(records, FEOS_COMPONENTS[FEOS_WATER])
    target = _feos_parameter_lookup(records, FEOS_COMPONENTS[FEOS_TARGET])
    counter = _feos_parameter_lookup(records, FEOS_COMPONENTS[FEOS_COUNTER])

    def _assoc_sites(record: dict[str, object]) -> int:
        sites = record.get("association_sites")
        if not sites:
            return 0
        total = 0
        for site in sites:
            total += int(round(float(site.get("na", 0.0)) + float(site.get("nb", 0.0))))
        return total

    def _assoc_energy(record: dict[str, object]) -> float:
        sites = record.get("association_sites") or []
        value = sites[0].get("epsilon_k_ab") if sites else None
        return float(value) if value is not None else 0.0

    def _assoc_volume(record: dict[str, object]) -> float:
        sites = record.get("association_sites") or []
        value = sites[0].get("kappa_ab") if sites else None
        return float(value) if value is not None else 0.0

    def _binary_value(name1: str, name2: str, key: str, default: float = 0.0) -> float:
        for entry in binary_entries:
            id1 = entry.get("id1", {}).get("name")
            id2 = entry.get("id2", {}).get("name")
            if {id1, id2} == {name1, name2}:
                value = entry.get("model_record", {}).get(key)
                return float(value) if value is not None else default
        return default

    common = _common_fields(
        package="feos",
        config="advanced_held2014",
        component_order="|".join(FEOS_COMPONENTS),
        target_component=FEOS_COMPONENTS[FEOS_TARGET],
        counter_component=FEOS_COMPONENTS[FEOS_COUNTER],
        water_component=FEOS_COMPONENTS[FEOS_WATER],
        x_target=float(x[FEOS_TARGET]),
        x_counter=float(x[FEOS_COUNTER]),
        x_water=float(x[FEOS_WATER]),
        density_mol_m3=density,
        compressibility_factor=z_total,
        sigma_target=float(target["sigma"]),
        sigma_counter=float(counter["sigma"]),
        sigma_water=float(water["sigma"]),
        epsilon_target=float(target["epsilon_k"]),
        epsilon_counter=float(counter["epsilon_k"]),
        epsilon_water=float(water["epsilon_k"]),
        born_radius_target=0.0,
        born_radius_counter=0.0,
        water_target_k=_binary_value(FEOS_COMPONENTS[FEOS_WATER], FEOS_COMPONENTS[FEOS_TARGET], "k_ij"),
        water_counter_k=_binary_value(FEOS_COMPONENTS[FEOS_WATER], FEOS_COMPONENTS[FEOS_COUNTER], "k_ij"),
        target_counter_k=_binary_value(FEOS_COMPONENTS[FEOS_TARGET], FEOS_COMPONENTS[FEOS_COUNTER], "k_ij"),
        dielectric_target=float((target.get("permittivity_record") or {}).get("ExperimentalData", {}).get("data", [[0.0, 0.0]])[0][1]),
        dielectric_counter=float((counter.get("permittivity_record") or {}).get("ExperimentalData", {}).get("data", [[0.0, 0.0]])[0][1]),
        dielectric_water=float((water.get("permittivity_record") or {}).get("ExperimentalData", {}).get("data", [[0.0, 0.0]])[1][1]),
        target_segment_m=float(target["m"]),
        counter_segment_m=float(counter["m"]),
        water_segment_m=float(water["m"]),
        target_assoc_sites=_assoc_sites(target),
        counter_assoc_sites=_assoc_sites(counter),
        water_assoc_sites=_assoc_sites(water),
        target_assoc_energy_k=_assoc_energy(target),
        counter_assoc_energy_k=_assoc_energy(counter),
        water_assoc_energy_k=_assoc_energy(water),
        target_assoc_volume=_assoc_volume(target),
        counter_assoc_volume=_assoc_volume(counter),
        water_assoc_volume=_assoc_volume(water),
        reference_pressure_pa=float(p_ref_pa),
        state_pressure_pa=pressure_pa,
        state_volume_m3=volume_m3,
    )

    def _a_contrib(x_new: np.ndarray) -> dict[str, float]:
        state = feos.State(
            FEOS_EOS,
            temperature=overlay.T_REF * si.KELVIN,
            density=state_inf.density,
            molefracs=x_new,
            total_moles=si.MOL,
        )
        pairs = state.residual_molar_helmholtz_energy_contributions()
        out: dict[str, float] = {}
        for label, quantity in pairs:
            out[_label_slug(label)] = _to_kj_per_mol(quantity)
        return out

    return state_inf, x, common, z_total, z_residual, _a_contrib


def _feos_raw_label_dict(pairs, converter) -> dict[str, float]:
    return {_label_slug(label): converter(quantity) for label, quantity in pairs}


def _feos_row(
    contribution: str,
    state_inf,
    x: np.ndarray,
    common: dict[str, object],
    z_total: float,
    z_residual: float,
    a_eval,
) -> dict[str, object]:
    a_raw = _feos_raw_label_dict(state_inf.residual_molar_helmholtz_energy_contributions(), _to_kj_per_mol)
    p_raw = _feos_raw_label_dict(state_inf.pressure_contributions(), _to_pa)
    mu_target_raw = _feos_raw_label_dict(
        state_inf.chemical_potential_contributions(FEOS_TARGET, feos.Contributions.Residual),
        _to_kj_per_mol,
    )
    mu_counter_raw = _feos_raw_label_dict(
        state_inf.chemical_potential_contributions(FEOS_COUNTER, feos.Contributions.Residual),
        _to_kj_per_mol,
    )
    mu_water_raw = _feos_raw_label_dict(
        state_inf.chemical_potential_contributions(FEOS_WATER, feos.Contributions.Residual),
        _to_kj_per_mol,
    )

    if contribution == "hc":
        a_term = a_raw["hard_sphere"] + a_raw["hard_chain"]
        p_alpha = p_raw["hard_sphere"] + p_raw["hard_chain"]
        target_mu = mu_target_raw["hard_sphere"] + mu_target_raw["hard_chain"]
        counter_mu = mu_counter_raw["hard_sphere"] + mu_counter_raw["hard_chain"]
        water_mu = mu_water_raw["hard_sphere"] + mu_water_raw["hard_chain"]
        extra = {
            "feos_raw_a_hard_sphere_kj_mol": a_raw["hard_sphere"],
            "feos_raw_a_hard_chain_kj_mol": a_raw["hard_chain"],
            "feos_raw_z_hard_sphere_kj_mol": overlay.R_GAS * overlay.T_REF * p_raw["hard_sphere"] * float(common["state_volume_m3"]) / (overlay.R_GAS * overlay.T_REF) / 1000.0,
            "feos_raw_z_hard_chain_kj_mol": overlay.R_GAS * overlay.T_REF * p_raw["hard_chain"] * float(common["state_volume_m3"]) / (overlay.R_GAS * overlay.T_REF) / 1000.0,
            "feos_raw_target_mu_hard_sphere_kj_mol": mu_target_raw["hard_sphere"],
            "feos_raw_target_mu_hard_chain_kj_mol": mu_target_raw["hard_chain"],
            "feos_raw_counter_mu_hard_sphere_kj_mol": mu_counter_raw["hard_sphere"],
            "feos_raw_counter_mu_hard_chain_kj_mol": mu_counter_raw["hard_chain"],
            "feos_raw_water_mu_hard_sphere_kj_mol": mu_water_raw["hard_sphere"],
            "feos_raw_water_mu_hard_chain_kj_mol": mu_water_raw["hard_chain"],
        }
        fd_keys = ("hard_sphere", "hard_chain")
    else:
        key = {
            "disp": "dispersion",
            "assoc": "association",
            "dh": "ionic",
            "born": "born",
        }[contribution]
        a_term = a_raw[key]
        p_alpha = p_raw[key]
        target_mu = mu_target_raw[key]
        counter_mu = mu_counter_raw[key]
        water_mu = mu_water_raw[key]
        extra = {
            f"feos_raw_a_{key}_kj_mol": a_raw[key],
            f"feos_raw_target_mu_{key}_kj_mol": target_mu,
            f"feos_raw_counter_mu_{key}_kj_mol": counter_mu,
            f"feos_raw_water_mu_{key}_kj_mol": water_mu,
        }
        fd_keys = (key,)

    volume_m3 = float(common["state_volume_m3"])
    z_alpha = p_alpha * volume_m3 / (overlay.R_GAS * overlay.T_REF)
    z_term = overlay.R_GAS * overlay.T_REF * z_alpha / 1000.0
    z_correction = overlay.R_GAS * overlay.T_REF * (-(z_alpha / z_residual) * math.log(z_total)) / 1000.0

    def _fd_eval(x_new: np.ndarray) -> float:
        values = a_eval(x_new)
        return float(sum(values[k] for k in fd_keys))

    row = dict(common)
    row.update(
        {
            "contribution": contribution,
            "a_kj_mol": float(a_term),
            "z_kj_mol": float(z_term),
            "a_plus_z_kj_mol": float(a_term + z_term),
            "target_mu_kj_mol": float(target_mu),
            "counter_mu_kj_mol": float(counter_mu),
            "water_mu_kj_mol": float(water_mu),
            "target_lnfug_kj_mol": float(target_mu + z_correction),
            "target_z_correction_kj_mol": float(z_correction),
            "weighted_mu_kj_mol": float(np.dot(x, np.asarray([water_mu, target_mu, counter_mu], dtype=float))),
            "weighted_mu_minus_a_plus_z_kj_mol": float(np.dot(x, np.asarray([water_mu, target_mu, counter_mu], dtype=float)) - (a_term + z_term)),
            "target_mu_minus_a_z_kj_mol": float(target_mu - a_term - z_term),
            "counter_mu_minus_a_z_kj_mol": float(counter_mu - a_term - z_term),
            "water_mu_minus_a_z_kj_mol": float(water_mu - a_term - z_term),
            "simplex_fd_target_da_dx_kj_mol": _forward_simplex_fd(x, FEOS_TARGET, FEOS_WATER, _fd_eval),
            "simplex_fd_counter_da_dx_kj_mol": _forward_simplex_fd(x, FEOS_COUNTER, FEOS_WATER, _fd_eval),
            "simplex_fd_water_da_dx_kj_mol": _forward_simplex_fd(x, FEOS_WATER, FEOS_TARGET, _fd_eval),
            "raw_label_count": len(a_raw),
        }
    )
    row.update(extra)
    return row


def _clapeyron_row(contribution: str, payload: dict) -> dict[str, object]:
    comp = payload["contributions"][contribution]
    x = np.asarray(payload["z"], dtype=float)
    mu = np.asarray(comp["mu_components_kj_mol"], dtype=float)
    mumz = np.asarray(comp["mu_minus_a_z_components_kj_mol"], dtype=float)
    k_matrix = np.asarray(payload["k_matrix"], dtype=float)
    sigma_diag = np.asarray(payload["sigma_diag_m"], dtype=float) * 1.0e10
    sigma_born = np.asarray(payload["sigma_born_m"], dtype=float) * 1.0e10
    eps_diag = np.asarray(payload["epsilon_diag_k"], dtype=float)
    diel = np.asarray(payload["dielectric_loaded"], dtype=float)
    row = _common_fields(
        package="Clapeyron.jl",
        config=payload["config"],
        component_order="|".join(payload["components"]),
        target_component=payload["target_component"],
        counter_component=payload["counter_component"],
        water_component=payload["components"][0],
        x_target=float(x[1]),
        x_counter=float(x[2]),
        x_water=float(x[0]),
        density_mol_m3=float(payload["density_mol_m3"]),
        compressibility_factor=float(payload["compressibility_factor"]),
        sigma_target=float(sigma_diag[1]),
        sigma_counter=float(sigma_diag[2]),
        sigma_water=float(sigma_diag[0]),
        epsilon_target=float(eps_diag[1]),
        epsilon_counter=float(eps_diag[2]),
        epsilon_water=float(eps_diag[0]),
        born_radius_target=float(sigma_born[1]),
        born_radius_counter=float(sigma_born[2]),
        water_target_k=float(k_matrix[0, 1]),
        water_counter_k=float(k_matrix[0, 2]),
        target_counter_k=float(k_matrix[1, 2]),
        dielectric_target=float(diel[1]),
        dielectric_counter=float(diel[2]),
        dielectric_water=float(diel[0]),
        state_volume_m3=float(payload["state_volume_m3"]),
        state_pressure_pa=float(payload["reference_pressure_pa"]),
        reference_pressure_pa=float(payload["reference_pressure_pa"]),
        target_segment_m=float(payload.get("segment_diag", [None, None, None])[1]) if payload.get("segment_diag") else None,
        counter_segment_m=float(payload.get("segment_diag", [None, None, None])[2]) if payload.get("segment_diag") else None,
        water_segment_m=float(payload.get("segment_diag", [None, None, None])[0]) if payload.get("segment_diag") else None,
    )
    row.update(
        {
            "contribution": contribution,
            "a_kj_mol": float(comp["a_kj_mol"]),
            "z_kj_mol": float(comp["z_kj_mol"]),
            "a_plus_z_kj_mol": float(comp["a_plus_z_kj_mol"]),
            "target_mu_kj_mol": float(mu[1]),
            "counter_mu_kj_mol": float(mu[2]),
            "water_mu_kj_mol": float(mu[0]),
            "target_lnfug_kj_mol": float(comp["lnfug_target_kj_mol"]),
            "target_z_correction_kj_mol": float(comp["z_correction_kj_mol"]),
            "weighted_mu_kj_mol": float(comp["weighted_mu_kj_mol"]),
            "weighted_mu_minus_a_plus_z_kj_mol": float(comp["weighted_mu_kj_mol"] - comp["a_plus_z_kj_mol"]),
            "target_mu_minus_a_z_kj_mol": float(mumz[1]),
            "counter_mu_minus_a_z_kj_mol": float(mumz[2]),
            "water_mu_minus_a_z_kj_mol": float(mumz[0]),
            "simplex_fd_target_da_dx_kj_mol": _float_or_nan(comp["simplex_fd_target_kj_mol"]),
            "simplex_fd_counter_da_dx_kj_mol": _float_or_nan(comp["simplex_fd_counter_kj_mol"]),
            "simplex_fd_water_da_dx_kj_mol": _float_or_nan(comp["simplex_fd_water_kj_mol"]),
        }
    )
    if contribution == "assoc":
        x_vals = comp.get("assoc_X_values", [])
        delta_vals = comp.get("assoc_delta_values", [])
        site_matrix = comp.get("assoc_site_matrix", [])
        row["clapeyron_assoc_pair_count"] = len(delta_vals)
        for idx, value in enumerate(x_vals, start=1):
            row[f"clapeyron_assoc_X_{idx}"] = float(value)
        for idx, value in enumerate(delta_vals, start=1):
            row[f"clapeyron_assoc_delta_{idx}"] = float(value)
        for i, matrix_row in enumerate(site_matrix, start=1):
            for j, value in enumerate(matrix_row, start=1):
                row[f"clapeyron_assoc_site_matrix_{i}_{j}"] = float(value)
    return row


def _summary_from_row(row: dict[str, object]) -> dict[str, object]:
    excluded = {
        "contribution",
        "a_kj_mol",
        "z_kj_mol",
        "a_plus_z_kj_mol",
        "target_mu_kj_mol",
        "counter_mu_kj_mol",
        "water_mu_kj_mol",
        "target_lnfug_kj_mol",
        "target_z_correction_kj_mol",
        "weighted_mu_kj_mol",
        "weighted_mu_minus_a_plus_z_kj_mol",
        "target_mu_minus_a_z_kj_mol",
        "counter_mu_minus_a_z_kj_mol",
        "water_mu_minus_a_z_kj_mol",
        "simplex_fd_target_da_dx_kj_mol",
        "simplex_fd_counter_da_dx_kj_mol",
        "simplex_fd_water_da_dx_kj_mol",
    }
    return {key: value for key, value in row.items() if key not in excluded}


def _write_report(rows_by_contribution: dict[str, list[dict[str, object]]]) -> None:
    lines = [
        "# Na+ Package Component Breakdown Notes",
        "",
        "These CSV files compare the current `PC-SAFT` package, `feos`, and `Clapeyron.jl` for the Na+/Cl-/water infinite-dilution state used in Figure 3 package comparisons.",
        "",
        "## Column meanings",
        "",
        "- `a_kj_mol`: contribution Helmholtz term $RT a^\\alpha$.",
        "- `z_kj_mol`: contribution compressibility term $RT Z^\\alpha$.",
        "- `target_mu_kj_mol`: package-exposed Na+ contribution on the $\\mu^\\alpha$ basis when available.",
        "- `target_lnfug_kj_mol`: Na+ fugacity-style contribution reconstructed as $RT\\left(\\mu^\\alpha - \\frac{Z^\\alpha}{Z-1}\\ln Z\\right)$.",
        "- `simplex_fd_*`: finite-difference proxy for $RT\\,\\partial a^\\alpha/\\partial x_k$ at fixed $T,\\rho$ or fixed $T,V$ depending on the package API.",
        "",
        "## Exposure differences",
        "",
        "- `pcsaft` exposes explicit analytical derivative pieces for `hc`, `disp`, and `assoc`, including the same structures documented in `equations_v2.tex` such as $\\partial a^{hs}/\\partial x_k$, $\\partial g^{hs}/\\partial x_k$, and the association-site derivatives solved from the linear system.",
        "- `feos` exposes contribution-resolved residual Helmholtz energies, pressures, and chemical-potential contributions, but not the sub-derivatives behind hard-sphere, hard-chain, or dispersion. The CSV therefore records raw contribution labels plus finite-difference proxies.",
        "- `Clapeyron.jl` exposes the per-contribution Helmholtz functions and obtains $\\mu$-style quantities through `VT_molar_gradient`. It exposes real association site fractions $X$, association strengths $\\Delta$, and the site matrix, but not named derivative subpieces like $dX/dx_k$ or $d\\Delta/dx_k` through a public API.",
        "",
        "## Files",
        "",
    ]
    for contribution, rows in rows_by_contribution.items():
        packages = ", ".join(str(row["package"]) for row in rows)
        lines.append(f"- `{contribution}`: rows for {packages}.")
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    pc_state, pc_geom, pc_terms, pc_x, pc_common, pc_term_state = _pcsaft_base()
    feos_state, feos_x, feos_common, feos_z_total, feos_z_residual, feos_a_eval = _feos_base()
    clapeyron_payload = _run_clapeyron_json()

    rows_by_contribution: dict[str, list[dict[str, object]]] = {key: [] for key in OUTPUTS}
    summary_rows: list[dict[str, object]] = []

    for contribution in OUTPUTS:
        pc_row = _pcsaft_row(contribution, pc_state, pc_geom, pc_terms, pc_x, pc_common, pc_term_state)
        feos_row = _feos_row(contribution, feos_state, feos_x, feos_common, feos_z_total, feos_z_residual, feos_a_eval)
        clapeyron_row = _clapeyron_row(contribution, clapeyron_payload)
        rows = [pc_row, feos_row, clapeyron_row]
        rows_by_contribution[contribution] = rows
        _write_csv(OUTPUTS[contribution], rows)
        summary_rows.extend(_summary_from_row(row) for row in rows if contribution == "hc")

    _write_csv(OUTPUT_COMMON, summary_rows)
    _write_report(rows_by_contribution)
    print(f"Wrote {OUTPUT_COMMON}", flush=True)
    for output in OUTPUTS.values():
        print(f"Wrote {output}", flush=True)
    print(f"Wrote {REPORT_PATH}", flush=True)


if __name__ == "__main__":
    main()
