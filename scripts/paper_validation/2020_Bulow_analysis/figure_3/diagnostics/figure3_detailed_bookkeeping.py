from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
FIGURE_DIR = SCRIPT_DIR.parent
ANALYSIS_ROOT = FIGURE_DIR.parent
REPO_ROOT = ANALYSIS_ROOT.parents[2]

if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts._env import require_epcsaft_install

require_epcsaft_install()

import _model_overlay as overlay
import _plot_common as common
from epcsaft.parameters import get_prop_dict
from scripts._epcsaft_oop import epcsaft_density, epcsaft_fugacity_coefficient_terms, epcsaft_pressure


DATA_PATH = FIGURE_DIR / "data" / "water_contributions.csv"
OUTPUT_BOOKKEEPING = SCRIPT_DIR / "figure3_detailed_bookkeeping.csv"
OUTPUT_HC = SCRIPT_DIR / "figure3_hc_dadx_components.csv"
OUTPUT_DISP = SCRIPT_DIR / "figure3_disp_dadx_components.csv"
OUTPUT_ASSOC = SCRIPT_DIR / "figure3_assoc_dadx_components.csv"
OUTPUT_DH = SCRIPT_DIR / "figure3_dh_dadx_components.csv"
OUTPUT_BORN = SCRIPT_DIR / "figure3_born_dadx_components.csv"

T_REF = overlay.T_REF
P_REF = overlay.P_REF
EPS = overlay.EPS
EPS_INF = overlay.EPS_INF
R_GAS = overlay.R_GAS
PI = math.pi
N_AV = 6.02214076e23
RT_KJMOL = R_GAS * T_REF / 1000.0

CONTRIBUTION_MAP = {
    "hc": {"paper_rows": ("hc avg", "hc"), "suffix": "hc"},
    "disp": {"paper_rows": ("disp avg", "disp"), "suffix": "disp"},
    "assoc": {"paper_rows": ("assoc avg", "assoc"), "suffix": "assoc"},
    "dh": {"paper_rows": (), "suffix": "ion"},
    "born": {"paper_rows": ("born avg", "born"), "suffix": "born"},
}

A0 = np.asarray([0.9105631445, 0.6361281449, 2.6861347891, -26.547362491, 97.759208784, -159.59154087, 91.297774084], dtype=float)
A1 = np.asarray([-0.3084016918, 0.1860531159, -2.5030047259, 21.419793629, -65.255885330, 83.318680481, -33.746922930], dtype=float)
A2 = np.asarray([-0.0906148351, 0.4527842806, 0.5962700728, -1.7241829131, -4.1302112531, 13.776631870, -8.6728470368], dtype=float)
B0 = np.asarray([0.7240946941, 2.2382791861, -4.0025849485, -21.003576815, 26.855641363, 206.55133841, -355.60235612], dtype=float)
B1 = np.asarray([-0.5755498075, 0.6995095521, 3.8925673390, -17.215471648, 192.67226447, -161.82646165, -165.20769346], dtype=float)
B2 = np.asarray([0.0976883116, -0.2557574982, -9.1558561530, 20.642075974, -38.804430052, 93.626774077, -29.666905585], dtype=float)


def _kjmol(value: float) -> float:
    return float(RT_KJMOL * float(value))


def _stable_logz_over_zminus1(z_total: float) -> float:
    if abs(z_total - 1.0) < 1.0e-10:
        return 1.0
    return math.log(z_total) / (z_total - 1.0)


def _paper_value(frame: common.Table, contribution: str, ion: str) -> float:
    if contribution == "dh":
        return 0.0
    for row_key in CONTRIBUTION_MAP[contribution]["paper_rows"]:
        if row_key in frame.index:
            return float(frame.scalar(row_key, ion))
    raise KeyError(f"Missing paper row for contribution '{contribution}'.")


def _vector_param(params: dict[str, object], key: str, ncomp: int, default: float = 0.0) -> np.ndarray:
    value = params.get(key)
    if value is None:
        return np.full(ncomp, default, dtype=float)
    arr = np.asarray(value, dtype=float)
    if arr.size == 0:
        return np.full(ncomp, default, dtype=float)
    return arr.astype(float)


def _matrix_param(params: dict[str, object], key: str, ncomp: int, default: float = 0.0) -> np.ndarray:
    value = params.get(key)
    if value is None:
        return np.full((ncomp, ncomp), default, dtype=float)
    arr = np.asarray(value, dtype=float)
    if arr.size == 0:
        return np.full((ncomp, ncomp), default, dtype=float)
    return arr.reshape((ncomp, ncomp)).astype(float)


def _effective_diameters(params: dict[str, object], t: float) -> np.ndarray:
    s = np.asarray(params["s"], dtype=float)
    e = np.asarray(params["e"], dtype=float)
    z = np.asarray(params.get("z", np.zeros_like(s)), dtype=float)
    d_ion_mode = int(params.get("d_ion_mode", 1))
    d = s * (1.0 - 0.12 * np.exp(-3.0 * e / t))
    for i in range(len(s)):
        if abs(z[i]) <= 1.0e-12:
            continue
        if d_ion_mode == 0:
            d[i] = s[i]
        elif d_ion_mode == 1:
            d[i] = s[i] * (1.0 - 0.12)
        elif d_ion_mode == 2:
            d[i] = s[i] * (1.0 - 0.12 * math.exp(-3.0 * e[i] / t))
        else:
            raise ValueError(f"Unsupported d_ion_mode '{d_ion_mode}'.")
    return d


def _state_for_ion(ion: str) -> dict[str, object]:
    species = overlay._species_for_ion(ion, "water")
    x = np.asarray([EPS, EPS, 1.0 - 2.0 * EPS], dtype=float)
    params = get_prop_dict("2020_Bulow", species, x, T_REF)
    rho = epcsaft_density(T_REF, P_REF, x, params, phase="liq")
    z = np.asarray(params.get("z", []), dtype=float)
    idx_ion = np.where(np.abs(z) > 1.0e-12)[0]
    idx_solv = np.where(np.abs(z) <= 1.0e-12)[0]
    x_ref = x.copy()
    x_ref[idx_ion] = 0.0
    solv_sum = float(np.sum(x_ref[idx_solv]))
    if solv_sum > 0.0:
        x_ref[idx_solv] /= solv_sum
    else:
        x_ref[idx_solv] = 1.0 / len(idx_solv)
    p_ref = epcsaft_pressure(T_REF, rho, x_ref, params)
    x_inf = x_ref.copy()
    ion_idx = species.index(ion)
    x_inf[ion_idx] = EPS_INF
    x_inf /= np.sum(x_inf)
    phase = "vap" if rho < 900.0 else "liq"
    rho_inf = epcsaft_density(T_REF, p_ref, x_inf, params, phase=phase)
    terms = epcsaft_fugacity_coefficient_terms(T_REF, rho_inf, x_inf, params)
    return {
        "ion": ion,
        "species": species,
        "x_inf": x_inf,
        "rho_inf": float(rho_inf),
        "params": params,
        "terms": terms,
        "ion_idx": ion_idx,
    }


def _mixture_geometry(state: dict[str, object]) -> dict[str, object]:
    params = state["params"]
    x = np.asarray(state["x_inf"], dtype=float)
    rho = float(state["rho_inf"])
    ncomp = len(x)
    m = np.asarray(params["m"], dtype=float)
    s = np.asarray(params["s"], dtype=float)
    e = np.asarray(params["e"], dtype=float)
    z = np.asarray(params.get("z", np.zeros(ncomp)), dtype=float)
    k_ij = _matrix_param(params, "k_ij", ncomp, default=0.0)
    l_ij = _matrix_param(params, "l_ij", ncomp, default=0.0)
    d = _effective_diameters(params, T_REF)
    den = rho * N_AV / 1.0e30
    zeta = np.zeros(4, dtype=float)
    for order in range(4):
        zeta[order] = PI / 6.0 * den * np.sum(x * m * (d ** order))
    eta = float(zeta[3])
    m_avg = float(np.sum(x * m))
    s_ij = np.zeros((ncomp, ncomp), dtype=float)
    e_ij = np.zeros((ncomp, ncomp), dtype=float)
    ghs = np.zeros((ncomp, ncomp), dtype=float)
    m2es3 = 0.0
    m2e2s3 = 0.0
    for i in range(ncomp):
        for j in range(ncomp):
            s_ij[i, j] = 0.5 * (s[i] + s[j]) * (1.0 - l_ij[i, j])
            if np.any(np.abs(z) > 1.0e-12):
                if z[i] * z[j] <= 0.0:
                    e_ij[i, j] = math.sqrt(e[i] * e[j]) * (1.0 - k_ij[i, j])
            else:
                e_ij[i, j] = math.sqrt(e[i] * e[j]) * (1.0 - k_ij[i, j])
            m2es3 += x[i] * x[j] * m[i] * m[j] * (e_ij[i, j] / T_REF) * (s_ij[i, j] ** 3)
            m2e2s3 += x[i] * x[j] * m[i] * m[j] * ((e_ij[i, j] / T_REF) ** 2) * (s_ij[i, j] ** 3)
            dij = d[i] * d[j] / (d[i] + d[j])
            ghs[i, j] = 1.0 / (1.0 - zeta[3]) + dij * 3.0 * zeta[2] / (1.0 - zeta[3]) ** 2 + (dij ** 2) * 2.0 * (zeta[2] ** 2) / (1.0 - zeta[3]) ** 3
    ares_hs = 1.0 / zeta[0] * (3.0 * zeta[1] * zeta[2] / (1.0 - zeta[3]) + (zeta[2] ** 3) / (zeta[3] * (1.0 - zeta[3]) ** 2) + (((zeta[2] ** 3) / (zeta[3] ** 2)) - zeta[0]) * math.log(1.0 - zeta[3]))
    zhs = zeta[3] / (1.0 - zeta[3]) + 3.0 * zeta[1] * zeta[2] / zeta[0] / (1.0 - zeta[3]) ** 2 + (3.0 * (zeta[2] ** 3) - zeta[3] * (zeta[2] ** 3)) / zeta[0] / (1.0 - zeta[3]) ** 3
    a = A0 + (m_avg - 1.0) / m_avg * A1 + (m_avg - 1.0) / m_avg * (m_avg - 2.0) / m_avg * A2
    b = B0 + (m_avg - 1.0) / m_avg * B1 + (m_avg - 1.0) / m_avg * (m_avg - 2.0) / m_avg * B2
    i1 = float(np.sum(a * (eta ** np.arange(7))))
    i2 = float(np.sum(b * (eta ** np.arange(7))))
    c1 = 1.0 / (1.0 + m_avg * (8.0 * eta - 2.0 * eta * eta) / (1.0 - eta) ** 4 + (1.0 - m_avg) * (20.0 * eta - 27.0 * eta * eta + 12.0 * eta**3 - 2.0 * eta**4) / ((1.0 - eta) * (2.0 - eta)) ** 2)
    c2 = -(c1**2) * (m_avg * (-4.0 * eta * eta + 20.0 * eta + 8.0) / (1.0 - eta) ** 5 + (1.0 - m_avg) * (2.0 * eta**3 + 12.0 * eta * eta - 48.0 * eta + 40.0) / ((1.0 - eta) * (2.0 - eta)) ** 3)
    return {
        "x": x, "rho": rho, "den": den, "m": m, "s": s, "e": e, "z": z, "d": d, "zeta": zeta, "eta": eta, "m_avg": m_avg,
        "s_ij": s_ij, "e_ij": e_ij, "ghs": ghs, "m2es3": m2es3, "m2e2s3": m2e2s3, "ares_hs": float(ares_hs),
        "zhs": float(zhs), "a_coeffs": a, "b_coeffs": b, "i1": i1, "i2": i2, "c1": float(c1), "c2": float(c2), "params": params,
    }


def _generic_contribution_rows(state: dict[str, object], frame: common.Table, contribution: str) -> list[dict[str, object]]:
    terms = state["terms"]
    z_weight = _stable_logz_over_zminus1(float(terms["z_total"]))
    rows: list[dict[str, object]] = []
    suffix = CONTRIBUTION_MAP[contribution]["suffix"]
    zvals = np.asarray(state["params"]["z"], dtype=float)
    for idx, ion in enumerate(state["species"]):
        paper_mu = 0.0 if contribution == "dh" else (_paper_value(frame, contribution, ion) if abs(zvals[idx]) > 1.0e-12 else 0.0)
        mu_red = float(np.asarray(terms[f"mu_{suffix}"], dtype=float)[idx])
        lnfug_red = float(np.asarray(terms[f"lnfugcoef_{suffix}"], dtype=float)[idx])
        a_term = float(terms[f"a_{suffix}"])
        z_term = float(terms[f"z_{suffix}"])
        dadx_term = float(np.asarray(terms[f"dadx_{suffix}"], dtype=float)[idx])
        sum_x_dadx_term = float(terms[f"sum_x_dadx_{suffix}"])
        z_quot_red = -z_term * z_weight
        rows.append({
            "ion": ion, "contr": contribution, "paper_mu": paper_mu, "epcsaft_mu": _kjmol(mu_red),
            "epcsaft_mu_manual": _kjmol(a_term + z_term + dadx_term - sum_x_dadx_term), "a": _kjmol(a_term),
            "z": _kjmol(z_term), "dadx": _kjmol(dadx_term), "sum_xj_dadx": _kjmol(sum_x_dadx_term),
            "z_quot": _kjmol(z_quot_red), "epcsaft_lnfug": _kjmol(lnfug_red), "epcsaft_lnfug_manual": _kjmol(mu_red + z_quot_red),
            "lnfug_red": lnfug_red, "lnact_star_red": 0.0, "gamma_star": 1.0, "Ghyd": _kjmol(lnfug_red),
            "manual_minus_epcsaft_mu": _kjmol(a_term + z_term + dadx_term - sum_x_dadx_term - mu_red),
            "manual_minus_epcsaft_lnfug": _kjmol((mu_red + z_quot_red) - lnfug_red), "paper_minus_epcsaft_mu": paper_mu - _kjmol(mu_red),
        })
    return rows


def _total_row(ion: str, rows: list[dict[str, object]], total_lnfug_red: float) -> dict[str, object]:
    sum_gamma_logs = float(np.sum([float(row["lnact_star_red"]) for row in rows]))
    return {
        "ion": ion, "contr": "total",
        "paper_mu": float(np.sum([float(row["paper_mu"]) for row in rows])),
        "epcsaft_mu": float(np.sum([float(row["epcsaft_mu"]) for row in rows])),
        "epcsaft_mu_manual": float(np.sum([float(row["epcsaft_mu_manual"]) for row in rows])),
        "a": float(np.sum([float(row["a"]) for row in rows])),
        "z": float(np.sum([float(row["z"]) for row in rows])),
        "dadx": float(np.sum([float(row["dadx"]) for row in rows])),
        "sum_xj_dadx": float(np.sum([float(row["sum_xj_dadx"]) for row in rows])),
        "z_quot": float(np.sum([float(row["z_quot"]) for row in rows])),
        "epcsaft_lnfug": float(np.sum([float(row["epcsaft_lnfug"]) for row in rows])),
        "epcsaft_lnfug_manual": float(np.sum([float(row["epcsaft_lnfug_manual"]) for row in rows])),
        "lnfug_red": float(np.sum([float(row["lnfug_red"]) for row in rows])),
        "lnact_star_red": sum_gamma_logs,
        "gamma_star": math.exp(sum_gamma_logs),
        "Ghyd": float(np.sum([float(row["Ghyd"]) for row in rows])),
        "manual_minus_epcsaft_mu": float(np.sum([float(row["manual_minus_epcsaft_mu"]) for row in rows])),
        "manual_minus_epcsaft_lnfug": float(np.sum([float(row["manual_minus_epcsaft_lnfug"]) for row in rows])),
        "paper_minus_epcsaft_mu": float(np.sum([float(row["paper_minus_epcsaft_mu"]) for row in rows])),
        "epcsaft_total_lnfug_from_api": _kjmol(total_lnfug_red),
    }


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _hc_breakdown(state: dict[str, object], geom: dict[str, object], frame: common.Table) -> list[dict[str, object]]:
    x = np.asarray(geom["x"], dtype=float)
    m = np.asarray(geom["m"], dtype=float)
    d = np.asarray(geom["d"], dtype=float)
    zeta = np.asarray(geom["zeta"], dtype=float)
    ghs = np.asarray(geom["ghs"], dtype=float)
    ares_hs = float(geom["ares_hs"])
    m_avg = float(geom["m_avg"])
    den = float(geom["den"])
    terms = state["terms"]
    ions = state["species"]
    ncomp = len(x)
    dghsii_dx = np.zeros((ncomp, ncomp), dtype=float)
    dahs_dx = np.zeros(ncomp, dtype=float)
    dzeta_dx = np.zeros(4, dtype=float)
    for i in range(ncomp):
        for order in range(4):
            dzeta_dx[order] = PI / 6.0 * den * m[i] * (d[i] ** order)
        for j in range(ncomp):
            djj = d[j] * d[j] / (d[j] + d[j])
            dghsii_dx[i, j] = (
                dzeta_dx[3] / (1.0 - zeta[3]) ** 2
                + djj * (3.0 * dzeta_dx[2] / (1.0 - zeta[3]) ** 2 + 6.0 * zeta[2] * dzeta_dx[3] / (1.0 - zeta[3]) ** 3)
                + (djj**2) * (4.0 * zeta[2] * dzeta_dx[2] / (1.0 - zeta[3]) ** 3 + 6.0 * zeta[2] * zeta[2] * dzeta_dx[3] / (1.0 - zeta[3]) ** 4)
            )
        dahs_dx[i] = (
            -dzeta_dx[0] / zeta[0] * ares_hs
            + 1.0 / zeta[0]
            * (
                3.0 * (dzeta_dx[1] * zeta[2] + zeta[1] * dzeta_dx[2]) / (1.0 - zeta[3])
                + 3.0 * zeta[1] * zeta[2] * dzeta_dx[3] / (1.0 - zeta[3]) ** 2
                + 3.0 * zeta[2] * zeta[2] * dzeta_dx[2] / zeta[3] / (1.0 - zeta[3]) ** 2
                + (zeta[2] ** 3) * dzeta_dx[3] * (3.0 * zeta[3] - 1.0) / (zeta[3] ** 2) / (1.0 - zeta[3]) ** 3
                + math.log(1.0 - zeta[3]) * (((3.0 * zeta[2] * zeta[2] * dzeta_dx[2] * zeta[3]) - 2.0 * (zeta[2] ** 3) * dzeta_dx[3]) / (zeta[3] ** 3) - dzeta_dx[0])
                + (zeta[0] - (zeta[2] ** 3) / (zeta[3] ** 2)) * dzeta_dx[3] / (1.0 - zeta[3])
            )
        )
    rows: list[dict[str, object]] = []
    zvals = np.asarray(state["params"]["z"], dtype=float)
    for i, ion in enumerate(ions):
        if abs(zvals[i]) <= 1.0e-12:
            continue
        term_mi_ares_hs = m[i] * ares_hs
        term_mavg_dahs = m_avg * dahs_dx[i]
        term_chain_sum = sum(x[j] * (m[j] - 1.0) / ghs[j, j] * dghsii_dx[i, j] for j in range(ncomp))
        term_log = (m[i] - 1.0) * math.log(ghs[i, i])
        dadx = term_mi_ares_hs + term_mavg_dahs - term_chain_sum - term_log
        mu = float(np.asarray(terms["mu_hc"], dtype=float)[i])
        paper_mu = _paper_value(frame, "hc", ion)
        rows.append({
            "ion": ion, "paper_mu": paper_mu, "epcsaft_mu": _kjmol(mu), "epcsaft_composition_derivative_residual_helmholtz": _kjmol(dadx),
            "resid_a_plus_z_minus_sumxj_dadx": _kjmol(float(terms["a_hc"]) + float(terms["z_hc"]) - float(terms["sum_x_dadx_hc"])),
            "m_i_ares_hs": _kjmol(term_mi_ares_hs), "m_avg_dahs_dx": _kjmol(term_mavg_dahs), "minus_chain_sum": _kjmol(-term_chain_sum),
            "minus_log_term": _kjmol(-term_log), "dadx_manual": _kjmol(dadx),
            "manual_minus_api_dadx": _kjmol(dadx - float(np.asarray(terms["dadx_hc"], dtype=float)[i])),
            "paper_minus_dadx": paper_mu - _kjmol(dadx), "paper_minus_mu": paper_mu - _kjmol(mu),
        })
    return rows


def _disp_breakdown(state: dict[str, object], geom: dict[str, object], frame: common.Table) -> list[dict[str, object]]:
    x = np.asarray(geom["x"], dtype=float)
    m = np.asarray(geom["m"], dtype=float)
    d = np.asarray(geom["d"], dtype=float)
    eta = float(geom["eta"])
    m_avg = float(geom["m_avg"])
    a = np.asarray(geom["a_coeffs"], dtype=float)
    b = np.asarray(geom["b_coeffs"], dtype=float)
    i1 = float(geom["i1"])
    i2 = float(geom["i2"])
    c1 = float(geom["c1"])
    m2es3 = float(geom["m2es3"])
    m2e2s3 = float(geom["m2e2s3"])
    s_ij = np.asarray(geom["s_ij"], dtype=float)
    e_ij = np.asarray(geom["e_ij"], dtype=float)
    den = float(geom["den"])
    terms = state["terms"]
    ions = state["species"]
    zvals = np.asarray(state["params"]["z"], dtype=float)
    rows: list[dict[str, object]] = []
    for i, ion in enumerate(ions):
        if abs(zvals[i]) <= 1.0e-12:
            continue
        dzeta3_dx = PI / 6.0 * den * m[i] * (d[i] ** 3)
        dI1_dx = 0.0
        dI2_dx = 0.0
        dm2es3_dx = 0.0
        dm2e2s3_dx = 0.0
        for l in range(7):
            daa_dx = m[i] / (m_avg**2) * A1[l] + m[i] / (m_avg**2) * (3.0 - 4.0 / m_avg) * A2[l]
            db_dx = m[i] / (m_avg**2) * B1[l] + m[i] / (m_avg**2) * (3.0 - 4.0 / m_avg) * B2[l]
            dI1_dx += a[l] * l * dzeta3_dx * (eta ** (l - 1)) + daa_dx * (eta**l)
            dI2_dx += b[l] * l * dzeta3_dx * (eta ** (l - 1)) + db_dx * (eta**l)
        for j in range(len(x)):
            dm2es3_dx += x[j] * m[j] * (e_ij[i, j] / T_REF) * (s_ij[i, j] ** 3)
            dm2e2s3_dx += x[j] * m[j] * ((e_ij[i, j] / T_REF) ** 2) * (s_ij[i, j] ** 3)
        dm2es3_dx *= 2.0 * m[i]
        dm2e2s3_dx *= 2.0 * m[i]
        dC1_dx = geom["c2"] * dzeta3_dx - (c1**2) * (m[i] * (8.0 * eta - 2.0 * eta * eta) / (1.0 - eta) ** 4 - m[i] * (20.0 * eta - 27.0 * eta * eta + 12.0 * eta**3 - 2.0 * eta**4) / ((1.0 - eta) * (2.0 - eta)) ** 2)
        term_dI1 = -2.0 * PI * den * dI1_dx * m2es3
        term_dm2es3 = -2.0 * PI * den * i1 * dm2es3_dx
        prefactor = -PI * den
        term_mi_c1_i2 = prefactor * (m[i] * c1 * i2) * m2e2s3
        term_mavg_dc1_i2 = prefactor * (m_avg * dC1_dx * i2) * m2e2s3
        term_mavg_c1_di2 = prefactor * (m_avg * c1 * dI2_dx) * m2e2s3
        term_mavg_c1_i2_dm2 = prefactor * (m_avg * c1 * i2) * dm2e2s3_dx
        dadx = term_dI1 + term_dm2es3 + term_mi_c1_i2 + term_mavg_dc1_i2 + term_mavg_c1_di2 + term_mavg_c1_i2_dm2
        mu = float(np.asarray(terms["mu_disp"], dtype=float)[i])
        paper_mu = _paper_value(frame, "disp", ion)
        rows.append({
            "ion": ion, "paper_mu": paper_mu, "epcsaft_mu": _kjmol(mu), "epcsaft_composition_derivative_residual_helmholtz": _kjmol(dadx),
            "resid_a_plus_z_minus_sumxj_dadx": _kjmol(float(terms["a_disp"]) + float(terms["z_disp"]) - float(terms["sum_x_dadx_disp"])),
            "term_dI1": _kjmol(term_dI1), "term_dm2es3": _kjmol(term_dm2es3), "term_mi_c1_i2": _kjmol(term_mi_c1_i2),
            "term_mavg_dc1_i2": _kjmol(term_mavg_dc1_i2), "term_mavg_c1_di2": _kjmol(term_mavg_c1_di2),
            "term_mavg_c1_i2_dm2e2s3": _kjmol(term_mavg_c1_i2_dm2), "dadx_manual": _kjmol(dadx),
            "manual_minus_api_dadx": _kjmol(dadx - float(np.asarray(terms["dadx_disp"], dtype=float)[i])),
            "paper_minus_dadx": paper_mu - _kjmol(dadx), "paper_minus_mu": paper_mu - _kjmol(mu),
        })
    return rows


def _assoc_breakdown(state: dict[str, object], geom: dict[str, object], frame: common.Table) -> list[dict[str, object]]:
    params = geom["params"]
    x = np.asarray(geom["x"], dtype=float)
    den = float(geom["den"])
    zeta = np.asarray(geom["zeta"], dtype=float)
    d = np.asarray(geom["d"], dtype=float)
    s_ij = np.asarray(geom["s_ij"], dtype=float)
    ghs = np.asarray(geom["ghs"], dtype=float)
    assoc_num = np.asarray(params.get("assoc_num", []), dtype=int)
    if assoc_num.size == 0 or not np.any(assoc_num):
        return []
    ncomp = len(x)
    e_assoc = np.asarray(params["e_assoc"], dtype=float)
    vol_a = np.asarray(params["vol_a"], dtype=float)
    k_hb = _matrix_param(params, "k_hb", ncomp, default=0.0)
    assoc_matrix = np.asarray(params["assoc_matrix"], dtype=int)
    num_sites = int(np.sum(assoc_num))
    assoc_matrix = assoc_matrix.reshape((num_sites, num_sites))
    iA: list[int] = []
    for comp_idx, site_count in enumerate(assoc_num):
        iA.extend([comp_idx] * int(site_count))
    x_assoc = np.asarray([x[idx] for idx in iA], dtype=float)
    delta_ij = np.zeros((num_sites, num_sites), dtype=float)
    ddelta_dx = np.zeros((ncomp, num_sites, num_sites), dtype=float)
    for si in range(num_sites):
        ci = iA[si]
        for sj in range(num_sites):
            if assoc_matrix[si, sj] == 0:
                continue
            cj = iA[sj]
            eABij = 0.5 * (e_assoc[ci] + e_assoc[cj])
            volABij = math.sqrt(vol_a[ci] * vol_a[cj]) * (math.sqrt(s_ij[ci, ci] * s_ij[cj, cj]) / (0.5 * (s_ij[ci, ci] + s_ij[cj, cj]))) ** 3
            if np.any(k_hb):
                volABij *= 1.0 - k_hb[ci, cj]
            delta_ij[si, sj] = ghs[ci, cj] * (math.exp(eABij / T_REF) - 1.0) * (s_ij[ci, cj] ** 3) * volABij
            for k in range(ncomp):
                dghsd_dx = PI / 6.0 * geom["m"][k] * (
                    (d[k] ** 3) / (1.0 - zeta[3]) ** 2
                    + 3.0 * d[ci] * d[cj] / (d[ci] + d[cj]) * ((d[k] ** 2) / (1.0 - zeta[3]) ** 2 + 2.0 * (d[k] ** 3) * zeta[2] / (1.0 - zeta[3]) ** 3)
                    + 2.0 * ((d[ci] * d[cj] / (d[ci] + d[cj])) ** 2) * (2.0 * (d[k] ** 2) * zeta[2] / (1.0 - zeta[3]) ** 3 + 3.0 * (d[k] ** 3) * (zeta[2] ** 2) / (1.0 - zeta[3]) ** 4)
                )
                ddelta_dx[k, si, sj] = dghsd_dx * (math.exp(eABij / T_REF) - 1.0) * (s_ij[ci, cj] ** 3) * volABij
    xa = np.zeros(num_sites, dtype=float)
    for si in range(num_sites):
        diag = delta_ij[si, si]
        if abs(diag) < 1.0e-30:
            xa[si] = 0.02
        else:
            xa[si] = (-1.0 + math.sqrt(1.0 + 8.0 * den * diag)) / (4.0 * den * diag)
            if not np.isfinite(xa[si]):
                xa[si] = 0.02
    xa_old = xa.copy()
    ctr = 0
    dif = 1000.0
    while ctr < 100 and dif > 1.0e-15:
        ctr += 1
        xa_new = xa_old.copy()
        for si in range(num_sites):
            xa_new[si] = 1.0 / (1.0 + sum(den * x_assoc[sj] * xa_old[sj] * delta_ij[si, sj] for sj in range(num_sites)))
        dif = float(np.sum(np.abs(xa_new - xa_old)))
        xa_old = 0.5 * (xa_new + xa_old)
        xa = xa_new
    dim = num_sites * ncomp
    amat = np.zeros((dim, dim), dtype=float)
    bvec = np.zeros(dim, dtype=float)
    idx1 = 0
    row = 0
    for comp in range(ncomp):
        for site_j in range(num_sites):
            sum1 = 0.0
            for site_k in range(num_sites):
                sum1 += den * x_assoc[site_k] * xa[site_k] * ddelta_dx[comp, site_j, site_k]
                amat[row, comp * num_sites + site_k] = (xa[site_j] ** 2) * den * x_assoc[site_k] * delta_ij[site_j, site_k]
            sum2 = 0.0
            for l in range(int(assoc_num[comp])):
                sum2 += xa[idx1 + l] * delta_ij[idx1 + l, site_j]
            amat[row, row] += 1.0
            bvec[row] = -(xa[site_j] ** 2) * (sum1 + sum2)
            row += 1
        idx1 += int(assoc_num[comp])
    dxa_dx = np.linalg.solve(amat, bvec).reshape((ncomp, num_sites))
    terms = state["terms"]
    rows: list[dict[str, object]] = []
    zvals = np.asarray(state["params"]["z"], dtype=float)
    for comp, ion in enumerate(state["species"]):
        if abs(zvals[comp]) <= 1.0e-12:
            continue
        chain_terms = [x[iA[site]] * den * dxa_dx[comp, site] * (1.0 / xa[site] - 0.5) for site in range(num_sites)]
        site_terms: list[float] = []
        running = 0
        for comp_idx, site_count in enumerate(assoc_num):
            if comp_idx == comp:
                for local_site in range(int(site_count)):
                    site_idx = running + local_site
                    site_terms.append(math.log(xa[site_idx]) - 0.5 * xa[site_idx] + 0.5)
            running += int(site_count)
        dadx = float(np.sum(chain_terms) + np.sum(site_terms))
        mu = float(np.asarray(terms["mu_assoc"], dtype=float)[comp])
        row_data: dict[str, object] = {
            "ion": ion, "paper_mu": _paper_value(frame, "assoc", ion), "epcsaft_mu": _kjmol(mu), "epcsaft_composition_derivative_residual_helmholtz": _kjmol(dadx),
            "resid_a_plus_z_minus_sumxj_dadx": _kjmol(float(terms["a_assoc"]) + float(terms["z_assoc"]) - float(terms["sum_x_dadx_assoc"])),
            "chain_term_sum": _kjmol(float(np.sum(chain_terms))), "site_term_sum": _kjmol(float(np.sum(site_terms))),
            "dadx_manual": _kjmol(dadx), "manual_minus_api_dadx": _kjmol(dadx - float(np.asarray(terms["dadx_assoc"], dtype=float)[comp])),
            "paper_minus_dadx": _paper_value(frame, "assoc", ion) - _kjmol(dadx), "paper_minus_mu": _paper_value(frame, "assoc", ion) - _kjmol(mu),
        }
        for site_idx in range(num_sites):
            row_data[f"XA_site_{site_idx+1}"] = float(xa[site_idx])
            row_data[f"dXA_dx_site_{site_idx+1}"] = float(dxa_dx[comp, site_idx])
            row_data[f"chain_site_{site_idx+1}"] = _kjmol(chain_terms[site_idx])
            row_data[f"self_site_{site_idx+1}"] = _kjmol(site_terms[site_idx]) if site_idx < len(site_terms) else 0.0
        rows.append(row_data)
    return rows


def main() -> None:
    frame = common.load_indexed_csv(DATA_PATH)
    bookkeeping_rows: list[dict[str, object]] = []
    hc_rows: list[dict[str, object]] = []
    disp_rows: list[dict[str, object]] = []
    assoc_rows: list[dict[str, object]] = []
    dh_rows: list[dict[str, object]] = []
    born_rows: list[dict[str, object]] = []
    for ion in frame.columns:
        state = _state_for_ion(ion)
        geom = _mixture_geometry(state)
        ion_rows = [next(row for row in _generic_contribution_rows(state, frame, contr) if row["ion"] == ion) for contr in CONTRIBUTION_MAP]
        bookkeeping_rows.extend(ion_rows)
        total_lnfug_red = float(np.asarray(state["terms"]["lnfugcoef_total"], dtype=float)[state["ion_idx"]])
        bookkeeping_rows.append(_total_row(ion, ion_rows, total_lnfug_red))
        hc_rows.extend([row for row in _hc_breakdown(state, geom, frame) if row["ion"] == ion])
        disp_rows.extend([row for row in _disp_breakdown(state, geom, frame) if row["ion"] == ion])
        assoc_rows.extend([row for row in _assoc_breakdown(state, geom, frame) if row["ion"] == ion])
        for contribution, target_rows in (("dh", dh_rows), ("born", born_rows)):
            suffix = CONTRIBUTION_MAP[contribution]["suffix"]
            idx = state["ion_idx"]
            paper_mu = 0.0 if contribution == "dh" else _paper_value(frame, contribution, ion)
            target_rows.append({
                "ion": ion, "paper_mu": paper_mu, "epcsaft_mu": _kjmol(float(np.asarray(state["terms"][f"mu_{suffix}"], dtype=float)[idx])),
                "a": _kjmol(float(state["terms"][f"a_{suffix}"])), "z": _kjmol(float(state["terms"][f"z_{suffix}"])),
                "dadx": _kjmol(float(np.asarray(state["terms"][f"dadx_{suffix}"], dtype=float)[idx])), "sum_xj_dadx": _kjmol(float(state["terms"][f"sum_x_dadx_{suffix}"])),
                "resid_a_plus_z_minus_sumxj_dadx": _kjmol(float(state["terms"][f"a_{suffix}"]) + float(state["terms"][f"z_{suffix}"]) - float(state["terms"][f"sum_x_dadx_{suffix}"])),
                "paper_minus_dadx": paper_mu - _kjmol(float(np.asarray(state["terms"][f"dadx_{suffix}"], dtype=float)[idx])), "paper_minus_mu": paper_mu - _kjmol(float(np.asarray(state["terms"][f"mu_{suffix}"], dtype=float)[idx])),
            })
    _write_csv(OUTPUT_BOOKKEEPING, ["ion", "contr", "paper_mu", "epcsaft_mu", "epcsaft_mu_manual", "a", "z", "dadx", "sum_xj_dadx", "z_quot", "epcsaft_lnfug", "epcsaft_lnfug_manual", "lnfug_red", "lnact_star_red", "gamma_star", "Ghyd", "manual_minus_epcsaft_mu", "manual_minus_epcsaft_lnfug", "paper_minus_epcsaft_mu", "epcsaft_total_lnfug_from_api"], bookkeeping_rows)
    _write_csv(OUTPUT_HC, ["ion", "paper_mu", "epcsaft_mu", "epcsaft_composition_derivative_residual_helmholtz", "resid_a_plus_z_minus_sumxj_dadx", "m_i_ares_hs", "m_avg_dahs_dx", "minus_chain_sum", "minus_log_term", "dadx_manual", "manual_minus_api_dadx", "paper_minus_dadx", "paper_minus_mu"], hc_rows)
    _write_csv(OUTPUT_DISP, ["ion", "paper_mu", "epcsaft_mu", "epcsaft_composition_derivative_residual_helmholtz", "resid_a_plus_z_minus_sumxj_dadx", "term_dI1", "term_dm2es3", "term_mi_c1_i2", "term_mavg_dc1_i2", "term_mavg_c1_di2", "term_mavg_c1_i2_dm2e2s3", "dadx_manual", "manual_minus_api_dadx", "paper_minus_dadx", "paper_minus_mu"], disp_rows)
    assoc_fields = ["ion", "paper_mu", "epcsaft_mu", "epcsaft_composition_derivative_residual_helmholtz", "resid_a_plus_z_minus_sumxj_dadx", "chain_term_sum", "site_term_sum", "dadx_manual", "manual_minus_api_dadx", "paper_minus_dadx", "paper_minus_mu"]
    max_assoc_sites = max((max(int(key.split("_")[-1]) for key in row if key.startswith("XA_site_")) for row in assoc_rows), default=0)
    for idx in range(1, max_assoc_sites + 1):
        assoc_fields.extend([f"XA_site_{idx}", f"dXA_dx_site_{idx}", f"chain_site_{idx}", f"self_site_{idx}"])
    _write_csv(OUTPUT_ASSOC, assoc_fields, assoc_rows)
    _write_csv(OUTPUT_DH, ["ion", "paper_mu", "epcsaft_mu", "a", "z", "dadx", "sum_xj_dadx", "resid_a_plus_z_minus_sumxj_dadx", "paper_minus_dadx", "paper_minus_mu"], dh_rows)
    _write_csv(OUTPUT_BORN, ["ion", "paper_mu", "epcsaft_mu", "a", "z", "dadx", "sum_xj_dadx", "resid_a_plus_z_minus_sumxj_dadx", "paper_minus_dadx", "paper_minus_mu"], born_rows)
    print(f"Wrote {OUTPUT_BOOKKEEPING}", flush=True)
    print(f"Wrote {OUTPUT_HC}", flush=True)
    print(f"Wrote {OUTPUT_DISP}", flush=True)
    print(f"Wrote {OUTPUT_ASSOC}", flush=True)
    print(f"Wrote {OUTPUT_DH}", flush=True)
    print(f"Wrote {OUTPUT_BORN}", flush=True)


if __name__ == "__main__":
    main()

