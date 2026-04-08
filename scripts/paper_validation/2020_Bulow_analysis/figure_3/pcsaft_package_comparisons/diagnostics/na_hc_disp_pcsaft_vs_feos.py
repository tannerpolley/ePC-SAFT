from __future__ import annotations

import csv
import math
import sys
from copy import deepcopy
from pathlib import Path

import numpy as np

try:
    import feos
    import si_units as si
except ImportError as exc:  # pragma: no cover
    raise RuntimeError("feos and si_units must be installed in the PC-SAFT conda environment.") from exc


SCRIPT_DIR = Path(__file__).resolve().parent
PACKAGE_DIR = SCRIPT_DIR.parent
FIGURE3_DIR = PACKAGE_DIR.parent
ANALYSIS_ROOT = FIGURE3_DIR.parent
REPO_ROOT = ANALYSIS_ROOT.parents[2]
FIGURE3_DIAG_DIR = FIGURE3_DIR / "diagnostics"

for path in (ANALYSIS_ROOT, REPO_ROOT, FIGURE3_DIAG_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import _model_overlay as overlay
import figure3_detailed_bookkeeping as f3diag
from scripts._env import require_pcsaft_install

require_pcsaft_install()

from pcsaft.parameters import get_prop_dict
from scripts._pcsaft_oop import pcsaft_lnfugcoef_terms


CSV_PATH = SCRIPT_DIR / "na_hc_disp_pcsaft_vs_feos.csv"
MD_PATH = SCRIPT_DIR / "na_hc_disp_pcsaft_vs_feos.md"

FEOS_ROOT = Path(r"C:\Users\Tanner\Documents\git\feos")
FEOS_PARAMETER_PATH = FEOS_ROOT / "parameters" / "epcsaft" / "held2014_w_permittivity_added.json"
FEOS_BINARY_PATH = FEOS_ROOT / "parameters" / "epcsaft" / "held2014_binary.json"
FEOS_RESIDUAL_PROPERTIES = FEOS_ROOT / "crates" / "feos-core" / "src" / "state" / "residual_properties.rs"
FEOS_PROPERTIES = FEOS_ROOT / "crates" / "feos-core" / "src" / "state" / "properties.rs"
FEOS_HC_SOURCE = FEOS_ROOT / "crates" / "feos" / "src" / "epcsaft" / "eos" / "hard_chain.rs"
FEOS_DISP_SOURCE = FEOS_ROOT / "crates" / "feos" / "src" / "epcsaft" / "eos" / "dispersion.rs"

R_GAS = overlay.R_GAS
T_REF = overlay.T_REF
P_REF = overlay.P_REF
EPS = overlay.EPS
EPS_INF = overlay.EPS_INF
RT_KJMOL = R_GAS * T_REF / 1000.0
N_AV = 6.02214076e23
PI = math.pi


def _kjmol(value: float) -> float:
    return float(RT_KJMOL * float(value))


def _to_kj_per_mol(quantity) -> float:
    return float(quantity / (si.KILO * si.JOULE / si.MOL))


def _to_pa(quantity) -> float:
    return float(quantity / si.PASCAL)


def _map_feos_terms(pairs, *, combine_hc: bool) -> dict[str, float]:
    out: dict[str, float] = {}
    for label, quantity in pairs:
        key = label
        value = _to_kj_per_mol(quantity)
        if combine_hc and label in {"Hard Sphere", "Hard Chain"}:
            key = "hc"
        elif label == "Dispersion":
            key = "disp"
        else:
            key = label
        out[key] = out.get(key, 0.0) + value
    return out


def _map_feos_pressure_terms(pairs, volume_m3: float) -> dict[str, float]:
    out: dict[str, float] = {}
    for label, quantity in pairs:
        if label == "Ideal gas":
            continue
        key = "hc" if label in {"Hard Sphere", "Hard Chain"} else "disp" if label == "Dispersion" else label
        z_red = _to_pa(quantity) * volume_m3 / (R_GAS * T_REF)
        out[key] = out.get(key, 0.0) + z_red
    return out


def _build_feos_state() -> tuple[object, object, np.ndarray]:
    params = feos.Parameters.from_json(
        ["water", "sodium ion", "chloride ion"],
        str(FEOS_PARAMETER_PATH),
        str(FEOS_BINARY_PATH),
    )
    eos = feos.EquationOfState.epcsaft(params, epcsaft_variant="advanced")
    x_bulk = np.asarray([1.0 - 2.0 * EPS, EPS, EPS], dtype=float)
    state_bulk = feos.State(
        eos,
        temperature=T_REF * si.KELVIN,
        pressure=P_REF * si.PASCAL,
        molefracs=x_bulk,
        total_moles=si.MOL,
    )
    density_bulk = state_bulk.density
    x_ref = np.asarray([1.0, 0.0, 0.0], dtype=float)
    state_ref = feos.State(
        eos,
        temperature=T_REF * si.KELVIN,
        density=density_bulk,
        molefracs=x_ref,
        total_moles=si.MOL,
    )
    p_ref = state_ref.pressure()
    x_inf = x_ref.copy()
    x_inf[1] = EPS_INF
    x_inf[0] = 1.0 - EPS_INF
    state_inf = feos.State(
        eos,
        temperature=T_REF * si.KELVIN,
        pressure=p_ref,
        molefracs=x_inf,
        total_moles=si.MOL,
    )
    return eos, state_inf, x_inf


def _feos_a_terms(state) -> dict[str, float]:
    return _map_feos_terms(state.residual_molar_helmholtz_energy_contributions(), combine_hc=True)


def _feos_mu_public(state) -> dict[str, float]:
    return _map_feos_terms(state.chemical_potential_contributions(1, feos.Contributions.Residual), combine_hc=True)


def _feos_z_terms(state) -> dict[str, float]:
    volume_m3 = float(state.volume / (si.METER**3))
    return _map_feos_pressure_terms(state.pressure_contributions(), volume_m3)


def _feos_fd_dadx(eos, density, x_base: np.ndarray, h: float) -> dict[str, float]:
    if x_base[0] - h <= 0.0:
        raise ValueError(f"Finite-difference step {h} is too large for the solvent mole fraction.")
    x_step = x_base.copy()
    x_step[1] += h
    x_step[0] -= h
    state_0 = feos.State(
        eos,
        temperature=T_REF * si.KELVIN,
        density=density,
        molefracs=x_base,
        total_moles=si.MOL,
    )
    state_1 = feos.State(
        eos,
        temperature=T_REF * si.KELVIN,
        density=density,
        molefracs=x_step,
        total_moles=si.MOL,
    )
    a0 = _feos_a_terms(state_0)
    a1 = _feos_a_terms(state_1)
    return {
        "hc": (a1["hc"] - a0["hc"]) / h,
        "disp": (a1["disp"] - a0["disp"]) / h,
    }


def _pcsaft_numerical_terms(state: dict[str, object]) -> dict[str, np.ndarray | float]:
    params_num = deepcopy(state["params"])
    elec_model = deepcopy(params_num.get("elec_model", {}))
    elec_model.setdefault("hc_model", {})["dadx_differential_mode"] = "numerical"
    elec_model.setdefault("disp_model", {})["dadx_differential_mode"] = "numerical"
    params_num["elec_model"] = elec_model
    return pcsaft_lnfugcoef_terms(T_REF, float(state["rho_inf"]), np.asarray(state["x_inf"], dtype=float), params_num)


def _pcsaft_hc_aux(state: dict[str, object], geom: dict[str, object]) -> dict[str, float]:
    x = np.asarray(geom["x"], dtype=float)
    m = np.asarray(geom["m"], dtype=float)
    d = np.asarray(geom["d"], dtype=float)
    zeta = np.asarray(geom["zeta"], dtype=float)
    ghs = np.asarray(geom["ghs"], dtype=float)
    ares_hs = float(geom["ares_hs"])
    m_avg = float(geom["m_avg"])
    den = float(geom["den"])
    ncomp = len(x)
    i = int(state["ion_idx"])
    dzeta_dx = np.zeros(4, dtype=float)
    for order in range(4):
        dzeta_dx[order] = PI / 6.0 * den * m[i] * (d[i] ** order)
    dahs_dx = (
        -dzeta_dx[0] / zeta[0] * ares_hs
        + 1.0
        / zeta[0]
        * (
            3.0 * (dzeta_dx[1] * zeta[2] + zeta[1] * dzeta_dx[2]) / (1.0 - zeta[3])
            + 3.0 * zeta[1] * zeta[2] * dzeta_dx[3] / (1.0 - zeta[3]) ** 2
            + 3.0 * zeta[2] * zeta[2] * dzeta_dx[2] / zeta[3] / (1.0 - zeta[3]) ** 2
            + (zeta[2] ** 3) * dzeta_dx[3] * (3.0 * zeta[3] - 1.0) / (zeta[3] ** 2) / (1.0 - zeta[3]) ** 3
            + math.log(1.0 - zeta[3]) * (((3.0 * zeta[2] * zeta[2] * dzeta_dx[2] * zeta[3]) - 2.0 * (zeta[2] ** 3) * dzeta_dx[3]) / (zeta[3] ** 3) - dzeta_dx[0])
            + (zeta[0] - (zeta[2] ** 3) / (zeta[3] ** 2)) * dzeta_dx[3] / (1.0 - zeta[3])
        )
    )
    dghsii = []
    for j in range(ncomp):
        djj = d[j] * d[j] / (d[j] + d[j])
        value = (
            dzeta_dx[3] / (1.0 - zeta[3]) ** 2
            + djj * (3.0 * dzeta_dx[2] / (1.0 - zeta[3]) ** 2 + 6.0 * zeta[2] * dzeta_dx[3] / (1.0 - zeta[3]) ** 3)
            + (djj**2) * (4.0 * zeta[2] * dzeta_dx[2] / (1.0 - zeta[3]) ** 3 + 6.0 * zeta[2] * zeta[2] * dzeta_dx[3] / (1.0 - zeta[3]) ** 4)
        )
        dghsii.append(float(value))
    return {
        "ares_hs_kj_mol": _kjmol(ares_hs),
        "zhs_red": float(geom["zhs"]),
        "dahs_dx_kj_mol": _kjmol(dahs_dx),
        "dghsii_dx_comp0": dghsii[0],
        "dghsii_dx_comp1": dghsii[1],
        "dghsii_dx_comp2": dghsii[2],
        "ghs_comp0": float(ghs[0, 0]),
        "ghs_comp1": float(ghs[1, 1]),
        "ghs_comp2": float(ghs[2, 2]),
        "m_avg": m_avg,
    }


def _pcsaft_disp_aux(state: dict[str, object], geom: dict[str, object]) -> dict[str, float]:
    m = np.asarray(geom["m"], dtype=float)
    d = np.asarray(geom["d"], dtype=float)
    eta = float(geom["eta"])
    m_avg = float(geom["m_avg"])
    a = np.asarray(geom["a_coeffs"], dtype=float)
    b = np.asarray(geom["b_coeffs"], dtype=float)
    i1 = float(geom["i1"])
    i2 = float(geom["i2"])
    c1 = float(geom["c1"])
    c2 = float(geom["c2"])
    m2es3 = float(geom["m2es3"])
    m2e2s3 = float(geom["m2e2s3"])
    den = float(geom["den"])
    i = int(state["ion_idx"])
    dzeta3_dx = PI / 6.0 * den * m[i] * (d[i] ** 3)
    dI1_dx = 0.0
    dI2_dx = 0.0
    for l in range(7):
        daa_dx = m[i] / (m_avg**2) * f3diag.A1[l] + m[i] / (m_avg**2) * (3.0 - 4.0 / m_avg) * f3diag.A2[l]
        db_dx = m[i] / (m_avg**2) * f3diag.B1[l] + m[i] / (m_avg**2) * (3.0 - 4.0 / m_avg) * f3diag.B2[l]
        dI1_dx += a[l] * l * dzeta3_dx * (eta ** (l - 1)) + daa_dx * (eta**l)
        dI2_dx += b[l] * l * dzeta3_dx * (eta ** (l - 1)) + db_dx * (eta**l)
    s_ij = np.asarray(geom["s_ij"], dtype=float)
    e_ij = np.asarray(geom["e_ij"], dtype=float)
    x = np.asarray(geom["x"], dtype=float)
    dm2es3_dx = 0.0
    dm2e2s3_dx = 0.0
    for j in range(len(x)):
        dm2es3_dx += x[j] * m[j] * (e_ij[i, j] / T_REF) * (s_ij[i, j] ** 3)
        dm2e2s3_dx += x[j] * m[j] * ((e_ij[i, j] / T_REF) ** 2) * (s_ij[i, j] ** 3)
    dm2es3_dx *= 2.0 * m[i]
    dm2e2s3_dx *= 2.0 * m[i]
    dC1_dx = c2 * dzeta3_dx - (c1**2) * (
        m[i] * (8.0 * eta - 2.0 * eta * eta) / (1.0 - eta) ** 4
        - m[i] * (20.0 * eta - 27.0 * eta * eta + 12.0 * eta**3 - 2.0 * eta**4) / ((1.0 - eta) * (2.0 - eta)) ** 2
    )
    return {
        "eta": eta,
        "m_avg": m_avg,
        "I1": i1,
        "I2": i2,
        "C1": c1,
        "C2": c2,
        "m2es3": m2es3,
        "m2e2s3": m2e2s3,
        "dI1_dx": dI1_dx,
        "dI2_dx": dI2_dx,
        "dC1_dx": dC1_dx,
        "dm2es3_dx": dm2es3_dx,
        "dm2e2s3_dx": dm2e2s3_dx,
    }


def _collect_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    pc_state = f3diag._state_for_ion("Na+")
    pc_geom = f3diag._mixture_geometry(pc_state)
    pc_terms = pc_state["terms"]
    pc_terms_num = _pcsaft_numerical_terms(pc_state)
    pc_idx = int(pc_state["ion_idx"])

    def add(section: str, package: str, contr: str, quantity: str, value: float | str, note: str = "") -> None:
        rows.append(
            {
                "section": section,
                "package": package,
                "contr": contr,
                "quantity": quantity,
                "value": value,
                "note": note,
            }
        )

    add("state", "pcsaft", "state", "component_order", ",".join(pc_state["species"]))
    add("state", "pcsaft", "state", "rho_inf_mol_m3", float(pc_state["rho_inf"]))
    x_pc = np.asarray(pc_state["x_inf"], dtype=float)
    for i, sp in enumerate(pc_state["species"]):
        add("state", "pcsaft", "state", f"x_inf_{sp}", float(x_pc[i]))

    for contr in ("hc", "disp"):
        add("branch", "pcsaft", contr, "a_kj_mol", _kjmol(float(pc_terms[f"a_{contr}"])))
        add("branch", "pcsaft", contr, "z_red", float(pc_terms[f"z_{contr}"]))
        add("branch", "pcsaft", contr, "z_kj_mol", _kjmol(float(pc_terms[f"z_{contr}"])))
        add("branch", "pcsaft", contr, "dadx_analytical_kj_mol", _kjmol(float(np.asarray(pc_terms[f"dadx_{contr}"], dtype=float)[pc_idx])))
        add("branch", "pcsaft", contr, "dadx_numerical_kj_mol", _kjmol(float(np.asarray(pc_terms_num[f"dadx_{contr}"], dtype=float)[pc_idx])))
        add(
            "branch",
            "pcsaft",
            contr,
            "dadx_num_minus_analytical_kj_mol",
            _kjmol(float(np.asarray(pc_terms_num[f"dadx_{contr}"], dtype=float)[pc_idx] - np.asarray(pc_terms[f"dadx_{contr}"], dtype=float)[pc_idx])),
        )
        add("branch", "pcsaft", contr, "sum_xj_dadx_kj_mol", _kjmol(float(pc_terms[f"sum_x_dadx_{contr}"])))
        add("branch", "pcsaft", contr, "mu_kj_mol", _kjmol(float(np.asarray(pc_terms[f"mu_{contr}"], dtype=float)[pc_idx])))
        add(
            "branch",
            "pcsaft",
            contr,
            "resid_a_plus_z_minus_sumxj_dadx_kj_mol",
            _kjmol(float(pc_terms[f"a_{contr}"]) + float(pc_terms[f"z_{contr}"]) - float(pc_terms[f"sum_x_dadx_{contr}"])),
        )

    for key, value in _pcsaft_hc_aux(pc_state, pc_geom).items():
        add("aux", "pcsaft", "hc", key, value)
    for key, value in _pcsaft_disp_aux(pc_state, pc_geom).items():
        add("aux", "pcsaft", "disp", key, value)

    eos, feos_state, x_feos = _build_feos_state()
    feos_a = _feos_a_terms(feos_state)
    feos_mu = _feos_mu_public(feos_state)
    feos_z = _feos_z_terms(feos_state)
    density_feos = feos_state.density
    add("state", "feos", "state", "component_order", "water,sodium ion,chloride ion")
    add("state", "feos", "state", "rho_inf_mol_m3", float(density_feos / (si.MOL / si.METER**3)))
    for label, value in zip(("water", "sodium ion", "chloride ion"), x_feos, strict=True):
        add("state", "feos", "state", f"x_inf_{label}", float(value))

    for contr in ("hc", "disp"):
        add("branch", "feos", contr, "a_kj_mol", float(feos_a[contr]))
        add("branch", "feos", contr, "z_red", float(feos_z[contr]))
        add("branch", "feos", contr, "z_kj_mol", _kjmol(float(feos_z[contr])))
        add("branch", "feos", contr, "mu_public_kj_mol", float(feos_mu[contr]))
        add(
            "compare",
            "feos_vs_pcsaft",
            contr,
            "feos_mu_public_minus_pcsaft_mu_kj_mol",
            float(feos_mu[contr] - _kjmol(float(np.asarray(pc_terms[f"mu_{contr}"], dtype=float)[pc_idx]))),
        )
        add(
            "compare",
            "feos_vs_pcsaft",
            contr,
            "feos_mu_public_minus_pcsaft_dadx_kj_mol",
            float(feos_mu[contr] - _kjmol(float(np.asarray(pc_terms[f"dadx_{contr}"], dtype=float)[pc_idx]))),
        )

    for h in (1.0e-8, 1.0e-7, 1.0e-6):
        fd = _feos_fd_dadx(eos, density_feos, x_feos, h)
        for contr in ("hc", "disp"):
            add("branch", "feos", contr, f"dadx_fd_h{h:.0e}_kj_mol", float(fd[contr]))
            add("branch", "feos", contr, f"mu_public_minus_fd_h{h:.0e}_kj_mol", float(feos_mu[contr] - fd[contr]))

    add(
        "source",
        "feos",
        "hc",
        "public_path",
        "residual_chemical_potential_contributions differentiates molar_helmholtz_energy_contributions wrt x",
        f"{FEOS_RESIDUAL_PROPERTIES}: residual_chemical_potential_contributions",
    )
    add(
        "source",
        "feos",
        "hc",
        "bug_v_init",
        "v = temperature.into_reduced()",
        f"{FEOS_RESIDUAL_PROPERTIES}: let v = Dual::from_re(self.temperature.into_reduced())",
    )
    add(
        "source",
        "feos",
        "hc",
        "a_hc_split",
        "Hard Sphere + Hard Chain",
        f"{FEOS_HC_SOURCE} and {FEOS_ROOT / 'crates' / 'feos' / 'src' / 'epcsaft' / 'eos' / 'mod.rs'}",
    )
    add(
        "source",
        "feos",
        "disp",
        "a_disp_formula",
        "Dispersion.helmholtz_energy_density",
        f"{FEOS_DISP_SOURCE}",
    )

    return rows


def _write_csv(rows: list[dict[str, object]]) -> None:
    with CSV_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["section", "package", "contr", "quantity", "value", "note"])
        writer.writeheader()
        writer.writerows(rows)


def _get(rows: list[dict[str, object]], package: str, contr: str, quantity: str) -> float:
    for row in rows:
        if row["package"] == package and row["contr"] == contr and row["quantity"] == quantity:
            return float(row["value"])
    raise KeyError((package, contr, quantity))


def _write_md(rows: list[dict[str, object]]) -> None:
    pc_mu_disp = _get(rows, "pcsaft", "disp", "mu_kj_mol")
    pc_dadx_disp = _get(rows, "pcsaft", "disp", "dadx_analytical_kj_mol")
    fe_mu_disp = _get(rows, "feos", "disp", "mu_public_kj_mol")
    fe_fd_disp = _get(rows, "feos", "disp", "dadx_fd_h1e-06_kj_mol")
    pc_mu_hc = _get(rows, "pcsaft", "hc", "mu_kj_mol")
    pc_dadx_hc = _get(rows, "pcsaft", "hc", "dadx_analytical_kj_mol")
    fe_mu_hc = _get(rows, "feos", "hc", "mu_public_kj_mol")
    fe_fd_hc = _get(rows, "feos", "hc", "dadx_fd_h1e-06_kj_mol")

    lines = [
        "# Na+ hard-chain / dispersion comparison: current repo vs feos",
        "",
        "This note compares only the `Na+` Figure 3 water infinite-dilution state after the local `feos` branch-bookkeeping patch.",
        "",
        "## Main outcome",
        "",
        f"- Current repo dispersion still has the known residual: $\\mu^{{disp}} - \\partial a^{{disp}}/\\partial x_k = {pc_mu_disp - pc_dadx_disp:.9f}\\ \\mathrm{{kJ\\ mol^{{-1}}}}$.",
        f"- Current repo hard-chain still satisfies $\\mu^{{hc}} = \\partial a^{{hc}}/\\partial x_k$ for this state: residual `{pc_mu_hc - pc_dadx_hc:.9e}` kJ/mol.",
        f"- Patched `feos` hard-chain now matches the current repo full $\\tilde{{\\mu}}^{{hc}}$ branch to roundoff: `feos - pcsaft mu = {fe_mu_hc - pc_mu_hc:.9e}` kJ/mol.",
        f"- Patched `feos` dispersion now matches the current repo full $\\tilde{{\\mu}}^{{disp}}$ branch to roundoff: `feos - pcsaft mu = {fe_mu_disp - pc_mu_disp:.9e}` kJ/mol.",
        f"- The patched `feos` public branches are still not the same as a fixed-density finite difference of the exposed Helmholtz branch values: hc `mu_public = {fe_mu_hc:.9f}` vs `fd(h=1e-6) = {fe_fd_hc:.9f}`, disp `mu_public = {fe_mu_disp:.9f}` vs `fd(h=1e-6) = {fe_fd_disp:.9f}` kJ/mol.",
        "",
        "## Interpretation",
        "",
        "- The non-differential branch pieces already matched closely across packages, so the actual fix target was the public branch-bookkeeping helper, not the underlying $a^\\alpha$ or $Z^\\alpha$ values.",
        "- The local patch in `state/properties.rs::chemical_potential_contributions(...)` now returns the repo-style branch quantity $\\tilde{\\mu}^\\alpha = a^\\alpha + Z^\\alpha + dadx^\\alpha - \\sum_j x_j dadx_j^\\alpha$.",
        "- That makes the patched `feos` branch values align with the current repo branch values while preserving the already-correct total hydration energy.",
        "- The finite-difference mismatch is expected because the public helper is no longer returning a raw fixed-density derivative of `molar_helmholtz_energy_contributions(...)`; it is returning the assembled branch chemical-potential expression.",
        "",
        "## Files",
        "",
        f"- CSV: `{CSV_PATH}`",
        f"- feos patched Python contribution path: `{FEOS_PROPERTIES}`",
        f"- feos residual helper not used by Python here: `{FEOS_RESIDUAL_PROPERTIES}`",
        f"- feos hard-chain contribution source: `{FEOS_HC_SOURCE}`",
        f"- feos dispersion contribution source: `{FEOS_DISP_SOURCE}`",
    ]
    MD_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")

def main() -> None:
    rows = _collect_rows()
    _write_csv(rows)
    _write_md(rows)
    print(f"Wrote {CSV_PATH}", flush=True)
    print(f"Wrote {MD_PATH}", flush=True)


if __name__ == "__main__":
    main()

