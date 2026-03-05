"""Reproduce Held 2014-style NaCl/KBr osmotic plot with 2008 and 2014 versions."""

import csv
import sys
from pathlib import Path

import matplotlib
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from pcsaft import pcsaft_den, pcsaft_fugcoef, pcsaft_p

matplotlib.use("Agg")
import matplotlib.pyplot as plt


T_REF = 298.15
P_REF = 1.0e5
DIELC_CONST = 78.09
MW_WATER = 18.0153e-3  # kg/mol

def _elec_model_no_born_constant():
    return {
        "rel_perm": {"rule": 0, "differential_mode": "analytical"},
        "DH_model": {"d_ion_mode": 1, "bjeruum_treatment": False},
        "include_born_model": False,
        "born_model": {
            "d_Born_mode": 0,
            "solvation_shell_model": False,
            "dielectric_saturation": False,
            "bulk_mode": "mix",
            "mu_born_model": {
                "differential_mode": "analytical",
                "comp_dep_rel_perm": True,
                "include_sum_term": True,
                "comp_dep_delta_d": False,
            },
        },
    }


def _water_sigma(t):
    return 2.7927 + 10.11 * np.exp(-0.01775 * t) - 1.417 * np.exp(-0.01146 * t)


def _mole_fraction_from_molality_11(molality):
    """1:1 electrolyte in water: [cation, anion, water] mole fractions."""
    n_water = 1.0 / MW_WATER
    n_cation = molality
    n_anion = molality
    n_total = n_water + n_cation + n_anion
    return np.asarray([n_cation, n_anion, n_water], dtype=float) / n_total


def _species_for_salt(salt):
    mapping = {
        "NaCl": ["Na+", "Cl-", "H2O"],
        "KBr": ["K+", "Br-", "H2O"],
    }
    if salt not in mapping:
        raise ValueError(f"Unsupported salt: {salt}")
    return mapping[salt]


def _load_osmotic_data(salt, m_min=0.0, m_max=4.0):
    root = REPO_ROOT
    candidates = [
        root / "data" / "osmotic" / "water" / f"{salt}.csv",
        root / "data" / "osmotic" / f"{salt}.csv",
    ]
    data_path = next((path for path in candidates if path.exists()), None)
    if data_path is None:
        raise FileNotFoundError(f"Data file not found in candidates: {candidates}")

    m_vals, phi_vals = [], []
    with data_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"No header found in {data_path}")
        m_key = next((k for k in reader.fieldnames if k and "m" in k.lower()), None)
        phi_key = next((k for k in reader.fieldnames if k and "osmotic" in k.lower()), None)
        if m_key is None or phi_key is None:
            raise KeyError(f"Expected molality/osmotic columns in {data_path}, got {reader.fieldnames}")
        for row in reader:
            m = float(row[m_key])
            phi = float(row[phi_key])
            if m_min <= m <= m_max:
                m_vals.append(m)
                phi_vals.append(phi)

    if not m_vals:
        raise ValueError(f"No osmotic data found for {salt} in {m_min}<=m<={m_max}.")

    m_arr = np.asarray(m_vals, dtype=float)
    phi_arr = np.asarray(phi_vals, dtype=float)
    order = np.argsort(m_arr)
    return m_arr[order], phi_arr[order]


def _build_params_2008(salt, t):
    """Strategy-1 style parameters from Held et al. (2008) values in Held 2014 Table 2 (left)."""
    sigma_w = _water_sigma(t)

    if salt == "NaCl":
        m = np.asarray([1.0, 1.0, 1.2047], dtype=float)
        s = np.asarray([2.4122, 3.0575, sigma_w], dtype=float)
        e = np.asarray([646.05, 47.29, 353.9449], dtype=float)
    elif salt == "KBr":
        m = np.asarray([1.0, 1.0, 1.2047], dtype=float)
        s = np.asarray([2.9698, 3.4573, sigma_w], dtype=float)
        e = np.asarray([271.05, 60.22, 353.9449], dtype=float)
    else:
        raise ValueError(f"Unsupported salt: {salt}")

    k_ij = np.zeros((3, 3), dtype=float)
    # Strategy-1 should not include cation-anion dispersion; kij=1 removes cross-dispersion.
    k_ij[0, 1] = 1.0
    k_ij[1, 0] = 1.0

    params = {
        "m": m,
        "s": s,
        "e": e,
        "e_assoc": np.asarray([0.0, 0.0, 2425.7], dtype=float),
        "vol_a": np.asarray([0.0, 0.0, 0.0451], dtype=float),
        "z": np.asarray([1.0, -1.0, 0.0], dtype=float),
        "k_ij": k_ij,
        "dielc": np.full(3, DIELC_CONST, dtype=float),
        "elec_model": _elec_model_no_born_constant(),
        "debug": False,
    }
    return params


def _build_params_2014(salt, t):
    """Strategy-2 style parameters from Held 2014 Table 2 (right) + Table 3 ion-ion kij."""
    sigma_w = _water_sigma(t)

    if salt == "NaCl":
        m = np.asarray([1.0, 1.0, 1.2047], dtype=float)
        s = np.asarray([2.8232, 2.7560, sigma_w], dtype=float)
        e = np.asarray([230.00, 170.00, 353.9449], dtype=float)
        k_ij = np.zeros((3, 3), dtype=float)
        k_na_w = -0.007981 * t + 2.37999
        k_ij[0, 2] = k_na_w
        k_ij[2, 0] = k_na_w
        k_ij[1, 2] = -0.25
        k_ij[2, 1] = -0.25
        k_ij[0, 1] = 0.317
        k_ij[1, 0] = 0.317
    elif salt == "KBr":
        m = np.asarray([1.0, 1.0, 1.2047], dtype=float)
        s = np.asarray([3.3417, 3.0707, sigma_w], dtype=float)
        e = np.asarray([200.00, 190.00, 353.9449], dtype=float)
        k_ij = np.zeros((3, 3), dtype=float)
        k_k_w = -0.004012 * t + 1.3959
        k_ij[0, 2] = k_k_w
        k_ij[2, 0] = k_k_w
        k_ij[1, 2] = -0.25
        k_ij[2, 1] = -0.25
        k_ij[0, 1] = -0.102
        k_ij[1, 0] = -0.102
    else:
        raise ValueError(f"Unsupported salt: {salt}")

    params = {
        "m": m,
        "s": s,
        "e": e,
        "e_assoc": np.asarray([0.0, 0.0, 2425.7], dtype=float),
        "vol_a": np.asarray([0.0, 0.0, 0.0451], dtype=float),
        "z": np.asarray([1.0, -1.0, 0.0], dtype=float),
        "k_ij": k_ij,
        "dielc": np.full(3, DIELC_CONST, dtype=float),
        "elec_model": _elec_model_no_born_constant(),
        "debug": False,
    }
    return params


def _calc_osmotic_curve(salt, m_values, strategy):
    builders = {
        "2008": _build_params_2008,
        "2014": _build_params_2014,
    }
    if strategy not in builders:
        raise ValueError(f"Unsupported strategy: {strategy}")
    params = builders[strategy](salt, T_REF)
    phi_calc = np.empty_like(m_values, dtype=float)

    for i, m_salt in enumerate(m_values):
        # Avoid singular m=0 in osmotic coefficient expression.
        m_eval = max(float(m_salt), 1e-12)
        x = _mole_fraction_from_molality_11(m_eval)
        rho = pcsaft_den(T_REF, P_REF, x, params, phase="liq")
        phi_calc[i] = _osmotic_molality_from_fugacity(T_REF, rho, x, params)

    if not np.all(np.isfinite(phi_calc)):
        raise ValueError(f"Non-finite osmotic results for {salt} ({strategy}).")
    return phi_calc


def _osmotic_molality_from_fugacity(t, rho, x, params):
    z = np.asarray(params["z"], dtype=float)
    idx_water = np.where(np.abs(z) <= 1e-12)[0]
    if idx_water.size != 1:
        raise ValueError(f"Expected exactly one neutral solvent for osmotic conversion, got {idx_water.size}.")
    iw = int(idx_water[0])

    molality = x / (x[iw] * 18.0153 / 1000.0)
    molality[iw] = 0.0

    x0 = np.zeros_like(x)
    x0[iw] = 1.0

    fugcoef = np.asarray(pcsaft_fugcoef(t, rho, x, params), dtype=float).reshape(-1)
    p_mix = pcsaft_p(t, rho, x, params)
    phase0 = "vap" if rho < 900.0 else "liq"
    rho0 = pcsaft_den(t, p_mix, x0, params, phase=phase0)
    fugcoef0 = np.asarray(pcsaft_fugcoef(t, rho0, x0, params), dtype=float).reshape(-1)
    gamma_w = fugcoef[iw] / fugcoef0[iw]

    phi_m = -1000.0 * np.log(x[iw] * gamma_w) / 18.0153 / np.sum(molality)
    return float(phi_m)


def run_validation_2014_repro():
    out_dir = REPO_ROOT / "data" / "osmotic" / "water" / "plot_fits"
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_path = out_dir / "validation_2014_repro_NaCl_KBr_fit.png"
    # Start above zero to avoid the m -> 0 osmotic singularity in model curves.
    m_plot = np.linspace(0.01, 4.0, 101)

    data = {}
    marker_map = {"NaCl": "o", "KBr": "s"}
    line_styles = {"NaCl": "--", "KBr": "--"}
    for salt in ("NaCl", "KBr"):
        m_exp, phi_exp = _load_osmotic_data(salt)
        data[salt] = {
            "m_exp": m_exp,
            "phi_exp": phi_exp,
            "phi_2008": _calc_osmotic_curve(salt, m_plot, "2008"),
            "phi_2014": _calc_osmotic_curve(salt, m_plot, "2014"),
        }

    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    for salt in ("NaCl", "KBr"):
        ax.scatter(
            data[salt]["m_exp"],
            data[salt]["phi_exp"],
            marker=marker_map[salt],
            s=36,
            color="black",
            facecolors="none",
            label=f"{salt} data",
        )
        ax.plot(
            m_plot,
            data[salt]["phi_2008"],
            color="gray",
            linewidth=1.6,
            linestyle=line_styles[salt],
            label=f"{salt} 2008",
        )
        ax.plot(
            m_plot,
            data[salt]["phi_2014"],
            color="black",
            linewidth=2.0,
            linestyle=line_styles[salt],
            label=f"{salt} 2014",
        )

    ax.set_xlim(0.0, 4.0)
    ax.set_ylim(0.8, 1.2)
    ax.set_xlabel("molality, m / mol kg$^{-1}$")
    ax.set_ylabel("molal osmotic coefficient, $\\phi_m$")
    ax.set_title("Aqueous NaCl and KBr molal osmotic coefficients at 298.15 K (2008/2014)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(plot_path, dpi=220)
    plt.close(fig)

    if not plot_path.exists():
        raise FileNotFoundError(f"Expected plot was not written: {plot_path}")


def test_validation_2014_repro():
    run_validation_2014_repro()

if __name__ == "__main__":
    run_validation_2014_repro()
