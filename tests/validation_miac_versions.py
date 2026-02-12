"""NaBr alcohol MIAC comparison across 2008/2014/2020 model versions."""

import csv
import sys
from pathlib import Path

import matplotlib
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pcsaft import pcsaft_den, pcsaft_miac_m
from data.epcsaft_properties import molality_to_molefraction

matplotlib.use("Agg")
import matplotlib.pyplot as plt


T_REF = 298.15
P_REF = 1.0e5


SOLVENT_PROPS = {
    "Methanol": {
        "MW": 32.04e-3,
        "m": 1.5255,
        "s": 3.2300,
        "e": 188.90,
        "e_assoc": 2899.5,
        "vol_a": 0.03518,
        "dielc": 33.05,
    },
    "Ethanol": {
        "MW": 46.068e-3,
        "m": 2.3827,
        "s": 3.1771,
        "e": 198.24,
        "e_assoc": 2653.4,
        "vol_a": 0.03238,
        "dielc": 24.88,
    },
}


ION_PROPS = {
    "2008": {"na_s": 2.4122, "na_e": 646.05, "br_s": 3.4573, "br_e": 60.22, "k_nabr": 1.0},
    "2014": {"na_s": 2.8232, "na_e": 230.00, "br_s": 3.0707, "br_e": 190.00, "k_nabr": 0.290},
    "2020": {"na_s": 2.8232, "na_e": 230.00, "br_s": 3.0707, "br_e": 190.00, "k_nabr": 0.290},
}


def _load_exp_data(solvent):
    if solvent == "Methanol":
        path = ROOT / "data" / "MIAC_m" / "methanol" / "methanol-NaBr.csv"
        m_key, g_key = "molality (kg/mol)", "gamma"
    elif solvent == "Ethanol":
        path = ROOT / "data" / "MIAC_m" / "ethanol" / "ethanol-NaBr.csv"
        m_key, g_key = "molality (kg/mol)", "gamma"
    else:
        raise ValueError(f"Unsupported solvent: {solvent}")

    if not path.exists():
        raise FileNotFoundError(f"Experimental data file not found: {path}")

    molal = []
    gamma = []
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {m_key, g_key}
        if not required.issubset(set(reader.fieldnames or [])):
            raise ValueError(f"Missing columns {sorted(required)} in {path}.")
        for row in reader:
            m = float(row[m_key])
            g = float(row[g_key])
            molal.append(m)
            gamma.append(g)

    molal = np.asarray(molal, dtype=float)
    gamma = np.asarray(gamma, dtype=float)
    species = ["Na+", "Br-", solvent]
    x_salt = np.empty_like(molal)
    for i, m in enumerate(molal):
        x = molality_to_molefraction(float(m), species=species)
        x_salt[i] = x[0] + x[1]

    mask = (x_salt >= 0.0) & (x_salt <= 0.05)
    if not np.any(mask):
        raise ValueError(f"No {solvent} experimental points in 0<=x_salt<=0.05.")

    order = np.argsort(x_salt[mask])
    return x_salt[mask][order], gamma[mask][order]


def _params_for_version(version, solvent):
    if version not in ION_PROPS:
        raise ValueError(f"Unsupported version: {version}")
    if solvent not in SOLVENT_PROPS:
        raise ValueError(f"Unsupported solvent: {solvent}")

    ion = ION_PROPS[version]
    solv = SOLVENT_PROPS[solvent]

    k_ij = np.zeros((3, 3), dtype=float)
    k_ij[0, 1] = ion["k_nabr"]
    k_ij[1, 0] = ion["k_nabr"]
    # Requested setup: no ion-solvent kij for alcohol comparisons.
    k_ij[0, 2] = 0.0
    k_ij[2, 0] = 0.0
    k_ij[1, 2] = 0.0
    k_ij[2, 1] = 0.0

    params = {
        "m": np.asarray([1.0, 1.0, solv["m"]], dtype=float),
        "s": np.asarray([ion["na_s"], ion["br_s"], solv["s"]], dtype=float),
        "e": np.asarray([ion["na_e"], ion["br_e"], solv["e"]], dtype=float),
        "MW": np.asarray([22.98e-3, 79.904e-3, solv["MW"]], dtype=float),
        "e_assoc": np.asarray([0.0, 0.0, solv["e_assoc"]], dtype=float),
        "vol_a": np.asarray([0.0, 0.0, solv["vol_a"]], dtype=float),
        "z": np.asarray([1.0, -1.0, 0.0], dtype=float),
        "k_ij": k_ij,
        "dielc": np.asarray([8.0, 8.0, solv["dielc"]], dtype=float),
        "DH_model": 1,
        "debug": False,
    }
    if version == "2020":
        params["born_model"] = 1
        params["dielc_rule"] = 1
    else:
        params["born_model"] = 0
        params["dielc_rule"] = 0
    return params


def _m_from_x_salt(x_salt, mw_solvent):
    n_solvent = 1.0 / mw_solvent
    return x_salt * n_solvent / (2.0 * (1.0 - x_salt))


def _calc_curve(solvent, version, x_salt_grid):
    params = _params_for_version(version, solvent)
    species = ["Na+", "Br-", solvent]
    molal_grid = _m_from_x_salt(x_salt_grid, SOLVENT_PROPS[solvent]["MW"])
    gamma = np.empty_like(x_salt_grid, dtype=float)
    for i, m in enumerate(molal_grid):
        x = molality_to_molefraction(float(m), species=species)
        rho = pcsaft_den(T_REF, P_REF, x, params, phase="liq")
        gamma[i] = pcsaft_miac_m(T_REF, rho, x, params, species=species)["Na+Br-"]
    if not np.all(np.isfinite(gamma)):
        raise ValueError(f"Non-finite MIAC values for {solvent} {version}.")
    return gamma


def run_validation_miac_versions():
    out_dir = ROOT / "tests" / "fit_plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    style = {
        "2008": {"color": "gray", "linestyle": "--", "lw": 1.8},
        "2014": {"color": "black", "linestyle": "-", "lw": 2.0},
        "2020": {"color": "tab:blue", "linestyle": "-.", "lw": 1.8},
    }
    x_grids = {
        "Methanol": np.linspace(1e-6, 0.05, 101),
        "Ethanol": np.linspace(1e-6, 0.05, 101),
    }
    outputs = {
        "Methanol": out_dir / "validation_miac_versions_nabr_methanol.png",
        "Ethanol": out_dir / "validation_miac_versions_nabr_ethanol.png",
    }

    for solvent in ("Methanol", "Ethanol"):
        x_salt_grid = x_grids[solvent]
        x_exp, g_exp = _load_exp_data(solvent)
        curves = {v: _calc_curve(solvent, v, x_salt_grid) for v in ("2008", "2014", "2020")}

        fig, ax = plt.subplots(figsize=(7.2, 5.0))
        ax.scatter(x_exp, g_exp, color="black", marker="o", s=34, facecolors="none", label=f"NaBr data ({solvent})")
        for version in ("2008", "2014", "2020"):
            ax.plot(
                x_salt_grid,
                curves[version],
                color=style[version]["color"],
                linestyle=style[version]["linestyle"],
                linewidth=style[version]["lw"],
                label=f"NaBr {version}",
            )
        ax.set_xlim(0.0, 0.05)
        ax.set_ylim(0.0, 1.0)
        ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
        ax.set_xlabel(r"total NaBr mole fraction, $x_{\mathrm{NaBr}} = x_{\mathrm{Na^+}} + x_{\mathrm{Br^-}}$")
        ax.set_ylabel(r"mean ionic activity coefficient, $\gamma_{\pm}^{m}$")
        ax.set_title(f"NaBr in {solvent.lower()} at 298.15 K (2008/2014/2020)")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(outputs[solvent], dpi=220)
        plt.close(fig)

        if not outputs[solvent].exists():
            raise FileNotFoundError(f"Expected plot was not written: {outputs[solvent]}")


def test_validation_miac_versions():
    run_validation_miac_versions()


if __name__ == "__main__":
    run_validation_miac_versions()
