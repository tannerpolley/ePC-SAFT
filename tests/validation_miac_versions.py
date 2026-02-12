"""Alcohol MIAC comparisons across 2008/2014/2020 model versions."""

import csv
import math
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
AXIS_LABEL_SIZE = 13
AXIS_TICK_SIZE = 11


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
    "2008": {
        "Li+": {"s": 1.8177, "e": 2697.28, "MW": 6.941e-3, "z": 1.0},
        "Na+": {"s": 2.4122, "e": 646.05, "MW": 22.98e-3, "z": 1.0},
        "K+": {"s": 2.9698, "e": 271.05, "MW": 39.0983e-3, "z": 1.0},
        "Cl-": {"s": 3.0575, "e": 47.29, "MW": 35.453e-3, "z": -1.0},
        "Br-": {"s": 3.4573, "e": 60.22, "MW": 79.904e-3, "z": -1.0},
        "I-": {"s": 3.9319, "e": 80.43, "MW": 126.90447e-3, "z": -1.0},
    },
    "2014": {
        "Li+": {"s": 2.8449, "e": 360.00, "MW": 6.941e-3, "z": 1.0},
        "Na+": {"s": 2.8232, "e": 230.00, "MW": 22.98e-3, "z": 1.0},
        "K+": {"s": 3.3417, "e": 200.00, "MW": 39.0983e-3, "z": 1.0},
        "Cl-": {"s": 2.7560, "e": 170.00, "MW": 35.453e-3, "z": -1.0},
        "Br-": {"s": 3.0707, "e": 190.00, "MW": 79.904e-3, "z": -1.0},
        "I-": {"s": 3.6672, "e": 200.00, "MW": 126.90447e-3, "z": -1.0},
    },
    "2020": {
        "Li+": {"s": 2.8449, "e": 360.00, "MW": 6.941e-3, "z": 1.0},
        "Na+": {"s": 2.8232, "e": 230.00, "MW": 22.98e-3, "z": 1.0},
        "K+": {"s": 3.3417, "e": 200.00, "MW": 39.0983e-3, "z": 1.0},
        "Cl-": {"s": 2.7560, "e": 170.00, "MW": 35.453e-3, "z": -1.0},
        "Br-": {"s": 3.0707, "e": 190.00, "MW": 79.904e-3, "z": -1.0},
        "I-": {"s": 3.6672, "e": 200.00, "MW": 126.90447e-3, "z": -1.0},
    },
}

SALT_SPECS = {
    "LiBr": {"cation": "Li+", "anion": "Br-"},
    "LiCl": {"cation": "Li+", "anion": "Cl-"},
    "NaBr": {"cation": "Na+", "anion": "Br-"},
    "NaCl": {"cation": "Na+", "anion": "Cl-"},
    "NaI": {"cation": "Na+", "anion": "I-"},
    "KBr": {"cation": "K+", "anion": "Br-"},
    "KI": {"cation": "K+", "anion": "I-"},
}

ION_PAIR_K = {
    "2008": {salt: 1.0 for salt in SALT_SPECS},
    "2014": {
        "LiBr": 0.591,
        "LiCl": 0.669,
        "NaBr": 0.290,
        "NaCl": 0.317,
        "NaI": 0.018,
        "KBr": -0.102,
        "KI": -0.312,
    },
    "2020": {
        "LiBr": 0.591,
        "LiCl": 0.669,
        "NaBr": 0.290,
        "NaCl": 0.317,
        "NaI": 0.018,
        "KBr": -0.102,
        "KI": -0.312,
    },
}

VERSIONS = ("2008", "2014", "2020_m0", "2020_m1", "2020_m2")

STYLE = {
    "2008": {"color": "gray", "linestyle": "--", "lw": 2.0},
    "2014": {"color": "black", "linestyle": "-", "lw": 2.0},
    "2020_m0": {"color": "tab:blue", "linestyle": "-.", "lw": 2.0},
    "2020_m1": {"color": "tab:red", "linestyle": ":", "lw": 2.2},
    "2020_m2": {"color": "tab:green", "linestyle": (0, (3, 1, 1, 1)), "lw": 2.0},
}


def _discover_combos():
    combos = []
    for solvent in ("Methanol", "Ethanol"):
        solvent_lower = solvent.lower()
        prefix = f"{solvent_lower}-"
        data_dir = ROOT / "data" / "MIAC_m" / solvent_lower
        if not data_dir.exists():
            continue

        for path in sorted(data_dir.glob(f"{prefix}*.csv")):
            salt = path.stem.replace(prefix, "", 1)
            if salt not in SALT_SPECS:
                raise ValueError(
                    f"Unsupported salt '{salt}' discovered in {path}. "
                    f"Add it to SALT_SPECS/ION_PROPS/ION_PAIR_K."
                )
            combos.append(
                {
                    "salt": salt,
                    "solvent": solvent,
                    "data_path": path,
                    "output": ROOT / "tests" / "fit_plots" / f"validation_miac_versions_{salt.lower()}_{solvent_lower}.png",
                }
            )

    if not combos:
        raise FileNotFoundError("No alcohol MIAC CSV files found in data/MIAC_m/methanol or data/MIAC_m/ethanol.")
    return combos


def _species_for_combo(salt, solvent):
    salt_spec = SALT_SPECS[salt]
    return [salt_spec["cation"], salt_spec["anion"], solvent]


def _pair_key(salt):
    salt_spec = SALT_SPECS[salt]
    return f"{salt_spec['cation']}{salt_spec['anion']}"


def _load_exp_data(combo):
    path = combo["data_path"]
    salt = combo["salt"]
    molal = []
    miac = []

    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])

        m_key = None
        for candidate in ("molality (kg/mol)", f"m_{salt}", "molality", "m"):
            if candidate in fields:
                m_key = candidate
                break
        if m_key is None:
            raise ValueError(
                f"Missing molality column in {path}. Tried: 'molality (kg/mol)', 'm_{salt}', 'molality', 'm'."
            )

        y_key = None
        for candidate in ("miac_m", "gamma"):
            if candidate in fields:
                y_key = candidate
                break
        if y_key is None:
            raise ValueError(f"Missing MIAC column in {path}. Tried: 'miac_m', 'gamma'.")

        for row in reader:
            molal.append(float(row[m_key]))
            miac.append(float(row[y_key]))

    molal = np.asarray(molal, dtype=float)
    miac = np.asarray(miac, dtype=float)
    if not np.all(np.isfinite(molal)) or not np.all(np.isfinite(miac)):
        raise ValueError(f"Non-finite experimental values in {path}.")
    if molal.size == 0:
        raise ValueError(f"No data rows in {path}.")

    order = np.argsort(molal)
    return molal[order], miac[order]


def _params_for_version(version, salt, solvent):
    base = "2020" if version.startswith("2020_") else version
    cation = SALT_SPECS[salt]["cation"]
    anion = SALT_SPECS[salt]["anion"]
    ion = ION_PROPS[base]
    solv = SOLVENT_PROPS[solvent]

    k_ij = np.zeros((3, 3), dtype=float)
    kij_pair = ION_PAIR_K[base][salt]
    k_ij[0, 1] = kij_pair
    k_ij[1, 0] = kij_pair

    params = {
        "m": np.asarray([1.0, 1.0, solv["m"]], dtype=float),
        "s": np.asarray([ion[cation]["s"], ion[anion]["s"], solv["s"]], dtype=float),
        "e": np.asarray([ion[cation]["e"], ion[anion]["e"], solv["e"]], dtype=float),
        "MW": np.asarray([ion[cation]["MW"], ion[anion]["MW"], solv["MW"]], dtype=float),
        "e_assoc": np.asarray([0.0, 0.0, solv["e_assoc"]], dtype=float),
        "vol_a": np.asarray([0.0, 0.0, solv["vol_a"]], dtype=float),
        "z": np.asarray([ion[cation]["z"], ion[anion]["z"], 0.0], dtype=float),
        "k_ij": k_ij,
        "dielc": np.asarray([8.0, 8.0, solv["dielc"]], dtype=float),
        "DH_model": 1,
        "debug": False,
    }
    if version.startswith("2020_"):
        params["born_model"] = 1
        params["dielc_rule"] = 1
        params["born_diff_mode"] = int(version.split("_m", 1)[1])
    else:
        params["born_model"] = 0
        params["dielc_rule"] = 0
        params["born_diff_mode"] = 0
    return params


def _calc_curve(combo, version, molal_grid):
    salt = combo["salt"]
    solvent = combo["solvent"]
    params = _params_for_version(version, salt, solvent)
    species = _species_for_combo(salt, solvent)
    pair_key = _pair_key(salt)

    gamma = np.empty_like(molal_grid, dtype=float)
    for i, m in enumerate(molal_grid):
        # Use a tiny positive molality for the left endpoint to keep the density/fugacity solve stable.
        m_eval = float(m) if m > 0.0 else 1e-12
        x = molality_to_molefraction(m_eval, species=species)
        rho = pcsaft_den(T_REF, P_REF, x, params, phase="liq")
        gamma[i] = pcsaft_miac_m(T_REF, rho, x, params, species=species)[pair_key]
    if not np.all(np.isfinite(gamma)):
        raise ValueError(f"Non-finite MIAC values for {salt}/{solvent} {version}.")
    return gamma


def run_validation_miac_versions():
    out_dir = ROOT / "tests" / "fit_plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    combos = _discover_combos()
    generated = []

    for combo in combos:
        salt = combo["salt"]
        solvent = combo["solvent"]
        molal_exp, miac_exp = _load_exp_data(combo)

        xmax = float(max(1, int(math.ceil(float(np.max(molal_exp))))))
        ymax = float(max(1, int(math.ceil(float(np.max(miac_exp))))))
        molal_grid = np.linspace(0.0, xmax, 701)
        curves = {version: _calc_curve(combo, version, molal_grid) for version in VERSIONS}

        fig, ax = plt.subplots(figsize=(7.2, 5.0))
        ax.scatter(
            molal_exp,
            miac_exp,
            color="black",
            marker="o",
            s=34,
            facecolors="none",
            label=f"{salt} data ({solvent})",
        )
        for version in VERSIONS:
            label = f"{salt} 2020 born_diff_mode={version[-1]}" if version.startswith("2020_") else f"{salt} {version}"
            ax.plot(
                molal_grid,
                curves[version],
                color=STYLE[version]["color"],
                linestyle=STYLE[version]["linestyle"],
                linewidth=STYLE[version]["lw"],
                label=label,
            )

        ax.set_xlim(0.0, xmax)
        ax.set_ylim(0.0, ymax)
        ax.set_xlabel(r"molality, $m$ / mol kg$^{-1}$", fontsize=AXIS_LABEL_SIZE)
        ax.set_ylabel(r"mean ionic activity coefficient, $\gamma_{\pm}^{m}$", fontsize=AXIS_LABEL_SIZE)
        ax.tick_params(axis="both", labelsize=AXIS_TICK_SIZE)
        ax.set_title(f"{salt} in {solvent.lower()} at 298.15 K (2008/2014/2020 born-diff modes)")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(combo["output"], dpi=220)
        plt.close(fig)

        if not combo["output"].exists():
            raise FileNotFoundError(f"Expected plot was not written: {combo['output']}")
        generated.append(combo["output"])

    print("Generated validation plots:")
    for path in generated:
        print(f"- {path}")


def test_validation_miac_versions():
    run_validation_miac_versions()


if __name__ == "__main__":
    run_validation_miac_versions()
