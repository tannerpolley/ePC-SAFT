"""MIAC comparisons across 2005/2008/2014 with selected 2020/2025 model variants."""

import csv
import math
import os
import sys
from pathlib import Path

import matplotlib
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from pcsaft import pcsaft_den, pcsaft_miac_m
from data.epcsaft_properties import get_prop_dict, molality_to_molefraction

matplotlib.use("Agg")
import matplotlib.pyplot as plt


T_REF = 298.15
P_REF = 1.0e5
AXIS_LABEL_SIZE = 13
AXIS_TICK_SIZE = 11


SALT_SPECS = {
    "LiBr": {"cation": "Li+", "anion": "Br-"},
    "LiCl": {"cation": "Li+", "anion": "Cl-"},
    "LiI": {"cation": "Li+", "anion": "I-"},
    "NaBr": {"cation": "Na+", "anion": "Br-"},
    "NaCl": {"cation": "Na+", "anion": "Cl-"},
    "NaI": {"cation": "Na+", "anion": "I-"},
    "KCl": {"cation": "K+", "anion": "Cl-"},
    "KBr": {"cation": "K+", "anion": "Br-"},
    "KI": {"cation": "K+", "anion": "I-"},
}

VERSIONS = ("2005", "2008", "2014", "2020_normal", "2020_no_dielc_dep", "2020_sum_off", "2025_num")

STYLE = {
    "2005": {"color": "black", "linestyle": "--", "lw": 1.5},
    "2008": {"color": "dimgray", "linestyle": "--", "lw": 1.5},
    "2014": {"color": "silver", "linestyle": "--", "lw": 1.5},
    "2020_normal": {"color": "tab:purple", "linestyle": "--", "lw": 1.5},
    "2020_no_dielc_dep": {"color": "tab:blue", "linestyle": "--", "lw": 1.5},
    "2020_sum_off": {"color": "tab:orange", "linestyle": "--", "lw": 1.5},
    "2025_num": {"color": "tab:red", "linestyle": "--", "lw": 1.5},
}
DISPLAY_LABEL = {
    "2005": "2005",
    "2008": "2008",
    "2014": "2014",
    "2020_normal": "2020 (normal)",
    "2020_no_dielc_dep": "2020 (no dielc dep)",
    "2020_sum_off": "2020 (sum term off)",
    "2025_num": "2025 numeric",
}


def _requested_scope():
    solvent = os.getenv("MIAC_SOLVENT", "all").strip().lower()
    if solvent not in {"all", "water", "methanol", "ethanol"}:
        raise ValueError("MIAC_SOLVENT must be one of: all, water, methanol, ethanol.")
    salt = os.getenv("MIAC_SALT", "").strip()
    return solvent, salt


def _canonical_salt_token(salt):
    salts = {"LiBr", "LiCl", "LiI", "NaBr", "NaCl", "NaI", "KCl", "KBr", "KI"}
    for item in salts:
        if item.lower() == salt.lower():
            return item
    return salt


def _discover_combos():
    selected_solvent, selected_salt = _requested_scope()
    combos = []
    for solvent in ("Water", "Methanol", "Ethanol"):
        solvent_lower = solvent.lower()
        data_dir = REPO_ROOT / "data" / "MIAC_m" / solvent_lower
        if not data_dir.exists():
            continue
        if selected_solvent != "all" and selected_solvent != solvent_lower:
            continue

        if solvent == "Water":
            candidates = sorted(data_dir.glob("*.csv"))
        else:
            prefix = f"{solvent_lower}-"
            candidates = sorted(data_dir.glob(f"{prefix}*.csv"))

        for path in candidates:
            if solvent == "Water":
                salt = path.stem
            else:
                salt = path.stem.replace(prefix, "", 1)
            if salt not in SALT_SPECS:
                raise ValueError(
                    f"Unsupported salt '{salt}' discovered in {path}. "
                    f"Add it to SALT_SPECS/ION_PROPS/ION_PAIR_K."
                )
            if selected_salt and salt.lower() != selected_salt.lower():
                continue
            salt_token = _canonical_salt_token(salt)
            output_dir = path.parent / "plot_fits"
            combos.append(
                {
                    "salt": salt,
                    "solvent": solvent,
                    "data_path": path,
                    "output_dir": output_dir,
                    "output": output_dir / f"maic_m_{solvent_lower}_{salt_token}_fit.png",
                }
            )

    if not combos:
        scope = f"solvent={selected_solvent}, salt={selected_salt or '<all>'}"
        raise FileNotFoundError(f"No MIAC CSV files found for {scope}.")
    return combos


def _species_for_combo(salt, solvent):
    salt_spec = SALT_SPECS[salt]
    if solvent == "Water":
        water_species = "H2O-2B-Li" if salt.startswith("Li") else "H2O-2B-NaCl"
        return [salt_spec["cation"], salt_spec["anion"], water_species]
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


def _params_for_version(version, species):
    if version == "2005":
        options = {"elec_model": {"preset": "2005", "dielc_rule": "constant"}, "debug": False}
    elif version == "2008":
        options = {"elec_model": {"preset": "2008", "dielc_rule": "constant"}, "debug": False}
    elif version == "2014":
        options = {"elec_model": {"preset": "2014_s2", "dielc_rule": "constant"}, "debug": False}
    elif version == "2020_normal":
        options = {
            "elec_model": {
                "preset": "2020",
                "dielc_rule": 3,
                "born_diff_model": "analytic",
                "born_diff_options": {
                    "include_dielc_conc_dep": True,
                    "include_sum_term": True,
                },
            },
            "debug": False,
        }
    elif version == "2020_no_dielc_dep":
        options = {
            "elec_model": {
                "preset": "2020",
                "dielc_rule": 3,
                "born_diff_model": "analytic",
                "born_diff_options": {"include_dielc_conc_dep": False},
            },
            "debug": False,
        }
    elif version == "2020_sum_off":
        options = {
            "elec_model": {
                "preset": "2020",
                "dielc_rule": 3,
                "born_diff_model": None,
                "born_diff_options": {"include_sum_term": False},
            },
            "debug": False,
        }
    elif version == "2025_num":
        options = {
            "elec_model": {
                "preset": "2025",
                "ssm_ds": True,
                "dielc_rule": "empirical",
                "dielc_diff_mode": "numeric",
                "eps_r_bulk": "solvent",
                "born_radius_model": 5,
                "born_diff_model": "numeric",
                "born_diff_options": {
                    "include_dielc_conc_dep": True,
                    "include_delta_d_i_conc_dep": True,
                },
            },
            "debug": False,
        }
    else:
        raise ValueError(f"Unsupported version '{version}'.")

    x_ref = molality_to_molefraction(1e-8, species=species)
    return get_prop_dict(species, x_ref, T_REF, user_options=options)


def _calc_curve(combo, version, molal_grid):
    salt = combo["salt"]
    solvent = combo["solvent"]
    species = _species_for_combo(salt, solvent)
    params = _params_for_version(version, species)
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


def run_validation_miac_versions2():
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
        combo["output_dir"].mkdir(parents=True, exist_ok=True)

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
            label = f"{salt} {DISPLAY_LABEL[version]}"
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
        ax.set_title(f"{salt} in {solvent.lower()} at 298.15 K (2005/2008/2014 + 2020/2025 variants)")
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


def test_validation_miac_versions2():
    run_validation_miac_versions2()


if __name__ == "__main__":
    run_validation_miac_versions2()


# Backward-compatible alias while this script transitions to the v2 name.
run_validation_miac_versions = run_validation_miac_versions2
