"""Figure 6a-style MIAC analysis for LiBr in ethanol.

This script compares four ePC-SAFT model variants against experimental MIAC data:
- ePC-SAFT advanced (2020)
- 2020 with Born disabled
- 2020 with constant dielectric rule
- ePC-SAFT revised baseline (2014 dataset)

Output is a single left-panel style plot with fixed axes:
- x: salt mole fraction (0 to 0.2)
- y: mean ionic activity coefficient (0 to 4)
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from data.epcsaft_properties import get_prop_dict
from pcsaft import pcsaft_den, pcsaft_miac_m

matplotlib.use("Agg")
import matplotlib.pyplot as plt


T_REF = 298.15
P_REF = 1.0e5
SALT = "LiBr"
SOLVENT = "ethanol"
CATION = "Li+"
ANION = "Br-"
PAIR_KEY = "Li+Br-"
SPECIES = [CATION, ANION, "Ethanol"]
MW_ETHANOL = 46.068e-3
SUM_NU = 2
SOURCE_COLOR_CYCLE = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown", "tab:pink", "tab:olive", "tab:cyan"]


def _read_csv_rows(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = [h.strip() for h in (reader.fieldnames or []) if h and h.strip()]
        rows = []
        for row in reader:
            clean: Dict[str, str] = {}
            for k, v in row.items():
                if not k:
                    continue
                ks = k.strip()
                if not ks:
                    continue
                clean[ks] = v.strip() if isinstance(v, str) else v
            rows.append(clean)
    return fields, rows


def _salt_mole_fraction_from_molality(molality: np.ndarray) -> np.ndarray:
    m = np.asarray(molality, dtype=float)
    n_solv = 1.0 / MW_ETHANOL
    denom = m + n_solv
    return np.where(denom > 0.0, m / denom, 0.0)


def _molality_for_salt_mole_fraction(x_salt: float) -> float:
    x = float(x_salt)
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        raise ValueError("Salt mole fraction target must be < 1.")
    n_solv = 1.0 / MW_ETHANOL
    return float((x * n_solv) / (1.0 - x))


def _molality_to_species_molefraction(molality: float) -> np.ndarray:
    m = max(float(molality), 0.0)
    n_solv = 1.0 / MW_ETHANOL
    n_cat = m
    n_an = m
    n_total = n_solv + n_cat + n_an
    if n_total <= 0.0:
        raise ValueError("Computed non-positive total moles while converting molality.")
    return np.asarray([n_cat / n_total, n_an / n_total, n_solv / n_total], dtype=float)


def _source_label(raw: str) -> str:
    s = str(raw or "").strip()
    return s if s else "Unspecified source"


def _color_for_source(source_label: str) -> str:
    if source_label == "Unspecified source":
        return "black"
    idx = sum(ord(ch) for ch in source_label) % len(SOURCE_COLOR_CYCLE)
    return SOURCE_COLOR_CYCLE[idx]


def _high_outlier_mask(values: np.ndarray) -> np.ndarray:
    """Return mask that removes only a single extreme high spike, if clearly separated."""
    y = np.asarray(values, dtype=float)
    finite = np.isfinite(y)
    positive = y > 0.0
    good_idx = np.where(finite & positive)[0]
    keep = finite.copy()

    # Conservative rule: only drop the single max point when it is both absolutely
    # large and strongly separated from the second-highest value.
    if good_idx.size >= 6:
        vals = y[good_idx]
        order = np.argsort(vals)
        v_max = float(vals[order[-1]])
        v_second = float(vals[order[-2]]) if vals.size >= 2 else 0.0
        if v_max > 5.0 and v_max > 3.0 * max(v_second, 1e-12):
            drop_idx = int(good_idx[order[-1]])
            keep[drop_idx] = False

    return keep


def _load_exp_data(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    fields, rows = _read_csv_rows(path)
    lookup = {f.lower(): f for f in fields}

    m_key = None
    for candidate in ("molality", "molality (kg/mol)", "m", "m (mol/kg)"):
        if candidate in lookup:
            m_key = lookup[candidate]
            break
    if m_key is None:
        raise ValueError(f"Missing molality column in {path}.")

    y_key = None
    for candidate in ("miac", "y"):
        if candidate in lookup:
            y_key = lookup[candidate]
            break
    if y_key is None:
        raise ValueError(f"Missing MIAC column in {path}. Tried: 'miac', 'y'.")

    x_key = lookup.get("mole_fraction")
    source_key = lookup.get("source")

    molality: List[float] = []
    x_salt: List[float] = []
    miac: List[float] = []
    sources: List[str] = []

    for row in rows:
        try:
            m_val = float(row.get(m_key, ""))
            y_val = float(row.get(y_key, ""))
        except (TypeError, ValueError):
            continue
        if not (math.isfinite(m_val) and math.isfinite(y_val)):
            continue

        if x_key is not None:
            try:
                x_val = float(row.get(x_key, ""))
                if not math.isfinite(x_val):
                    x_val = float(_salt_mole_fraction_from_molality(np.asarray([m_val], dtype=float))[0])
            except (TypeError, ValueError):
                x_val = float(_salt_mole_fraction_from_molality(np.asarray([m_val], dtype=float))[0])
        else:
            x_val = float(_salt_mole_fraction_from_molality(np.asarray([m_val], dtype=float))[0])

        src_raw = row.get(source_key, "") if source_key else ""
        sources.append(_source_label(src_raw))
        molality.append(m_val)
        x_salt.append(x_val)
        miac.append(y_val)

    if not molality:
        raise ValueError(f"No usable rows in {path}.")

    m_arr = np.asarray(molality, dtype=float)
    x_arr = np.asarray(x_salt, dtype=float)
    y_arr = np.asarray(miac, dtype=float)
    s_arr = np.asarray(sources, dtype=object)

    order = np.argsort(m_arr)
    return m_arr[order], x_arr[order], y_arr[order], s_arr[order]


def _build_params(dataset_name: str, user_options: dict | None = None) -> Dict[str, object]:
    x_ref = _molality_to_species_molefraction(1e-8)
    return get_prop_dict(dataset_name, SPECIES, x_ref, T_REF, user_options=user_options)


def _calc_curve(m_grid: np.ndarray, dataset_name: str, user_options: dict | None = None) -> np.ndarray:
    params = _build_params(dataset_name, user_options=user_options)
    gamma_m = np.empty_like(m_grid, dtype=float)

    for idx, m in enumerate(m_grid):
        m_eval = float(m) if m > 0.0 else 1e-12
        x = _molality_to_species_molefraction(m_eval)
        rho = pcsaft_den(T_REF, P_REF, x, params, phase="liq")
        gamma_m[idx] = pcsaft_miac_m(T_REF, rho, x, params, species=SPECIES)[PAIR_KEY]

    if not np.all(np.isfinite(gamma_m)):
        raise ValueError(f"Non-finite MIAC_m curve for dataset={dataset_name}, options={user_options}.")

    factor = 1.0 + MW_ETHANOL * m_grid * SUM_NU
    gamma = gamma_m * factor
    if not np.all(np.isfinite(gamma)):
        raise ValueError(f"Non-finite MIAC curve for dataset={dataset_name}, options={user_options}.")
    return gamma


def _curve_configs() -> List[Dict[str, object]]:
    return [
        {
            "label": "ePC-SAFT advanced (2020)",
            "dataset": "bulow_2020",
            "color": "green",
            "lw": 2.2,
            "linestyle": "-",
            "zorder": 3,
            "user_options": {},
        },
        {
            "label": "2020 w/o Born (born_model=0)",
            "dataset": "bulow_2020",
            "color": "gray",
            "lw": 2.0,
            "linestyle": "-",
            "zorder": 2,
            "user_options": {"elec_model": {"include_born_model": False}},
        },
        {
            "label": "2020 with constant dielectric rule",
            "dataset": "bulow_2020",
            "color": "orange",
            "lw": 2.2,
            "linestyle": "--",
            "zorder": 5,
            "user_options": {"elec_model": {"rel_perm": {"rule": "constant"}}},
        },
        {
            "label": "ePC-SAFT revised (2014)",
            "dataset": "held_2014",
            "color": "black",
            "lw": 2.0,
            "linestyle": "-",
            "zorder": 4,
            "user_options": {},
        },
    ]


def run_analysis(
    data_path: Path,
    output_path: Path,
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    grid_points: int,
    max_molality: float | None,
) -> Path:
    m_exp_raw, x_exp_raw, y_exp_raw, source_exp_raw = _load_exp_data(data_path)

    keep = _high_outlier_mask(y_exp_raw)
    removed = int(np.count_nonzero(~keep))
    if np.any(keep):
        m_exp = m_exp_raw[keep]
        x_exp = x_exp_raw[keep]
        y_exp = y_exp_raw[keep]
        source_exp = source_exp_raw[keep]
    else:
        m_exp = m_exp_raw
        x_exp = x_exp_raw
        y_exp = y_exp_raw
        source_exp = source_exp_raw
        removed = 0

    if removed > 0:
        print(f"[outlier-filter] figure6a removed {removed} high outlier experimental point(s).")

    m_upper = float(np.max(m_exp)) if max_molality is None else float(max_molality)
    if m_upper <= 0.0:
        m_upper = float(np.max(m_exp[m_exp >= 0.0]))

    # Ensure model curves span the full requested mole-fraction axis (e.g., x=0.2).
    m_upper = max(m_upper, _molality_for_salt_mole_fraction(x_max))

    m_grid = np.linspace(0.0, m_upper, int(grid_points))
    x_grid = _salt_mole_fraction_from_molality(m_grid)

    curves: List[Tuple[Dict[str, object], np.ndarray]] = []
    for cfg in _curve_configs():
        y_curve = _calc_curve(m_grid, str(cfg["dataset"]), user_options=dict(cfg["user_options"]))
        curves.append((cfg, y_curve))

    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    unique_sources = []
    for src in source_exp.tolist():
        if src not in unique_sources:
            unique_sources.append(src)

    if len(unique_sources) <= 1 and unique_sources[0] == "Unspecified source":
        ax.scatter(
            x_exp,
            y_exp,
            color="black",
            marker="o",
            s=38,
            facecolors="none",
            linewidths=1.1,
            label="Experimental data",
            zorder=7,
        )
    else:
        for src in unique_sources:
            mask = source_exp == src
            color = _color_for_source(str(src))
            lbl = f"Experimental data - {src}" if src != "Unspecified source" else "Experimental data - Unspecified source"
            ax.scatter(
                x_exp[mask],
                y_exp[mask],
                color=color,
                marker="o",
                s=40,
                facecolors="none",
                linewidths=1.15,
                label=lbl,
                zorder=7,
            )

    for cfg, y_curve in curves:
        ax.plot(
            x_grid,
            y_curve,
            color=str(cfg["color"]),
            linewidth=float(cfg["lw"]),
            linestyle=str(cfg.get("linestyle", "-")),
            zorder=float(cfg.get("zorder", 3)),
            label=str(cfg["label"]),
        )

    ax.set_xlim(float(x_min), float(x_max))
    ax.set_ylim(float(y_min), float(y_max))
    ax.set_xlabel(r"salt mole fraction, $x_{salt}$")
    ax.set_ylabel(r"mean ionic activity coefficient, $\gamma_{\pm}^{*}$")
    ax.set_title("LiBr in ethanol at 298.15 K and 1 bar (Figure 6a-style)")
    ax.grid(True, alpha=0.3, color="0.7")
    ax.tick_params(colors="black")
    for spine in ax.spines.values():
        spine.set_color("black")
        spine.set_linewidth(1.0)

    legend = ax.legend(fontsize=8)
    frame = legend.get_frame()
    frame.set_facecolor("white")
    frame.set_edgecolor("black")
    frame.set_alpha(1.0)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=220)
    plt.close(fig)

    if not output_path.exists():
        raise FileNotFoundError(f"Expected plot was not written: {output_path}")

    print(f"Loaded rows: {len(m_exp)} from {data_path}")
    src_counts = {s: int(np.sum(source_exp == s)) for s in np.unique(source_exp)}
    print(f"Data sources: {src_counts}")
    print(f"Molality grid points: {len(m_grid)} (0 to {m_upper:.6g} mol/kg)")
    print("Curves:")
    for cfg, _ in curves:
        print(f"- {cfg['label']} [{cfg['dataset']}]")
    print(f"Wrote: {output_path}")
    return output_path


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Figure 6a-style LiBr/ethanol MIAC analysis")
    parser.add_argument(
        "--data",
        type=Path,
        default=REPO_ROOT / "data" / "MIAC" / "ethanol" / "ethanol-LiBr.csv",
        help="Input CSV with columns including molality, mole_fraction, miac.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=REPO_ROOT / "scripts" / "2020_model_analysis" / "output" / "figure6a_libr_ethanol_2020_variants.png",
        help="Output PNG path.",
    )
    parser.add_argument("--x-min", type=float, default=0.0)
    parser.add_argument("--x-max", type=float, default=0.2)
    parser.add_argument("--y-min", type=float, default=0.0)
    parser.add_argument("--y-max", type=float, default=4.0)
    parser.add_argument("--grid-points", type=int, default=1201)
    parser.add_argument(
        "--max-molality",
        type=float,
        default=None,
        help="Optional upper molality bound for dense curve generation. Default: max experimental molality.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    run_analysis(
        data_path=Path(args.data),
        output_path=Path(args.out),
        x_min=float(args.x_min),
        x_max=float(args.x_max),
        y_min=float(args.y_min),
        y_max=float(args.y_max),
        grid_points=int(args.grid_points),
        max_molality=None if args.max_molality is None else float(args.max_molality),
    )


if __name__ == "__main__":
    main()
