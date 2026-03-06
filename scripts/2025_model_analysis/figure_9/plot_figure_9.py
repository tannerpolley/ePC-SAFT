from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import _common as common

OUTPUT = Path(__file__).with_name("figure_9.png")
PANELS = [
    ("a)", "NaBr", "water-methanol", 3.0),
    ("b)", "NaBr", "water-ethanol", 4.0),
    ("c)", "NaCl", "water-methanol", 2.0),
    ("d)", "NaCl", "water-ethanol", 2.0),
]
TARGETS = [
    (0.8, common.ORGANIC_COLOR, "^", "0.25"),
    (0.4, common.GREEN_COLOR, "s", "black"),
]


def _data_path(solvent_system: str, salt: str) -> Path:
    return common.REPO_ROOT / "data" / "MIAC" / solvent_system / f"{solvent_system}-{salt}.csv"


def _weight_to_mole_comp(solvent_system: str, weight_comp: dict[str, float]) -> dict[str, float]:
    solvents = [s for s in solvent_system.split("-") if s]
    numerators = {s: weight_comp[s] / common.SOLVENT_MW[s] for s in solvents}
    total = sum(numerators.values())
    return {s: numerators[s] / total for s in solvents}


def _read_weight_fraction_dataset(path: Path, solvent_system: str):
    fields, rows = common.read_csv_rows(path)
    lookup = {field.lower(): field for field in fields}
    molality_key = next((lookup[c] for c in ("molality", "m") if c in lookup), None)
    gamma_key = next((lookup[c] for c in ("miac_m", "gamma") if c in lookup), None)
    if molality_key is None or gamma_key is None:
        raise ValueError(f"Missing columns in {path}")
    organic = [s for s in solvent_system.split("-") if s != "water"][0]
    org_key = lookup.get(f"x_{organic}".lower())
    if org_key is None:
        org_key = lookup.get("x_methanol") if organic == "methanol" else lookup.get("x_ethanol")
    water_key = lookup.get("x_h2o") or lookup.get("x_water")
    out = []
    for row in rows:
        m = common.parse_float(row.get(molality_key))
        y = common.parse_float(row.get(gamma_key))
        w_org = common.parse_float(row.get(org_key)) if org_key else None
        w_water = common.parse_float(row.get(water_key)) if water_key else None
        if m is None or y is None or w_org is None:
            continue
        if w_water is None:
            w_water = 1.0 - w_org
        weight_comp = common.normalized_comp(solvent_system, {"water": w_water, organic: w_org})
        mole_comp = _weight_to_mole_comp(solvent_system, weight_comp)
        out.append(
            {
                "molality": m,
                "miac_m": y,
                "w_org": weight_comp[organic],
                "weight_signature": tuple((s, round(weight_comp[s], 6)) for s in weight_comp),
                "mole_comp": mole_comp,
            }
        )
    out.sort(key=lambda item: item["molality"])
    return out


def _group_by_weight(entries):
    grouped = {}
    for entry in entries:
        grouped.setdefault(entry["weight_signature"], []).append(entry)
    for rows in grouped.values():
        rows.sort(key=lambda item: item["molality"])
    return grouped


def _closest_group(entries, target_w_org: float):
    grouped = _group_by_weight(entries)
    candidates = []
    for rows in grouped.values():
        w_org = float(rows[0]["w_org"])
        candidates.append((abs(w_org - target_w_org), rows, w_org, rows[0]["mole_comp"]))
    candidates.sort(key=lambda item: item[0])
    return candidates[0]


def _plot_panel(ax, label, salt, solvent_system, m_max):
    path = _data_path(solvent_system, salt)
    entries = _read_weight_fraction_dataset(path, solvent_system) if path.exists() else []
    used_signatures = set()
    for target_w, marker_color, marker, line_color in TARGETS:
        comp_model = common.target_weight_fraction_to_comp(solvent_system, target_w)
        m_grid, y_model = common.mean_ionic_activity_curve("figiel_2025", salt, solvent_system, comp_model, m_max, points=600)
        ax.plot(m_grid, y_model, color=line_color, linewidth=1.3)
        if entries:
            _, rows, w_actual, _ = _closest_group(entries, target_w)
            signature = rows[0]["weight_signature"]
            if signature not in used_signatures:
                used_signatures.add(signature)
                ax.scatter([r["molality"] for r in rows], [r["miac_m"] for r in rows], s=24, marker=marker, facecolor=marker_color, edgecolor=marker_color, linewidth=0.8, label=common.safe_label_for_weight_fraction(target_w) if label == "a)" else None)
        else:
            print(f"[figure_9] missing experimental data: {path}")
    organic = [s for s in solvent_system.split("-") if s != "water"][0]
    common.panel_label(ax, label)
    ax.set_xlim(0.0, m_max)
    ax.set_ylim(0.0, 1.125)
    ax.set_title(f"{salt} in {organic}", fontsize=10)
    ax.set_xlabel(r"$\bar{m}_{salt}$ / mol kg$^{-1}$")
    ax.set_ylabel(r"$\gamma_{\pm}^{m,*}$ / -")


def main() -> None:
    common.configure_style()
    fig, axes = plt.subplots(2, 2, figsize=(6.6, 5.8))
    for ax, cfg in zip(axes.flat, PANELS):
        _plot_panel(ax, *cfg)
    handles = [
        plt.Line2D([0], [0], marker="^", linestyle="None", color=common.ORGANIC_COLOR, markerfacecolor=common.ORGANIC_COLOR, markeredgecolor=common.ORGANIC_COLOR, label=r"$w_{org}^{salt-free}=0.8$"),
        plt.Line2D([0], [0], marker="s", linestyle="None", color=common.GREEN_COLOR, markerfacecolor=common.GREEN_COLOR, markeredgecolor=common.GREEN_COLOR, label=r"$w_{org}^{salt-free}=0.4$"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.02), fontsize=9)
    common.save_figure(fig, OUTPUT)
    plt.close(fig)


if __name__ == "__main__":
    main()
