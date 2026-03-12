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
    (0.8, common.ORGANIC_COLOR, "^", common.ORGANIC_COLOR, "-", r"Model, $w_{org}^{salt-free}=0.8$"),
    (0.4, common.GREEN_COLOR, "s", common.GREEN_COLOR, "-", r"Model, $w_{org}^{salt-free}=0.4$"),
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
    w_org_key = lookup.get(f"w_{organic}_salt_free".lower()) or lookup.get(f"w_{organic}".lower())
    w_water_key = lookup.get("w_h2o_salt_free") or lookup.get("w_water_salt_free") or lookup.get("w_h2o") or lookup.get("w_water")
    org_key = lookup.get(f"x_{organic}".lower())
    if org_key is None:
        org_key = lookup.get("x_methanol") if organic == "methanol" else lookup.get("x_ethanol")
    water_key = lookup.get("x_h2o") or lookup.get("x_water")
    out = []
    for row in rows:
        m = common.parse_float(row.get(molality_key))
        y = common.parse_float(row.get(gamma_key))
        w_org = common.parse_float(row.get(w_org_key)) if w_org_key else None
        w_water = common.parse_float(row.get(w_water_key)) if w_water_key else None
        if w_org is None:
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
        candidates.append((abs(w_org - target_w_org), rows))
    candidates.sort(key=lambda item: item[0])
    return candidates[0][1]


def _plot_panel(ax, label: str, salt: str, solvent_system: str, m_max: float) -> None:
    path = _data_path(solvent_system, salt)
    entries = _read_weight_fraction_dataset(path, solvent_system) if path.exists() else []

    for target_w, marker_color, marker, line_color, line_style, line_label in TARGETS:
        comp_model = common.target_weight_fraction_to_comp(solvent_system, target_w)
        m_grid, y_model = common.mean_ionic_activity_curve("figiel_2025", salt, solvent_system, comp_model, m_max, points=600)
        ax.plot(m_grid, y_model, color=line_color, linestyle=line_style, linewidth=2.2, zorder=5)

        if entries:
            rows = _closest_group(entries, target_w)
            data_label = common.safe_label_for_weight_fraction(target_w) if label == "a)" else None
            ax.scatter(
                [r["molality"] for r in rows],
                [r["miac_m"] for r in rows],
                s=24,
                marker=marker,
                facecolor="none",
                edgecolor=marker_color,
                linewidth=1.0,
                zorder=6,
                label=data_label,
            )
        else:
            print(f"[figure_9] missing experimental data: {path}")

        if label == "a)":
            ax.plot([], [], color=line_color, linestyle=line_style, linewidth=2.2, label=line_label)

    organic = [s for s in solvent_system.split("-") if s != "water"][0]
    common.panel_label(ax, label)
    ax.set_xlim(0.0, m_max)
    ax.set_ylim(0.0, 1.125)
    ax.set_title(f"{salt} in {organic}", fontsize=10, pad=8)
    ax.set_xlabel(r"$\bar{m}_{salt}$ / mol kg$^{-1}$", labelpad=4)
    ax.set_ylabel(r"$\gamma_{\pm}^{m,*}$ / -", labelpad=4)


def main() -> None:
    fig, axes = plt.subplots(2, 2, figsize=(8.8, 7.8))
    for ax, cfg in zip(axes.flat, PANELS):
        _plot_panel(ax, *cfg)
    axes.flat[0].legend(loc="upper right", fontsize=8, frameon=False)
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.09, top=0.96, wspace=0.20, hspace=0.32)
    common.save_figure(fig, OUTPUT)
    plt.close(fig)


if __name__ == "__main__":
    main()
