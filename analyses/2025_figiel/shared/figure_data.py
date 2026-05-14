from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
ANALYSIS_ROOT = Path(__file__).resolve().parents[1]
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))

from shared import _common as common

SOURCE_ROOT = Path(__file__).resolve().parent
FIGURES_ROOT = ANALYSIS_ROOT / "figures"
DATASET = "2025_Figiel"

FIELDNAMES = [
    "figure_id",
    "panel_id",
    "series_id",
    "series_label",
    "source_type",
    "dataset",
    "parameter_set",
    "salt",
    "ion",
    "solvent_system",
    "composition",
    "composition_basis",
    "target_w_org",
    "x_name",
    "x_value",
    "x_units",
    "y_name",
    "y_value",
    "y_units",
    "curve_points",
    "user_options",
    "source_file",
    "category",
    "display_order",
]
NUMERIC_FIELDS = {"target_w_org", "x_value", "y_value", "curve_points", "display_order"}

FIG5_PANELS = [
    ("a)", ["LiCl", "NaCl", "KCl"], "Cl$^-$ salts in water"),
    ("b)", ["LiBr", "NaBr", "KBr"], "Br$^-$ salts in water"),
]
FIG8_PANELS = [
    ("a)", "LiBr", 5.0, 3.0),
    ("b)", "NaI", 1.5, 1.125),
    ("c)", "NaBr", 1.5, 1.125),
]
FIG8_CURVE_MAX_OVERRIDES = {
    ("NaI", "methanol"): 0.77,
    ("NaBr", "ethanol"): 0.195,
}
FIG9_PANELS = [
    ("a)", "NaBr", "water-methanol", 3.0),
    ("b)", "NaBr", "water-ethanol", 4.0),
    ("c)", "NaCl", "water-methanol", 2.0),
    ("d)", "NaCl", "water-ethanol", 2.0),
]
FIG9_TARGETS = (0.8, 0.4)
FIG9_DATA_TOL = 0.08
FIG9_CURVE_MAX_BY_PANEL = {
    "b)": {0.8: 1.5},
    "c)": {0.8: 0.6},
    "d)": {0.8: 0.3},
}


def payload_path(figure_id: str) -> Path:
    return FIGURES_ROOT / figure_id / "output" / f"{figure_id}_series.csv"


def figure_input_dir(figure_id: str) -> Path:
    return FIGURES_ROOT / figure_id / "input"


def figure_input_path(figure_id: str, *parts: str) -> Path:
    return figure_input_dir(figure_id).joinpath(*parts)


def _format_value(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if math.isnan(value):
            return "nan"
        return format(value, ".17g")
    if isinstance(value, int):
        return str(value)
    return str(value)


def _composition_json(comp: dict[str, float] | None) -> str:
    if not comp:
        return ""
    return json.dumps({k: float(v) for k, v in sorted(comp.items())}, sort_keys=True, separators=(",", ":"))


def _options_json(options: dict | None) -> str:
    return json.dumps(options or {}, sort_keys=True, separators=(",", ":"))


def _source_ref(path: Path | None) -> str:
    return path.relative_to(REPO_ROOT).as_posix() if path is not None else ""


def _row(
    *,
    figure_id: str,
    panel_id: str,
    series_id: str,
    series_label: str,
    source_type: str,
    x_name: str,
    x_value: float,
    y_name: str,
    y_value: float,
    display_order: int,
    dataset: str = "",
    parameter_set: str = "",
    salt: str = "",
    ion: str = "",
    solvent_system: str = "",
    composition: dict[str, float] | None = None,
    composition_basis: str = "",
    target_w_org: float | None = None,
    x_units: str = "",
    y_units: str = "",
    curve_points: int | None = None,
    user_options: dict | None = None,
    source_file: Path | None = None,
    category: str = "",
) -> dict[str, str]:
    data = {
        "figure_id": figure_id,
        "panel_id": panel_id,
        "series_id": series_id,
        "series_label": series_label,
        "source_type": source_type,
        "dataset": dataset,
        "parameter_set": parameter_set,
        "salt": salt,
        "ion": ion,
        "solvent_system": solvent_system,
        "composition": _composition_json(composition),
        "composition_basis": composition_basis,
        "target_w_org": target_w_org,
        "x_name": x_name,
        "x_value": x_value,
        "x_units": x_units,
        "y_name": y_name,
        "y_value": y_value,
        "y_units": y_units,
        "curve_points": curve_points,
        "user_options": _options_json(user_options),
        "source_file": _source_ref(source_file),
        "category": category,
        "display_order": display_order,
    }
    return {field: _format_value(data.get(field, "")) for field in FIELDNAMES}


def _write_payload(figure_id: str, rows: list[dict[str, str]]) -> Path:
    path = payload_path(figure_id)
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=FIELDNAMES).to_csv(path, index=False)
    return path


def read_payload(figure_id: str) -> list[dict[str, str]]:
    path = payload_path(figure_id)
    if not path.exists():
        raise FileNotFoundError(f"Missing generated payload {path}. Run generate_figure_data.py first.")
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def select_rows(rows: Iterable[dict[str, str]], **criteria: str) -> list[dict[str, str]]:
    out = []
    for row in rows:
        if all(row.get(key, "") == str(value) for key, value in criteria.items()):
            out.append(row)
    return sorted(out, key=lambda row: (float(row.get("display_order") or 0.0), float(row.get("x_value") or 0.0)))


def xy(rows: Iterable[dict[str, str]]) -> tuple[list[float], list[float]]:
    selected = list(rows)
    return [float(row["x_value"]) for row in selected], [float(row["y_value"]) for row in selected]


def rows_by_key(rows: Iterable[dict[str, str]], key: str) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(row.get(key, ""), []).append(row)
    return {name: select_rows(values) for name, values in grouped.items()}


def generate_figure_4() -> list[dict[str, str]]:
    figure_id = "figure_4"
    source = figure_input_path(figure_id, "water.csv")
    _, raw_rows = common.read_csv_rows(source)
    literature: dict[str, float] = {}
    for raw in raw_rows:
        ion = str(raw.get("Ion", "")).strip()
        value = common.parse_float(raw.get("Gsolv (kJ/mol)"))
        if ion and value is not None:
            key = f"{ion}+" if ion in {"H", "Li", "Na", "K"} else f"{ion}-" if ion in {"Cl", "Br", "I"} else ion
            literature[key] = value

    rows: list[dict[str, str]] = []
    ions = ["Li+", "Na+", "K+", "Cl-", "Br-", "I-"]
    for order, ion in enumerate(ions):
        rows.append(
            _row(
                figure_id=figure_id,
                panel_id="combined",
                series_id="literature",
                series_label="Literature",
                source_type="literature",
                ion=ion,
                solvent_system="water",
                x_name="ion_index",
                x_value=float(order),
                y_name="negative_gsolv",
                y_value=-float(literature.get(ion, float("nan"))),
                y_units="kJ/mol",
                source_file=source,
                category=ion,
                display_order=order,
            )
        )
        for dataset, label in (("2025_Figiel", "ePC-SAFT 2025"), ("2020_Bulow", "ePC-SAFT 2020")):
            try:
                y_value = -common.gsolv_ion(dataset, ion, "water", {"water": 1.0})
            except Exception:
                y_value = float("nan")
            rows.append(
                _row(
                    figure_id=figure_id,
                    panel_id="combined",
                    series_id=f"model_{dataset}",
                    series_label=label,
                    source_type="model",
                    dataset=dataset,
                    parameter_set=dataset,
                    ion=ion,
                    solvent_system="water",
                    composition={"water": 1.0},
                    composition_basis="mole_fraction",
                    x_name="ion_index",
                    x_value=float(order),
                    y_name="negative_gsolv",
                    y_value=float(y_value),
                    y_units="kJ/mol",
                    curve_points=1,
                    category=ion,
                    display_order=order,
                )
            )
    return rows


def generate_figure_5() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    source_root = figure_input_dir("figure_5") / "water"
    order = 0
    for panel_id, salts, _title in FIG5_PANELS:
        for salt in salts:
            source = source_root / f"water-{salt}.csv"
            for entry in common.read_miac_dataset(source, "water"):
                rows.append(
                    _row(
                        figure_id="figure_5",
                        panel_id=panel_id,
                        series_id=f"data_{salt}",
                        series_label=f"{salt} data",
                        source_type="literature",
                        salt=salt,
                        solvent_system="water",
                        composition={"water": 1.0},
                        composition_basis="mole_fraction",
                        x_name="molality",
                        x_value=float(entry["molality"]),
                        x_units="mol/kg",
                        y_name="miac_m",
                        y_value=float(entry["miac_m"]),
                        y_units="-",
                        source_file=source,
                        category=salt,
                        display_order=order,
                    )
                )
                order += 1
            m_grid, y_model = common.mean_ionic_activity_curve(DATASET, salt, "water", {"water": 1.0}, 6.0, points=600)
            for m, y in zip(m_grid, y_model):
                rows.append(
                    _row(
                        figure_id="figure_5",
                        panel_id=panel_id,
                        series_id=f"model_{salt}",
                        series_label=f"{salt} model",
                        source_type="model",
                        dataset=DATASET,
                        parameter_set=DATASET,
                        salt=salt,
                        solvent_system="water",
                        composition={"water": 1.0},
                        composition_basis="mole_fraction",
                        x_name="molality",
                        x_value=float(m),
                        x_units="mol/kg",
                        y_name="miac_m",
                        y_value=float(y),
                        y_units="-",
                        curve_points=600,
                        category=salt,
                        display_order=order,
                    )
                )
                order += 1
    return rows


def _load_transfer_xy(path: Path) -> tuple[np.ndarray, np.ndarray]:
    fields, raw_rows = common.read_csv_rows(path)
    x_key, y_key = fields[0], fields[1]
    xs, ys = [], []
    for row in raw_rows:
        x = common.parse_float(row.get(x_key))
        y = common.parse_float(row.get(y_key))
        if x is not None and y is not None:
            xs.append(x)
            ys.append(y)
    return np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)


def generate_figure_6() -> list[dict[str, str]]:
    data_root = figure_input_dir("figure_6") / "G_trans" / "water"
    panels = [
        ("a)", "K+", "methanol", data_root / "methanol" / "K.csv"),
        ("b)", "Br-", "methanol", data_root / "methanol" / "Br.csv"),
        ("c)", "Na+", "ethanol", data_root / "ethanol" / "Na.csv"),
        ("d)", "Cl-", "ethanol", data_root / "ethanol" / "Cl.csv"),
    ]
    rows: list[dict[str, str]] = []
    order = 0
    for panel_id, ion, organic, source in panels:
        x_data, y_data = _load_transfer_xy(source)
        solvent_system = f"water-{organic}"
        for x_value, y_value in zip(x_data, y_data):
            rows.append(
                _row(
                    figure_id="figure_6",
                    panel_id=panel_id,
                    series_id=f"data_{ion}_{organic}",
                    series_label="Literature data",
                    source_type="literature",
                    ion=ion,
                    solvent_system=solvent_system,
                    x_name=f"x_{organic}",
                    x_value=float(x_value),
                    x_units="-",
                    y_name="gtrans",
                    y_value=float(y_value),
                    y_units="kJ/mol",
                    source_file=source,
                    category=f"{ion}_{organic}",
                    display_order=order,
                )
            )
            order += 1
        x_grid = np.linspace(0.0, 1.0, 401)
        y_model = common.transfer_curve(DATASET, ion, organic, x_grid)
        for x_value, y_value in zip(x_grid, y_model):
            rows.append(
                _row(
                    figure_id="figure_6",
                    panel_id=panel_id,
                    series_id=f"model_{ion}_{organic}",
                    series_label="ePC-SAFT fit",
                    source_type="model",
                    dataset=DATASET,
                    parameter_set=DATASET,
                    ion=ion,
                    solvent_system=solvent_system,
                    composition_basis="mole_fraction",
                    x_name=f"x_{organic}",
                    x_value=float(x_value),
                    x_units="-",
                    y_name="gtrans",
                    y_value=float(y_value),
                    y_units="kJ/mol",
                    curve_points=401,
                    category=f"{ion}_{organic}",
                    display_order=order,
                )
            )
            order += 1
    return rows


def generate_figure_7() -> list[dict[str, str]]:
    source = figure_input_path("figure_7", "methanol", "methanol-NaBr.csv")
    data = common.read_miac_dataset(source, "methanol")
    m_max = max(float(row["molality"]) for row in data)
    rows: list[dict[str, str]] = []
    order = 0
    for entry in data:
        rows.append(
            _row(
                figure_id="figure_7",
                panel_id="main",
                series_id="data_NaBr_methanol",
                series_label="Literature",
                source_type="literature",
                salt="NaBr",
                solvent_system="methanol",
                composition={"methanol": 1.0},
                composition_basis="mole_fraction",
                x_name="molality",
                x_value=float(entry["molality"]),
                x_units="mol/kg",
                y_name="miac_m",
                y_value=float(entry["miac_m"]),
                y_units="-",
                source_file=source,
                category="NaBr",
                display_order=order,
            )
        )
        order += 1
    for series_id, label, options in (
        ("model_default", "Figiel 2025", None),
        ("model_rule1", "Rule 1", {"elec_model": {"rel_perm": {"rule": 1}}}),
    ):
        m_grid, y_model = common.mean_ionic_activity_curve(
            DATASET, "NaBr", "methanol", {"methanol": 1.0}, m_max, points=500, user_options=options
        )
        for m, y in zip(m_grid, y_model):
            rows.append(
                _row(
                    figure_id="figure_7",
                    panel_id="main",
                    series_id=series_id,
                    series_label=label,
                    source_type="model",
                    dataset=DATASET,
                    parameter_set=DATASET,
                    salt="NaBr",
                    solvent_system="methanol",
                    composition={"methanol": 1.0},
                    composition_basis="mole_fraction",
                    x_name="molality",
                    x_value=float(m),
                    x_units="mol/kg",
                    y_name="miac_m",
                    y_value=float(y),
                    y_units="-",
                    curve_points=500,
                    user_options=options,
                    category="NaBr",
                    display_order=order,
                )
            )
            order += 1
    return rows


def generate_figure_8() -> list[dict[str, str]]:
    source_root = figure_input_dir("figure_8")
    rows: list[dict[str, str]] = []
    order = 0
    for panel_id, salt, m_max, _y_max in FIG8_PANELS:
        for solvent, label in (("methanol", "Methanol"), ("ethanol", "Ethanol")):
            source = source_root / solvent / f"{solvent}-{salt}.csv"
            for entry in common.read_miac_dataset(source, solvent):
                rows.append(
                    _row(
                        figure_id="figure_8",
                        panel_id=panel_id,
                        series_id=f"data_{salt}_{solvent}",
                        series_label=f"{label} data",
                        source_type="literature",
                        salt=salt,
                        solvent_system=solvent,
                        composition={solvent: 1.0},
                        composition_basis="mole_fraction",
                        x_name="molality",
                        x_value=float(entry["molality"]),
                        x_units="mol/kg",
                        y_name="miac_m",
                        y_value=float(entry["miac_m"]),
                        y_units="-",
                        source_file=source,
                        category=f"{salt}_{solvent}",
                        display_order=order,
                    )
                )
                order += 1
            curve_max = FIG8_CURVE_MAX_OVERRIDES.get((salt, solvent), m_max)
            m_grid, y_model = common.mean_ionic_activity_curve(
                DATASET, salt, solvent, {solvent: 1.0}, curve_max, points=600
            )
            for m, y in zip(m_grid, y_model):
                rows.append(
                    _row(
                        figure_id="figure_8",
                        panel_id=panel_id,
                        series_id=f"model_{salt}_{solvent}",
                        series_label=f"{label} model",
                        source_type="model",
                        dataset=DATASET,
                        parameter_set=DATASET,
                        salt=salt,
                        solvent_system=solvent,
                        composition={solvent: 1.0},
                        composition_basis="mole_fraction",
                        x_name="molality",
                        x_value=float(m),
                        x_units="mol/kg",
                        y_name="miac_m",
                        y_value=float(y),
                        y_units="-",
                        curve_points=600,
                        category=f"{salt}_{solvent}",
                        display_order=order,
                    )
                )
                order += 1
    return rows


def _read_weight_fraction_dataset(path: Path, solvent_system: str) -> list[dict[str, object]]:
    fields, raw_rows = common.read_csv_rows(path)
    lookup = {field.lower(): field for field in fields}
    molality_key = next((lookup[c] for c in ("molality", "m") if c in lookup), None)
    gamma_key = next((lookup[c] for c in ("miac_m", "gamma") if c in lookup), None)
    if molality_key is None or gamma_key is None:
        raise ValueError(f"Missing columns in {path}")
    organic = [s for s in solvent_system.split("-") if s != "water"][0]
    w_org_key = lookup.get(f"w_{organic}_salt_free".lower()) or lookup.get(f"w_{organic}".lower())
    w_water_key = (
        lookup.get("w_h2o_salt_free") or lookup.get("w_water_salt_free") or lookup.get("w_h2o") or lookup.get("w_water")
    )
    org_key = lookup.get(f"x_{organic}".lower()) or (
        lookup.get("x_methanol") if organic == "methanol" else lookup.get("x_ethanol")
    )
    water_key = lookup.get("x_h2o") or lookup.get("x_water")
    out = []
    for row in raw_rows:
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
        out.append(
            {
                "molality": m,
                "miac_m": y,
                "w_org": weight_comp[organic],
                "weight_signature": tuple((s, round(weight_comp[s], 6)) for s in weight_comp),
                "weight_comp": weight_comp,
            }
        )
    out.sort(key=lambda item: float(item["molality"]))
    return out


def _group_by_weight(
    entries: Iterable[dict[str, object]],
) -> dict[tuple[tuple[str, float], ...], list[dict[str, object]]]:
    grouped: dict[tuple[tuple[str, float], ...], list[dict[str, object]]] = {}
    for entry in entries:
        grouped.setdefault(entry["weight_signature"], []).append(entry)
    for values in grouped.values():
        values.sort(key=lambda item: float(item["molality"]))
    return grouped


def _closest_group(
    entries: list[dict[str, object]], target_w_org: float, tol: float = FIG9_DATA_TOL
) -> list[dict[str, object]] | None:
    candidates = []
    for values in _group_by_weight(entries).values():
        w_org = float(values[0]["w_org"])
        candidates.append((abs(w_org - target_w_org), values))
    if not candidates:
        return None
    candidates.sort(key=lambda item: item[0])
    return candidates[0][1] if candidates[0][0] <= tol else None


def _figure9_curve_m_max(panel_id: str, target_w_org: float, default_max: float) -> float:
    return min(default_max, FIG9_CURVE_MAX_BY_PANEL.get(panel_id, {}).get(target_w_org, default_max))


def generate_figure_9() -> list[dict[str, str]]:
    source_root = figure_input_dir("figure_9")
    rows: list[dict[str, str]] = []
    order = 0
    for panel_id, salt, solvent_system, m_max in FIG9_PANELS:
        source = source_root / solvent_system / f"{solvent_system}-{salt}.csv"
        entries = _read_weight_fraction_dataset(source, solvent_system) if source.exists() else []
        for target_w_org in FIG9_TARGETS:
            selected = _closest_group(entries, target_w_org) if entries else None
            if selected:
                for entry in selected:
                    rows.append(
                        _row(
                            figure_id="figure_9",
                            panel_id=panel_id,
                            series_id=f"data_{salt}_{solvent_system}_w{target_w_org:.1f}",
                            series_label=f"w_org={target_w_org:.1f} data",
                            source_type="literature",
                            salt=salt,
                            solvent_system=solvent_system,
                            composition=entry["weight_comp"],
                            composition_basis="salt_free_weight_fraction",
                            target_w_org=target_w_org,
                            x_name="molality",
                            x_value=float(entry["molality"]),
                            x_units="mol/kg",
                            y_name="miac_m",
                            y_value=float(entry["miac_m"]),
                            y_units="-",
                            source_file=source,
                            category=f"{salt}_{solvent_system}_w{target_w_org:.1f}",
                            display_order=order,
                        )
                    )
                    order += 1
            curve_max = _figure9_curve_m_max(panel_id, target_w_org, m_max)
            comp_model = common.target_weight_fraction_to_comp(solvent_system, target_w_org)
            m_grid, y_model = common.mean_ionic_activity_curve(
                DATASET, salt, solvent_system, comp_model, curve_max, points=600
            )
            for m, y in zip(m_grid, y_model):
                rows.append(
                    _row(
                        figure_id="figure_9",
                        panel_id=panel_id,
                        series_id=f"model_{salt}_{solvent_system}_w{target_w_org:.1f}",
                        series_label=f"w_org={target_w_org:.1f} model",
                        source_type="model",
                        dataset=DATASET,
                        parameter_set=DATASET,
                        salt=salt,
                        solvent_system=solvent_system,
                        composition=comp_model,
                        composition_basis="salt_free_weight_fraction",
                        target_w_org=target_w_org,
                        x_name="molality",
                        x_value=float(m),
                        x_units="mol/kg",
                        y_name="miac_m",
                        y_value=float(y),
                        y_units="-",
                        curve_points=600,
                        category=f"{salt}_{solvent_system}_w{target_w_org:.1f}",
                        display_order=order,
                    )
                )
                order += 1

    for solvent_system, salts, x_max in (
        ("water-methanol", ["NaBr", "NaCl"], 3.3),
        ("water-ethanol", ["NaBr", "NaCl"], 4.4),
    ):
        comp_model = common.target_weight_fraction_to_comp(solvent_system, 0.4)
        for salt in salts:
            source = source_root / solvent_system / f"{solvent_system}-{salt}.csv"
            entries = _read_weight_fraction_dataset(source, solvent_system) if source.exists() else []
            selected = _closest_group(entries, 0.4) if entries else None
            if not selected:
                continue
            for entry in selected:
                rows.append(
                    _row(
                        figure_id="figure_9",
                        panel_id=f"40wt_{solvent_system}",
                        series_id=f"data_{salt}_{solvent_system}_40wt",
                        series_label=f"{salt} data",
                        source_type="literature",
                        salt=salt,
                        solvent_system=solvent_system,
                        composition=entry["weight_comp"],
                        composition_basis="salt_free_weight_fraction",
                        target_w_org=0.4,
                        x_name="molality",
                        x_value=float(entry["molality"]),
                        x_units="mol/kg",
                        y_name="miac_m",
                        y_value=float(entry["miac_m"]),
                        y_units="-",
                        source_file=source,
                        category=f"{salt}_{solvent_system}_40wt",
                        display_order=order,
                    )
                )
                order += 1
            curve_max = min(x_max, max(float(entry["molality"]) for entry in selected) * 1.1)
            m_grid, y_model = common.mean_ionic_activity_curve(
                DATASET, salt, solvent_system, comp_model, curve_max, points=600
            )
            for m, y in zip(m_grid, y_model):
                rows.append(
                    _row(
                        figure_id="figure_9",
                        panel_id=f"40wt_{solvent_system}",
                        series_id=f"model_{salt}_{solvent_system}_40wt",
                        series_label=f"{salt} fit",
                        source_type="model",
                        dataset=DATASET,
                        parameter_set=DATASET,
                        salt=salt,
                        solvent_system=solvent_system,
                        composition=comp_model,
                        composition_basis="salt_free_weight_fraction",
                        target_w_org=0.4,
                        x_name="molality",
                        x_value=float(m),
                        x_units="mol/kg",
                        y_name="miac_m",
                        y_value=float(y),
                        y_units="-",
                        curve_points=600,
                        category=f"{salt}_{solvent_system}_40wt",
                        display_order=order,
                    )
                )
                order += 1
    return rows


GENERATORS = {
    "figure_4": generate_figure_4,
    "figure_5": generate_figure_5,
    "figure_6": generate_figure_6,
    "figure_7": generate_figure_7,
    "figure_8": generate_figure_8,
    "figure_9": generate_figure_9,
}


def write_figure(figure_id: str) -> Path:
    return _write_payload(figure_id, GENERATORS[figure_id]())


def write_all() -> list[Path]:
    return [write_figure(figure_id) for figure_id in GENERATORS]


def _compare_values(expected: str, actual: str, field: str, rtol: float, atol: float) -> bool:
    if field not in NUMERIC_FIELDS:
        return expected == actual
    if expected == actual:
        return True
    if not expected and not actual:
        return True
    try:
        exp = float(expected)
        act = float(actual)
    except ValueError:
        return expected == actual
    if math.isnan(exp) and math.isnan(act):
        return True
    return math.isclose(exp, act, rel_tol=rtol, abs_tol=atol)


def compare_all(rtol: float = 1e-9, atol: float = 1e-10) -> list[str]:
    failures: list[str] = []
    for figure_id, generator in GENERATORS.items():
        baseline_path = payload_path(figure_id)
        if not baseline_path.exists():
            failures.append(f"{figure_id}: missing baseline {baseline_path}")
            continue
        baseline = read_payload(figure_id)
        current = generator()
        if len(baseline) != len(current):
            failures.append(f"{figure_id}: row count changed, baseline={len(baseline)} current={len(current)}")
            continue
        for index, (expected, actual) in enumerate(zip(baseline, current), start=2):
            for field in FIELDNAMES:
                if not _compare_values(expected.get(field, ""), actual.get(field, ""), field, rtol, atol):
                    failures.append(
                        f"{figure_id}:{index}: {field} changed from {expected.get(field, '')!r} to {actual.get(field, '')!r}"
                    )
                    break
    return failures


def main() -> None:
    paths = write_all()
    for path in paths:
        print(path.relative_to(REPO_ROOT))


if __name__ == "__main__":
    main()
