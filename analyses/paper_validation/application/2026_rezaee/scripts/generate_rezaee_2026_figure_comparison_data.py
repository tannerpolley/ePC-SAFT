from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from _paths import ANALYSIS_DIR

PROCESSED_DIR = ANALYSIS_DIR / "data" / "processed"
RESULTS_DIR = ANALYSIS_DIR / "results" / "reaction_equilibrium"
FIGURES_DIR = ANALYSIS_DIR / "figures"

DIGITIZED_POINTS_CSV = PROCESSED_DIR / "rezaee_2026_paper_figure_digitized_points.csv"
SECTION32_ROWS_CSV = PROCESSED_DIR / "rezaee_2026_section32_equilibrium_replication_rows.csv"
SUMMARY_JSON = RESULTS_DIR / "rezaee_2026_package_figure_comparison_data_summary.json"

FIGURE_SPECS: dict[str, dict[str, Any]] = {
    "fig7": {
        "folder": "figure_7",
        "paper_label": "Fig. 7",
        "caption": "Deviation of calculated extraction percentage from experimental data [17].",
        "axis_max": 60.0,
        "package_case_id": "held_2014_s2_no_born_no_kij_pH_stoich",
        "package_x_column": "Li_extraction_pct_exp",
        "package_y_column": "Li_extraction_pct_calc",
        "quantity": "lithium_extraction_pct",
        "x_label": "Experimental Li extraction / %",
        "y_label": "Calculated Li extraction / %",
    },
    "fig8": {
        "folder": "figure_8",
        "paper_label": "Fig. 8",
        "caption": "Deviation of calculated selectivity from experimental data [17].",
        "axis_max": 6.0,
        "package_case_id": "held_2014_s2_no_born_no_kij_pH_stoich",
        "package_x_column": "selectivity_exp",
        "package_y_column": "selectivity_calc",
        "quantity": "li_na_selectivity",
        "x_label": "Experimental Li/Na selectivity / -",
        "y_label": "Calculated Li/Na selectivity / -",
    },
    "fig10": {
        "folder": "figure_10",
        "paper_label": "Fig. 10",
        "caption": "Deviation of calculated lithium extraction percentage from experimental data [17] using k_ij.",
        "axis_max": 60.0,
        "package_case_id": "held_2014_s2_no_born_table9_kij_pH_stoich",
        "package_x_column": "Li_extraction_pct_exp",
        "package_y_column": "Li_extraction_pct_calc",
        "quantity": "lithium_extraction_pct_table9_kij",
        "x_label": "Experimental Li extraction / %",
        "y_label": "Calculated Li extraction / %",
    },
    "fig11": {
        "folder": "figure_11",
        "paper_label": "Fig. 11",
        "caption": "Deviation of calculated selectivity from experimental data [17] using k_ij.",
        "axis_max": 6.0,
        "package_case_id": "held_2014_s2_no_born_table9_kij_pH_stoich",
        "package_x_column": "selectivity_exp",
        "package_y_column": "selectivity_calc",
        "quantity": "li_na_selectivity_table9_kij",
        "x_label": "Experimental Li/Na selectivity / -",
        "y_label": "Calculated Li/Na selectivity / -",
    },
}


def _jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_jsonable(item) for item in value]
    return value


def _aard_pct(x: pd.Series, y: pd.Series) -> float:
    x_values = x.to_numpy(dtype=float)
    y_values = y.to_numpy(dtype=float)
    return float((100.0 * np.abs(y_values - x_values) / np.clip(x_values, 1.0e-300, None)).mean())


def _paper_rows(digitized: pd.DataFrame, figure_id: str, spec: dict[str, Any]) -> pd.DataFrame:
    rows = digitized.loc[digitized["figure_id"] == figure_id].copy()
    if rows.empty:
        raise ValueError(f"No digitized paper rows found for {figure_id} in {DIGITIZED_POINTS_CSV}")
    rows = rows.rename(columns={"x": "experimental_value", "y": "calculated_value"})
    rows["series"] = "digitized_paper_model"
    rows["model_case_id"] = "digitized_from_published_figure"
    rows["point_id"] = [f"{figure_id}_paper_{index + 1:03d}" for index in range(len(rows))]
    rows["quantity"] = str(spec["quantity"])
    rows["source_role"] = "paper_figure_digitization"
    return rows[
        [
            "figure_id",
            "paper_label",
            "caption",
            "axis_max",
            "quantity",
            "series",
            "model_case_id",
            "point_id",
            "experimental_value",
            "calculated_value",
            "source_role",
        ]
    ]


def _package_rows(package_rows: pd.DataFrame, figure_id: str, spec: dict[str, Any]) -> pd.DataFrame:
    case_id = str(spec["package_case_id"])
    rows = package_rows.loc[package_rows["case_id"] == case_id].copy()
    if rows.empty:
        raise ValueError(f"No current package model rows found for case_id={case_id!r} in {SECTION32_ROWS_CSV}")
    return pd.DataFrame(
        {
            "figure_id": figure_id,
            "paper_label": str(spec["paper_label"]),
            "caption": str(spec["caption"]),
            "axis_max": float(spec["axis_max"]),
            "quantity": str(spec["quantity"]),
            "series": "current_epcsaft_package_model",
            "model_case_id": case_id,
            "point_id": [f"{figure_id}_package_{int(row.experiment_no):03d}" for row in rows.itertuples(index=False)],
            "experimental_value": rows[str(spec["package_x_column"])].to_numpy(dtype=float),
            "calculated_value": rows[str(spec["package_y_column"])].to_numpy(dtype=float),
            "source_role": "current_worktree_epcsaft_section32_replication",
        }
    )


def _write_metadata(folder: Path, figure_id: str, spec: dict[str, Any], data_path: Path) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    text = "\n".join(
        [
            f"figure_id: {figure_id}",
            f"paper_label: {spec['paper_label']}",
            f"description: {spec['caption']}",
            "series:",
            "  - digitized_paper_model",
            "  - current_epcsaft_package_model",
            f"package_case_id: {spec['package_case_id']}",
            f"quantity: {spec['quantity']}",
            f"x_label: {spec['x_label']}",
            f"y_label: {spec['y_label']}",
            f"axis_max: {spec['axis_max']}",
            f"data: {data_path.relative_to(ANALYSIS_DIR).as_posix()}",
            "source_images:",
            f"  source: figures/{spec['folder']}/input/source_image.jpg",
            f"  digitization_overlay: figures/{spec['folder']}/input/digitization_overlay.png",
        ]
    )
    (folder / f"{spec['folder']}.mpl.yaml").write_text(text + "\n", encoding="utf-8")


def main() -> int:
    digitized = pd.read_csv(DIGITIZED_POINTS_CSV)
    package_rows = pd.read_csv(SECTION32_ROWS_CSV)
    summary_entries: list[dict[str, Any]] = []

    for figure_id, spec in FIGURE_SPECS.items():
        figure_folder = FIGURES_DIR / str(spec["folder"])
        data_folder = figure_folder / "output" / "data"
        data_folder.mkdir(parents=True, exist_ok=True)
        combined = pd.concat(
            [_paper_rows(digitized, figure_id, spec), _package_rows(package_rows, figure_id, spec)],
            ignore_index=True,
        )
        out_csv = data_folder / f"{spec['folder']}_package_vs_paper_points.csv"
        combined.to_csv(out_csv, index=False)
        _write_metadata(figure_folder, figure_id, spec, out_csv)

        paper = combined.loc[combined["series"] == "digitized_paper_model"]
        package = combined.loc[combined["series"] == "current_epcsaft_package_model"]
        summary_entries.append(
            {
                "figure_id": figure_id,
                "paper_label": str(spec["paper_label"]),
                "data": str(out_csv.relative_to(ANALYSIS_DIR)),
                "package_case_id": str(spec["package_case_id"]),
                "digitized_paper_points": int(len(paper)),
                "package_model_points": int(len(package)),
                "digitized_paper_aard_pct": _aard_pct(paper["experimental_value"], paper["calculated_value"]),
                "current_epcsaft_package_aard_pct": _aard_pct(
                    package["experimental_value"],
                    package["calculated_value"],
                ),
            }
        )

    summary = {
        "status": "package_figure_comparison_data_generated",
        "source_digitized_points": str(DIGITIZED_POINTS_CSV.relative_to(ANALYSIS_DIR)),
        "source_package_rows": str(SECTION32_ROWS_CSV.relative_to(ANALYSIS_DIR)),
        "figures": summary_entries,
    }
    SUMMARY_JSON.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY_JSON.write_text(json.dumps(_jsonable(summary), indent=2) + "\n", encoding="utf-8")
    print(json.dumps(_jsonable(summary), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
