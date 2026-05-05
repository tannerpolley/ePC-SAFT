from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
ETHANOL_WATER_JCED2021_100KPA = REPO_ROOT / "data" / "regression" / "binary" / "ethanol_water_jced2021_vle_100kpa.csv"
ETHANOL_WATER_PAPER_PCSAFT_KIJ_100KPA = -0.0269
ETHANOL_WATER_HELD_2012_KIJ = -0.049


def ethanol_water_jced2021_vle_records(*, smoke_only: bool = False) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    with ETHANOL_WATER_JCED2021_100KPA.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            if smoke_only and str(row["use_in_smoke"]).strip().lower() != "true":
                continue
            records.append(
                {
                    "T": float(row["T"]),
                    "P": float(row["P"]),
                    "x_Ethanol": float(row["x_Ethanol"]),
                    "x_H2O": float(row["x_H2O"]),
                    "y_Ethanol": float(row["y_Ethanol"]),
                    "y_H2O": float(row["y_H2O"]),
                    "source_id": row["source_id"],
                    "source_doi": row["source_doi"],
                }
            )
    return records
