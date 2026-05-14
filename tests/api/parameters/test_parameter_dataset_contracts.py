from __future__ import annotations

import csv
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
PARAMETER_ROOT = REPO_ROOT / "data" / "reference" / "epcsaft_parameters"


def _assoc_parameter_is_enabled(value: str | None) -> bool:
    if value is None or not str(value).strip():
        return False
    try:
        return float(str(value).strip()) != 0.0
    except ValueError:
        return True


def test_parameter_csvs_use_2b_for_assoc_scheme_tokens() -> None:
    offenders: list[str] = []
    for path in PARAMETER_ROOT.rglob("*.csv"):
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames or "assoc_scheme" not in reader.fieldnames:
                continue
            for line_number, row in enumerate(reader, start=2):
                component = str(row.get("component", "")).strip()
                token = str(row.get("assoc_scheme") or "").strip()
                assoc_enabled = _assoc_parameter_is_enabled(row.get("e_assoc")) or _assoc_parameter_is_enabled(
                    row.get("vol_a")
                )
                if assoc_enabled:
                    if token != "2B":
                        offenders.append(
                            f"{path.relative_to(REPO_ROOT)}:{line_number}:{component}: assoc_scheme={token!r}"
                        )
                elif token:
                    offenders.append(
                        f"{path.relative_to(REPO_ROOT)}:{line_number}:{component}: non-associating token={token!r}"
                    )

    assert offenders == []
