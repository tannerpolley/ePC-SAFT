"""Helpers for creating user-owned ePC-SAFT parameter templates."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Iterable


PURE_TEMPLATE_COLUMNS = [
    "component",
    "m",
    "s",
    "e",
    "e_assoc",
    "vol_a",
    "assoc_scheme",
    "dipm",
    "dip_num",
    "z",
    "dielc",
    "d_born",
    "f_solv",
    "MW",
]

MATRIX_TEMPLATE_COLUMNS = ["component"]
REL_PERM_TEMPLATE_COLUMNS = ["organic", "a", "b", "c"]


def _prompt(prompt: str) -> str:
    value = input(f"{prompt}: ").strip()
    if not value:
        raise ValueError(f"{prompt} is required.")
    return value


def _normalize_species(species: Iterable[str] | str) -> list[str]:
    if isinstance(species, str):
        items = species.split(",")
    else:
        items = species
    normalized = [str(item).strip() for item in items if str(item).strip()]
    if not normalized:
        raise ValueError("At least one species must be provided.")
    if len(set(normalized)) != len(normalized):
        raise ValueError("Species names must be unique.")
    return normalized


def _infer_pure_template_name(species: list[str]) -> str:
    solvent_aliases = {
        "h2o": "water",
        "water": "water",
        "methanol": "methanol",
        "meoh": "methanol",
        "ethanol": "ethanol",
        "etoh": "ethanol",
        "propanol": "propanol",
        "butanol": "butanol",
        "isobutanol": "butanol",
    }
    neutrals = [name.strip().lower() for name in species if not name.strip().endswith("+") and not name.strip().endswith("-")]
    if len(neutrals) == 1 and neutrals[0] in solvent_aliases:
        return f"{solvent_aliases[neutrals[0]]}.csv"
    return "any_solvent.csv"


def _write_csv(path: Path, header: list[str], rows: list[list[str | float | int | None]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for row in rows:
            writer.writerow(["" if value is None else value for value in row])


def _write_matrix_template(path: Path, species: list[str]) -> None:
    header = MATRIX_TEMPLATE_COLUMNS + species
    rows: list[list[str | float | int | None]] = []
    for row_species in species:
        row = [row_species]
        for col_species in species:
            row.append(0.0 if row_species == col_species else None)
        rows.append(row)
    _write_csv(path, header, rows)


def _write_pure_template(path: Path, species: list[str]) -> None:
    rows = [[name] + [None] * (len(PURE_TEMPLATE_COLUMNS) - 1) for name in species]
    _write_csv(path, PURE_TEMPLATE_COLUMNS, rows)


def create_parameter_template(
    location: str | Path | None = None,
    folder_name: str | None = None,
    species: Iterable[str] | str | None = None,
) -> Path:
    """Create a user-owned dataset scaffold and return its root path.

    If any of the inputs are omitted, the function prompts for them.
    The generated layout matches the loader expectations used by
    ``ePCSAFTMixture.from_dataset(...)``.
    """

    if location is None:
        location = _prompt("Template location")
    if folder_name is None:
        folder_name = _prompt("Template folder name")
    if species is None:
        species = _prompt("Comma-separated species list")

    root = Path(location).expanduser() / str(folder_name).strip()
    if root.exists():
        raise FileExistsError(f"Template folder already exists: {root}")

    species_list = _normalize_species(species)

    (root / "pure").mkdir(parents=True, exist_ok=False)
    (root / "mixed" / "binary_interaction").mkdir(parents=True, exist_ok=False)
    (root / "mixed" / "rel_perm").mkdir(parents=True, exist_ok=False)

    pure_name = _infer_pure_template_name(species_list)
    _write_pure_template(root / "pure" / pure_name, species_list)
    _write_matrix_template(root / "mixed" / "binary_interaction" / "k_ij.csv", species_list)
    _write_matrix_template(root / "mixed" / "binary_interaction" / "l_ij.csv", species_list)
    _write_matrix_template(root / "mixed" / "binary_interaction" / "k_hb_ij.csv", species_list)
    _write_csv(root / "mixed" / "rel_perm" / "parameters.csv", REL_PERM_TEMPLATE_COLUMNS, [])

    with (root / "user_options.json").open("w", encoding="utf-8") as handle:
        json.dump({}, handle, indent=2)
        handle.write("\n")

    return root
