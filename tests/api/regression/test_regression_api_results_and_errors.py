"""Regression API contract tests outside the hydrocarbon benchmark."""

from __future__ import annotations

import csv

import pytest

import epcsaft
from epcsaft import FitProblem, FitResult, create_parameter_template, write_fit_result
from epcsaft._types import InputError


def test_provenance_validation_rejects_indirect_dborn_without_electrostatic_data():
    with pytest.raises(InputError, match=r"d_born.*dielectric|d_born.*ion-activity"):
        epcsaft.validate_regression_provenance(
            [
                epcsaft.FitParameter(
                    "MEAH+",
                    "d_born",
                    source="mixed_reactive_vle",
                )
            ],
            strict=True,
        )

    report = epcsaft.validate_regression_provenance(
        [
            epcsaft.FitParameter(
                "MEAH+",
                "d_born",
                source="explicit_override",
            )
        ],
        strict=True,
    )
    assert report["warnings"] == []
    assert report["parameter_sources"]["MEAH+.d_born"] == "explicit_override"

def test_provenance_validation_gates_ion_involving_binary_interactions():
    with pytest.raises(InputError, match="same-sign ionic pair"):
        epcsaft.validate_regression_provenance(
            [
                epcsaft.BinaryInteraction(
                    ("MEACOO-", "HCO3-"),
                    parameter="k_ij",
                    source="direct_binary_vle",
                )
            ],
            species=["MEACOO-", "HCO3-"],
            charges=[-1.0, -1.0],
            strict=True,
        )

    with pytest.raises(InputError, match=r"neutral-ion.*direct"):
        epcsaft.validate_regression_provenance(
            [
                epcsaft.BinaryInteraction(
                    ("CO2", "MEACOO-"),
                    parameter="k_ij",
                    source="mixed_reactive_vle",
                )
            ],
            species=["CO2", "MEACOO-"],
            charges=[0.0, -1.0],
            strict=True,
        )

    report = epcsaft.validate_regression_provenance(
        [
            epcsaft.BinaryInteraction(
                ("MEAH+", "MEACOO-"),
                parameter="k_ij",
                source="direct_electrolyte_activity",
            )
        ],
        species=["MEAH+", "MEACOO-"],
        charges=[1.0, -1.0],
        strict=True,
    )
    assert report["warnings"] == []
    assert report["data_sources_by_parameter"]["MEAH+:MEACOO-.k_ij"] == ["direct_electrolyte_activity"]

def test_write_fit_result_updates_only_target_pure_row(tmp_path):
    dataset_root = create_parameter_template(tmp_path, "pure_case", ["H2O", "Na+"])
    result = FitResult(
        problem=FitProblem(
            mode="pure_neutral",
            component="H2O",
            fit_targets=("m", "s", "e"),
            optimization_parameters=("m", "s", "e"),
        ),
        fitted_values={"m": 1.2047, "s": 2.7927, "e": 353.95},
        rendered_values={"m": 1.2047, "s": 2.7927, "e": 353.95},
        success=True,
    )

    written = write_fit_result(result, dataset_root, overwrite=False)
    assert written == [dataset_root / "pure" / "water.csv"]

    with (dataset_root / "pure" / "water.csv").open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))

    target_row = next(row for row in rows if row["component"] == "H2O")
    other_row = next(row for row in rows if row["component"] == "Na+")
    assert target_row["m"] == "1.2047"
    assert target_row["s"] == "2.7927"
    assert target_row["e"] == "353.95"
    assert other_row["m"] == ""
    assert other_row["s"] == ""
    assert other_row["e"] == ""

def test_write_fit_result_updates_ion_row_and_binary_matrix_symmetrically(tmp_path):
    ion_root = create_parameter_template(tmp_path, "ion_case", ["H2O", "Na+"])
    ion_result = FitResult(
        problem=FitProblem(
            mode="pure_ion",
            component="Na+",
            fit_targets=("s", "e"),
            optimization_parameters=("s", "e"),
        ),
        fitted_values={"s": 2.84, "e": 231.2},
        rendered_values={"s": 2.84, "e": 231.2},
        success=True,
    )

    written = write_fit_result(ion_result, ion_root, overwrite=False)
    assert written == [ion_root / "pure" / "water.csv"]
    with (ion_root / "pure" / "water.csv").open("r", encoding="utf-8-sig", newline="") as handle:
        ion_rows = list(csv.DictReader(handle))
    ion_row = next(row for row in ion_rows if row["component"] == "Na+")
    assert ion_row["s"] == "2.84"
    assert ion_row["e"] == "231.2"

    binary_root = create_parameter_template(tmp_path, "binary_case", ["H2O", "Ethanol"])
    binary_result = FitResult(
        problem=FitProblem(
            mode="binary_pair",
            pair=("H2O", "Ethanol"),
            fit_targets=("k_ij", "l_ij", "k_hb_ij"),
            optimization_parameters=("k_ij", "l_ij", "k_hb_ij"),
        ),
        fitted_values={"k_ij": -0.06167, "l_ij": 0.0123, "k_hb_ij": -0.0345},
        rendered_values={"k_ij": -0.06167, "l_ij": 0.0123, "k_hb_ij": -0.0345},
        success=True,
    )

    written = write_fit_result(binary_result, binary_root, overwrite=False)
    assert written == [
        binary_root / "mixed" / "binary_interaction" / "k_ij.csv",
        binary_root / "mixed" / "binary_interaction" / "l_ij.csv",
        binary_root / "mixed" / "binary_interaction" / "k_hb_ij.csv",
    ]
    for filename, expected in (("k_ij.csv", "-0.06167"), ("l_ij.csv", "0.0123"), ("k_hb_ij.csv", "-0.0345")):
        with (binary_root / "mixed" / "binary_interaction" / filename).open(
            "r", encoding="utf-8-sig", newline=""
        ) as handle:
            rows = list(csv.reader(handle))
        header = rows[0]
        h2o_col = header.index("H2O")
        ethanol_col = header.index("Ethanol")
        h2o_row = next(row for row in rows if row[0] == "H2O")
        ethanol_row = next(row for row in rows if row[0] == "Ethanol")
        assert h2o_row[ethanol_col] == expected
        assert ethanol_row[h2o_col] == expected
    with pytest.raises(InputError, match="overwrite"):
        write_fit_result(binary_result, binary_root, overwrite=False)
