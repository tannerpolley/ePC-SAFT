# -*- coding: utf-8 -*-
"""Fast hydrocarbon neutral-regression contracts against Gross/Sadowski (2001)."""

from __future__ import annotations

import pytest

import numpy as np

from epcsaft.regression import _debug_native_pure_neutral_objective
from tests.helpers.regression_cases import HYDROCARBON_REFERENCE
from tests.helpers.regression_cases import _load_workbook_reference_rows
from tests.helpers.regression_cases import _neutral_fixed_parameters
from tests.helpers.regression_cases import _real_saturation_records


def test_hydrocarbon_reference_csv_matches_gross_2001_table2():
    csv_rows = _load_workbook_reference_rows()
    assert set(csv_rows) == set(HYDROCARBON_REFERENCE)
    for component, expected in HYDROCARBON_REFERENCE.items():
        for field, expected_value in expected.items():
            assert csv_rows[component][field] == pytest.approx(expected_value, rel=0.0, abs=1.0e-12)


def test_methane_reference_parameters_keep_native_objective_pinned():
    csv_rows = _load_workbook_reference_rows()
    reference = csv_rows["Methane"]

    debug = _debug_native_pure_neutral_objective(
        _real_saturation_records("Methane"),
        "Methane",
        assoc_scheme="",
        fixed_parameters=_neutral_fixed_parameters("Methane"),
        initial_guess=reference,
        x=reference,
    )

    assert debug["objective"] == pytest.approx(9.701615164740784e-06, rel=2.0e-4)
    assert debug["jacobian_shape"] == (8, 3)
    assert np.asarray(debug["density_raw_residuals"], dtype=float).shape == (4,)
    assert np.asarray(debug["pure_vle_raw_residuals"], dtype=float).shape == (4,)
    assert debug["density_solves"] >= 4
    assert debug["fused_state_evaluations"] >= 4
