# -*- coding: utf-8 -*-
"""Hydrocarbon neutral-regression benchmark against Gross/Sadowski (2001) Table 2."""

from __future__ import annotations

import pytest

from epcsaft import fit_pure_neutral
from tests.helpers.regression_cases import HYDROCARBON_REFERENCE
from tests.helpers.regression_cases import _load_workbook_reference_rows
from tests.helpers.regression_cases import _neutral_fixed_parameters
from tests.helpers.regression_cases import _real_saturation_records


def _print_benchmark_table(rows: list[dict[str, float]]) -> None:
    print("\nHydrocarbon neutral-component regression benchmark")
    print("Source: Gross/Sadowski (2001) Table 2 targets with NIST saturation pressure and liquid-density data")
    print(
        f"{'Species':<10} {'m_fit':>12} {'m_ref':>12} {'dm':>12} "
        f"{'s_fit':>12} {'s_ref':>12} {'ds':>12} "
        f"{'e_fit':>12} {'e_ref':>12} {'de':>12} "
        f"{'rho_rms':>12} {'vle_rms':>12}"
    )
    for row in rows:
        print(
            f"{row['species']:<10} "
            f"{row['m_fit']:>12.6f} {row['m_ref']:>12.6f} {row['m_delta']:>12.3e} "
            f"{row['s_fit']:>12.6f} {row['s_ref']:>12.6f} {row['s_delta']:>12.3e} "
            f"{row['e_fit']:>12.6f} {row['e_ref']:>12.6f} {row['e_delta']:>12.3e} "
            f"{row['rho_rms']:>12.3e} {row['vle_rms']:>12.3e}"
        )


def test_hydrocarbon_reference_csv_matches_gross_2001_table2():
    csv_rows = _load_workbook_reference_rows()
    assert set(csv_rows) == set(HYDROCARBON_REFERENCE)
    for component, expected in HYDROCARBON_REFERENCE.items():
        for field, expected_value in expected.items():
            assert csv_rows[component][field] == pytest.approx(expected_value, rel=0.0, abs=1.0e-12)


def test_hydrocarbon_neutral_regression_matches_literature_from_real_saturation_data():
    csv_rows = _load_workbook_reference_rows()
    benchmark_rows: list[dict[str, float]] = []

    for component in ("Methane", "Ethane", "Propane"):
        records = _real_saturation_records(component)
        reference = csv_rows[component]
        result = fit_pure_neutral(
            records,
            component,
            assoc_scheme="",
            fixed_parameters=_neutral_fixed_parameters(component),
            initial_guess={
                "m": reference["m"] * 1.08,
                "s": reference["s"] * 0.96,
                "e": reference["e"] * 1.05,
            },
            bounds={
                "m": (0.5, 3.5),
                "s": (2.0, 5.0),
                "e": (50.0, 400.0),
            },
        )

        assert result.success, result.message
        assert result.backend == "least_squares_native"
        assert result.problem.mode == "pure_neutral"
        assert result.metrics_by_term["density"] < 0.01
        assert result.metrics_by_term["pure_vle_fugacity_balance"] < 0.01
        assert result.fitted_values["m"] == pytest.approx(reference["m"], rel=0.0, abs=0.03)
        assert result.fitted_values["s"] == pytest.approx(reference["s"], rel=0.0, abs=0.05)
        assert result.fitted_values["e"] == pytest.approx(reference["e"], rel=0.0, abs=2.0)

        benchmark_rows.append(
            {
                "species": component,
                "m_fit": result.fitted_values["m"],
                "m_ref": reference["m"],
                "m_delta": result.fitted_values["m"] - reference["m"],
                "s_fit": result.fitted_values["s"],
                "s_ref": reference["s"],
                "s_delta": result.fitted_values["s"] - reference["s"],
                "e_fit": result.fitted_values["e"],
                "e_ref": reference["e"],
                "e_delta": result.fitted_values["e"] - reference["e"],
                "rho_rms": result.metrics_by_term["density"],
                "vle_rms": result.metrics_by_term["pure_vle_fugacity_balance"],
            }
        )

    _print_benchmark_table(benchmark_rows)
