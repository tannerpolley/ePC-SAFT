from __future__ import annotations

import numpy as np

import epcsaft
from epcsaft.regression import _fit_pure_neutral_least_squares_internal
from tests.helpers.regression_cases import _load_workbook_reference_rows
from tests.helpers.regression_cases import _methane_like_records
from tests.helpers.regression_cases import _minimal_neutral_metadata
from tests.helpers.regression_cases import _neutral_fixed_parameters
from tests.helpers.regression_cases import _real_saturation_records
from tests.plots.plot_helpers import finite_difference_gradient_values
from tests.plots.plot_helpers import save_comparison_plot
from tests.plots.plot_helpers import save_parity_plot


def test_regression_parameter_reference_comparison_plot() -> None:
    result = _fit_pure_neutral_least_squares_internal(
        _methane_like_records(),
        "Methane",
        assoc_scheme="",
        fixed_parameters=_minimal_neutral_metadata(16.043e-3),
        initial_guess={"m": 1.08, "s": 3.55, "e": 155.0},
        bounds={"m": (0.5, 3.5), "s": (2.0, 5.0), "e": (50.0, 400.0)},
    )
    expected_values = {"m": 1.0, "sigma": 3.7039, "epsilon": 150.03}
    actual_values = {
        "m": result.fitted_values["m"],
        "sigma": result.fitted_values["s"],
        "epsilon": result.fitted_values["e"],
    }

    labels = list(expected_values)
    save_comparison_plot(
        "regression_methane_parameter_reference_comparison.png",
        "Methane fitted parameters vs literature reference values",
        labels,
        np.asarray([actual_values[label] for label in labels], dtype=float),
        np.asarray([expected_values[label] for label in labels], dtype=float),
        category=("regression", "hydrocarbon"),
    )


def test_full_hydrocarbon_basis_parameter_comparison_plot() -> None:
    csv_rows = _load_workbook_reference_rows()
    labels: list[str] = []
    actual: list[float] = []
    expected: list[float] = []

    for component in ("Methane", "Ethane", "Propane"):
        reference = csv_rows[component]
        result = epcsaft.fit_pure_neutral(
            _real_saturation_records(component),
            component,
            assoc_scheme="",
            fixed_parameters=_neutral_fixed_parameters(component),
            initial_guess={
                "m": reference["m"] * 1.08,
                "s": reference["s"] * 0.96,
                "e": reference["e"] * 1.05,
            },
            bounds={"m": (0.5, 3.5), "s": (2.0, 5.0), "e": (50.0, 400.0)},
        )
        for rendered_name, key in (("m", "m"), ("sigma", "s"), ("epsilon", "e")):
            labels.append(f"{component} {rendered_name}")
            actual.append(float(result.fitted_values[key]))
            expected.append(float(reference[key]))

    save_comparison_plot(
        "regression_hydrocarbon_basis_parameter_comparison.png",
        "Hydrocarbon fitted parameters vs workbook/literature references",
        labels,
        np.asarray(actual, dtype=float),
        np.asarray(expected, dtype=float),
        category=("regression", "hydrocarbon"),
    )


def test_regression_gradient_finite_difference_parity_plot() -> None:
    actual, expected = finite_difference_gradient_values()
    save_parity_plot(
        "regression_gradient_finite_difference_parity.png",
        "Native regression gradient vs finite difference",
        ["m", "sigma", "epsilon"],
        actual,
        expected,
        category=("regression", "gradients"),
    )
