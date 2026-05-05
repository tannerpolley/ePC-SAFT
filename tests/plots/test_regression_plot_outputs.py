from __future__ import annotations

import numpy as np

import epcsaft
from tests.helpers.binary_regression_cases import ETHANOL_WATER_HELD_2012_KIJ
from tests.helpers.binary_regression_cases import ETHANOL_WATER_PAPER_PCSAFT_KIJ_100KPA
from tests.helpers.binary_regression_cases import ethanol_water_jced2021_vle_records
from epcsaft.regression import _fit_pure_neutral_least_squares_internal
from scripts import plot_outputs
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


def test_ethanol_water_binary_vle_real_data_kij_plot() -> None:
    result = epcsaft.fit_binary_pair(
        ethanol_water_jced2021_vle_records(smoke_only=True),
        ("Ethanol", "H2O"),
        dataset="2012_Held",
        initial_guess={"k_ij": ETHANOL_WATER_HELD_2012_KIJ},
        bounds={"k_ij": (-0.15, 0.10)},
        multistart=0,
    )
    assert result.success, result.message

    labels = ("2012 Held start", "native fit", "2021 paper PC-SAFT")
    values = np.asarray(
        [
            ETHANOL_WATER_HELD_2012_KIJ,
            result.fitted_values["k_ij"],
            ETHANOL_WATER_PAPER_PCSAFT_KIJ_100KPA,
        ],
        dtype=float,
    )
    colors = ("#6b7280", "#2563eb", "#111827")

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.bar(labels, values, color=colors)
    ax.axhline(0.0, color="0.25", linewidth=0.8)
    ax.set_ylabel(r"$k_{ij}$")
    ax.set_title("Ethanol/water native binary regression from real JCED VLE records")
    ax.text(
        0.02,
        0.95,
        f"native residual RMS: {result.metrics_by_term['binary_vle_fugacity_balance']:.4f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize="small",
    )
    ax.tick_params(axis="x", labelrotation=15)
    fig.tight_layout()

    out = plot_outputs.test_plot_path(
        __file__,
        "ethanol_water_binary_vle_real_data_kij.png",
        category=("regression", "binary_ethanol_water"),
    )
    plot_outputs.save_plot_figure(fig, out, dpi=160, svg_companion=True)
    plt.close(fig)
