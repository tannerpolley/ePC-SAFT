from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import epcsaft
from epcsaft import ePCSAFTMixture
from epcsaft.regression import _fit_pure_neutral_least_squares_internal
from scripts import plot_outputs
from tests.test_regression_api import _methane_like_records
from tests.test_regression_api import _minimal_neutral_metadata
from tests.test_runtime import _ionic_state
from tests.test_runtime import _neutral_state


def _hydrocarbon_basis_mixture() -> ePCSAFTMixture:
    params = {
        "m": np.asarray([1.0, 1.6069, 2.0020]),
        "s": np.asarray([3.7039, 3.5206, 3.6184]),
        "e": np.asarray([150.03, 191.42, 208.11]),
        "k_ij": np.asarray(
            [
                [0.0, 3.0e-4, 1.15e-2],
                [3.0e-4, 0.0, 5.10e-3],
                [1.15e-2, 5.10e-3, 0.0],
            ]
        ),
    }
    return ePCSAFTMixture.from_params(params, species=["Methane", "Ethane", "Propane"])


def _methanol_cyclohexane_mixture() -> ePCSAFTMixture:
    params = {
        "MW": np.asarray([32.042e-3, 84.147e-3]),
        "m": np.asarray([1.5255, 2.5303]),
        "s": np.asarray([3.2300, 3.8499]),
        "e": np.asarray([188.90, 278.11]),
        "e_assoc": np.asarray([2899.5, 0.0]),
        "vol_a": np.asarray([0.035176, 0.0]),
        "assoc_scheme": ["2B", None],
        "k_ij": np.asarray([[0.0, 0.051], [0.051, 0.0]]),
        "z": np.asarray([0.0, 0.0]),
        "dielc": np.asarray([33.05, 2.02]),
    }
    return ePCSAFTMixture.from_params(params, species=["Methanol", "Cyclohexane"])


def _assert_plot_with_data(path: Path) -> None:
    csv_path = path.parent / "data" / f"{path.stem}_plot_data.csv"
    assert path.exists()
    assert csv_path.exists()
    assert csv_path.stat().st_size > 0


def _save_comparison_plot(
    filename: str,
    title: str,
    labels: list[str],
    actual: np.ndarray,
    expected: np.ndarray,
    *,
    ylabel: str = "Value",
    relative_error: bool = True,
) -> Path:
    x = np.arange(len(labels), dtype=float)
    fig_height = 5.6 if relative_error else 4.4
    if relative_error:
        fig, axes = plt.subplots(2, 1, figsize=(max(7.0, len(labels) * 0.68), fig_height), height_ratios=[3, 1.35])
        ax, err_ax = axes
    else:
        fig, ax = plt.subplots(figsize=(max(7.0, len(labels) * 0.68), fig_height))
        err_ax = None

    ax.bar(x - 0.18, actual, width=0.36, label="Actual")
    ax.bar(x + 0.18, expected, width=0.36, label="Expected")
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xticks(x, labels, rotation=35, ha="right")
    ax.legend()
    ax.axhline(0.0, color="0.82", linewidth=0.8)

    if err_ax is not None:
        scale = np.maximum(np.abs(expected), 1.0e-30)
        err_ax.bar(x, (actual - expected) / scale, width=0.5, label="Relative error")
        err_ax.axhline(0.0, color="0.25", linewidth=0.8)
        err_ax.set_ylabel("Rel. err.")
        err_ax.set_xticks(x, labels, rotation=35, ha="right")

    output_path = plot_outputs.test_plot_path(__file__, filename)
    try:
        plot_outputs.save_plot_figure(fig, output_path, dpi=120)
    finally:
        plt.close(fig)
    _assert_plot_with_data(output_path)
    return output_path


def test_runtime_neutral_scalar_reference_comparison_plot() -> None:
    state, _species = _neutral_state()
    expected_values = {
        "rho": 14330.417110,
        "P": 1276374.1152948933,
        "Z": 0.04594621208078564,
        "ares": -3.54988545131505,
        "dadt": 0.03077401856781036,
        "hres": -15758.229958475444,
        "sres": -55.751451436621096,
        "gres": -2759.779056027235,
    }
    actual_values = {
        "rho": state.density(),
        "P": state.pressure(),
        "Z": state.compressibility_factor(),
        "ares": state.residual_helmholtz(),
        "dadt": state.temperature_derivative_residual_helmholtz(),
        "hres": state.residual_enthalpy(),
        "sres": state.residual_entropy(),
        "gres": state.residual_gibbs(),
    }

    labels = list(expected_values)
    _save_comparison_plot(
        "runtime_neutral_scalar_reference_comparison.png",
        "Neutral runtime scalar outputs vs pinned expected values",
        labels,
        np.asarray([actual_values[label] for label in labels], dtype=float),
        np.asarray([expected_values[label] for label in labels], dtype=float),
    )


def test_runtime_ionic_diagnostics_reference_comparison_plot() -> None:
    state, species = _ionic_state()
    component_activity = state.activity_coefficient(species=species)
    solvation = state.solvation_free_energy(species=species)
    expected_values = {
        "water gamma": 1.0000051724037697,
        "Na gamma": 0.9222113778654043,
        "Cl gamma": 0.9222258090371313,
        "mean gamma x": 0.9222185934230398,
        "mean gamma m": 0.9220341497043553,
        "Na gsolv": -475461.4260703414,
        "Cl gsolv": -489572.50284416083,
        "rho": 55344.274540081075,
        "rho mass": 997.1665703121223,
        "Z": 0.000728884077611683,
    }
    actual_values = {
        "water gamma": component_activity["water"],
        "Na gamma": component_activity["Na+"],
        "Cl gamma": component_activity["Cl-"],
        "mean gamma x": state.activity_coefficient(species=species, mean_ionic_form=True, basis="mole")["Na+Cl-"],
        "mean gamma m": state.activity_coefficient(species=species, mean_ionic_form=True, basis="molality")["Na+Cl-"],
        "Na gsolv": solvation["Na+"],
        "Cl gsolv": solvation["Cl-"],
        "rho": state.density(),
        "rho mass": state.mass_density(),
        "Z": state.compressibility_factor(),
    }

    labels = list(expected_values)
    _save_comparison_plot(
        "runtime_ionic_diagnostics_reference_comparison.png",
        "Ionic runtime diagnostics vs pinned expected values",
        labels,
        np.asarray([actual_values[label] for label in labels], dtype=float),
        np.asarray([expected_values[label] for label in labels], dtype=float),
    )


def test_native_branch_and_contribution_reference_comparison_plot() -> None:
    mix = _hydrocarbon_basis_mixture()
    composition = np.asarray([0.1, 0.3, 0.6])
    vapor_extreme = mix.state(T=600.0, x=composition, P=1.0, phase="vap")
    liquid_extreme = mix.state(T=220.0, x=composition, P=5.0e7, phase="liq")
    vapor_branch = mix.state(T=300.0, x=composition, P=1.0e3, phase="vap")
    liquid_branch = mix.state(T=300.0, x=composition, P=1.0e3, phase="liq")
    neutral_contract = mix.state(T=233.15, x=composition, rho=14330.417110)
    ares = neutral_contract.ares(return_contribution_terms=True)
    z_terms = neutral_contract.z(return_contribution_terms=True)

    expected_values = {
        "vap rho extreme": 2.0045400150430712e-4,
        "liq rho extreme": 16076.977238412512,
        "vap rho branch": 0.4009505832238275,
        "liq rho branch": 10700.137898056397,
        "ares total": -3.54988545131505,
        "Z total": 0.04594621208078564,
        "ares hc": 3.774229851214634,
        "ares disp": -7.324115302529684,
        "Z hc": 7.122473867439451,
        "Z disp": -8.076527655358666,
    }
    actual_values = {
        "vap rho extreme": vapor_extreme.density(),
        "liq rho extreme": liquid_extreme.density(),
        "vap rho branch": vapor_branch.density(),
        "liq rho branch": liquid_branch.density(),
        "ares total": ares["total"],
        "Z total": z_terms["total"],
        "ares hc": ares["terms"]["hc"],
        "ares disp": ares["terms"]["disp"],
        "Z hc": z_terms["terms"]["hc"],
        "Z disp": z_terms["terms"]["disp"],
    }

    labels = list(expected_values)
    _save_comparison_plot(
        "native_branch_contribution_reference_comparison.png",
        "Native branch and contribution outputs vs pinned expected values",
        labels,
        np.asarray([actual_values[label] for label in labels], dtype=float),
        np.asarray([expected_values[label] for label in labels], dtype=float),
    )


def test_vle_reference_comparison_plot() -> None:
    mix = _hydrocarbon_basis_mixture()
    result = mix.equilibrium(kind="tp_flash", T=220.0, P=1.0e5, z=[0.1, 0.3, 0.6])
    liquid, vapor = result.phases
    expected_values = {
        "beta liq": 0.0717673735624358,
        "beta vap": 0.9282326264375642,
        "liq methane": 0.0012963789214619132,
        "liq ethane": 0.06534426759935694,
        "liq propane": 0.9333593534791812,
        "vap methane": 0.10763138403472723,
        "vap ethane": 0.3181426779517501,
        "vap propane": 0.5742259380135226,
        "material tolerance": 1.0e-10,
        "fugacity tolerance": 1.0e-6,
    }
    actual_values = {
        "beta liq": liquid.phase_fraction,
        "beta vap": vapor.phase_fraction,
        "liq methane": liquid.composition[0],
        "liq ethane": liquid.composition[1],
        "liq propane": liquid.composition[2],
        "vap methane": vapor.composition[0],
        "vap ethane": vapor.composition[1],
        "vap propane": vapor.composition[2],
        "material tolerance": result.diagnostics["material_balance_error"],
        "fugacity tolerance": result.diagnostics["fugacity_residual_norm"],
    }

    labels = list(expected_values)
    _save_comparison_plot(
        "equilibrium_vle_reference_comparison.png",
        "Hydrocarbon VLE benchmark vs pinned expected values",
        labels,
        np.asarray([actual_values[label] for label in labels], dtype=float),
        np.asarray([expected_values[label] for label in labels], dtype=float),
        ylabel="Composition, fraction, or residual",
    )


def test_lle_reference_comparison_plot() -> None:
    mix = _methanol_cyclohexane_mixture()
    result = mix.equilibrium(
        kind="lle_flash",
        T=298.15,
        P=1.013e5,
        z=[0.45, 0.55],
        options=epcsaft.EquilibriumOptions(max_iterations=240, tolerance=1.0e-10, damping=0.5),
    )
    liq1, liq2 = result.phases
    expected_values = {
        "liq1 methanol": 0.1175783826,
        "liq1 cyclohexane": 0.8824216174,
        "liq2 methanol": 0.7985874309,
        "liq2 cyclohexane": 0.2014125691,
        "beta liq2": 0.4881309848,
        "material tolerance": 1.0e-10,
        "fugacity tolerance": 1.0e-9,
    }
    actual_values = {
        "liq1 methanol": liq1.composition[0],
        "liq1 cyclohexane": liq1.composition[1],
        "liq2 methanol": liq2.composition[0],
        "liq2 cyclohexane": liq2.composition[1],
        "beta liq2": liq2.phase_fraction,
        "material tolerance": result.diagnostics["material_balance_error"],
        "fugacity tolerance": result.diagnostics["fugacity_residual_norm"],
    }

    labels = list(expected_values)
    _save_comparison_plot(
        "equilibrium_lle_reference_comparison.png",
        "Methanol/cyclohexane LLE benchmark vs pinned expected values",
        labels,
        np.asarray([actual_values[label] for label in labels], dtype=float),
        np.asarray([expected_values[label] for label in labels], dtype=float),
        ylabel="Composition, fraction, or residual",
    )


def test_stability_reference_case_plot() -> None:
    lle_mix = _methanol_cyclohexane_mixture()
    unstable = lle_mix.equilibrium(
        kind="stability",
        T=298.15,
        P=1.013e5,
        z=[0.45, 0.55],
        parent_phase="liq",
        trial_phases=("liq",),
    )
    vle_mix = _hydrocarbon_basis_mixture()
    stable = vle_mix.equilibrium(
        kind="stability",
        T=300.0,
        P=1.0e5,
        z=[0.1, 0.3, 0.6],
        parent_phase="vap",
        trial_phases=("vap",),
    )

    labels = ["unstable min TPD", "unstable methanol", "unstable cyclohexane", "stable min TPD"]
    actual = np.asarray(
        [
            unstable.min_tpd,
            unstable.trial_composition[0],
            unstable.trial_composition[1],
            stable.min_tpd,
        ],
        dtype=float,
    )
    expected = np.asarray([-1.0e-4, 0.75, 0.25, -1.0e-6], dtype=float)
    _save_comparison_plot(
        "equilibrium_stability_reference_cases.png",
        "Neutral TPD stability cases vs decision thresholds",
        labels,
        actual,
        expected,
        ylabel="TPD or mole fraction",
        relative_error=False,
    )


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
    _save_comparison_plot(
        "regression_methane_parameter_reference_comparison.png",
        "Methane fitted parameters vs literature reference values",
        labels,
        np.asarray([actual_values[label] for label in labels], dtype=float),
        np.asarray([expected_values[label] for label in labels], dtype=float),
    )
