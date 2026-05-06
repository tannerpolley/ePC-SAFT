from __future__ import annotations

import numpy as np

from tests.helpers.native_cases import _ionic_state as _native_ionic_state
from tests.helpers.native_cases import _neutral_state as _native_neutral_state
from analyses.package_plot_smokes.tests.plots.plot_helpers import hydrocarbon_basis_mixture
from analyses.package_plot_smokes.tests.plots.plot_helpers import save_comparison_plot
from analyses.package_plot_smokes.tests.plots.plot_helpers import save_parity_plot


def test_native_branch_and_contribution_reference_comparison_plot() -> None:
    mix = hydrocarbon_basis_mixture()
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
    save_comparison_plot(
        "native_branch_contribution_reference_comparison.png",
        "Native branch and contribution outputs vs pinned expected values",
        labels,
        np.asarray([actual_values[label] for label in labels], dtype=float),
        np.asarray([expected_values[label] for label in labels], dtype=float),
        category=("native", "branches"),
    )


def test_runtime_pressure_density_constructor_parity_plot() -> None:
    labels: list[str] = []
    pressure_values: list[float] = []
    density_values: list[float] = []
    for state_name, state_factory in (("neutral", _native_neutral_state), ("ionic", _native_ionic_state)):
        mix, _species, pressure, density, temperature, composition = state_factory()
        from_pressure = mix.state(T=temperature, x=composition, P=pressure, phase="liq")
        from_density = mix.state(T=temperature, x=composition, rho=density)
        for label, pressure_value, density_value in (
            (f"{state_name} rho", from_pressure.density(), from_density.density()),
            (f"{state_name} P", from_pressure.pressure(), from_density.pressure()),
            (f"{state_name} Z", from_pressure.z(), from_density.z()),
            (f"{state_name} ares", from_pressure.ares(), from_density.ares()),
        ):
            labels.append(label)
            pressure_values.append(float(pressure_value))
            density_values.append(float(density_value))

    save_parity_plot(
        "runtime_pressure_density_constructor_parity.png",
        "Pressure-constructed vs density-constructed state parity",
        labels,
        np.asarray(pressure_values, dtype=float),
        np.asarray(density_values, dtype=float),
        category=("native", "branches"),
        xlabel="Density constructor",
        ylabel="Pressure constructor",
    )


def test_native_temperature_derivative_finite_difference_parity_plot() -> None:
    mix, _species, _pressure, density, temperature, composition = _native_neutral_state()
    states = [
        ("rho state", mix.state(T=temperature, x=composition, rho=density)),
        ("vap branch", mix.state(T=300.0, x=composition, P=1.0e3, phase="vap")),
        ("liq branch", mix.state(T=300.0, x=composition, P=1.0e3, phase="liq")),
    ]
    labels: list[str] = []
    actual: list[float] = []
    expected: list[float] = []
    delta_t = 1.0e-3
    for label, state in states:
        plus = mix.state(T=state.T + delta_t, x=composition, rho=state.density(), phase="liq")
        minus = mix.state(T=state.T - delta_t, x=composition, rho=state.density(), phase="liq")
        finite_difference = (plus.ares() - minus.ares()) / (2.0 * delta_t)
        derivative = state.temperature_derivative_residual_helmholtz(return_contribution_terms=True)
        labels.append(label)
        actual.append(float(derivative["total"]))
        expected.append(float(finite_difference))

    save_parity_plot(
        "native_temperature_derivative_finite_difference_parity.png",
        "Native temperature derivative vs finite difference",
        labels,
        np.asarray(actual, dtype=float),
        np.asarray(expected, dtype=float),
        category=("native", "derivatives"),
    )


def test_native_composition_derivative_finite_difference_parity_plot() -> None:
    labels: list[str] = []
    actual: list[float] = []
    expected: list[float] = []
    for state_name, state_factory in (("neutral", _native_neutral_state), ("ionic", _native_ionic_state)):
        mix, _species, _pressure, density, temperature, composition = state_factory()
        state = mix.state(T=temperature, x=composition, rho=density)
        derivative = np.asarray(state.composition_derivative_residual_helmholtz()["total"], dtype=float)
        for i, j in ((0, 1), (1, 2), (0, 2)):
            delta_x = min(1.0e-6, 0.25 * float(composition[i]), 0.25 * float(composition[j]))
            plus = composition.copy()
            minus = composition.copy()
            plus[i] += delta_x
            plus[j] -= delta_x
            minus[i] -= delta_x
            minus[j] += delta_x
            finite_difference = (
                mix.state(T=temperature, x=plus, rho=density).ares()
                - mix.state(T=temperature, x=minus, rho=density).ares()
            ) / (2.0 * delta_x)
            labels.append(f"{state_name} d{i}-d{j}")
            actual.append(float(derivative[i] - derivative[j]))
            expected.append(float(finite_difference))

    save_parity_plot(
        "native_composition_derivative_finite_difference_parity.png",
        "Native composition derivative vs constrained finite difference",
        labels,
        np.asarray(actual, dtype=float),
        np.asarray(expected, dtype=float),
        category=("native", "derivatives"),
    )
