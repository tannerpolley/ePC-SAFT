from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import epcsaft
from epcsaft import ePCSAFTMixture
from epcsaft.regression import _debug_native_pure_neutral_objective
from epcsaft.regression import _fit_pure_neutral_least_squares_internal
from scripts import plot_outputs
from tests.helpers.native_cases import _ionic_state as _native_ionic_state
from tests.helpers.native_cases import _neutral_state as _native_neutral_state
from tests.helpers.regression_cases import _load_workbook_reference_rows
from tests.helpers.regression_cases import _methane_like_records
from tests.helpers.regression_cases import _minimal_neutral_metadata
from tests.helpers.regression_cases import _neutral_fixed_parameters
from tests.helpers.regression_cases import _real_saturation_records
from tests.helpers.runtime_cases import _ionic_state
from tests.helpers.runtime_cases import _neutral_state


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


def _save_parity_plot(
    filename: str,
    title: str,
    labels: list[str],
    actual: np.ndarray,
    expected: np.ndarray,
    *,
    xlabel: str = "Expected",
    ylabel: str = "Actual",
) -> Path:
    actual = np.asarray(actual, dtype=float)
    expected = np.asarray(expected, dtype=float)
    finite = np.isfinite(actual) & np.isfinite(expected)
    plot_actual = actual[finite]
    plot_expected = expected[finite]

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8), width_ratios=[2.2, 1.5])
    ax, err_ax = axes
    ax.scatter(plot_expected, plot_actual, label="Comparison points")
    if plot_expected.size:
        lo = float(min(np.min(plot_expected), np.min(plot_actual)))
        hi = float(max(np.max(plot_expected), np.max(plot_actual)))
        pad = max((hi - lo) * 0.08, 1.0e-12)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color="0.25", linewidth=1.0, label="Parity")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend()

    x = np.arange(len(labels), dtype=float)
    scale = np.maximum(np.abs(expected), 1.0e-30)
    err_ax.bar(x, (actual - expected) / scale, width=0.55)
    err_ax.axhline(0.0, color="0.25", linewidth=0.8)
    err_ax.set_title("Relative error")
    err_ax.set_xticks(x, labels, rotation=45, ha="right")

    output_path = plot_outputs.test_plot_path(__file__, filename)
    try:
        plot_outputs.save_plot_figure(fig, output_path, dpi=120)
    finally:
        plt.close(fig)
    _assert_plot_with_data(output_path)
    return output_path


def _save_contribution_closure_plot(filename: str, title: str, rows: list[dict[str, object]]) -> Path:
    labels = [str(row["label"]) for row in rows]
    term_names = sorted({name for row in rows for name in row["terms"]})
    x = np.arange(len(rows), dtype=float)

    fig, axes = plt.subplots(2, 1, figsize=(max(8.0, len(rows) * 0.72), 7.0), height_ratios=[3.0, 1.4])
    ax, err_ax = axes
    positive_bottom = np.zeros(len(rows), dtype=float)
    negative_bottom = np.zeros(len(rows), dtype=float)
    for term_name in term_names:
        values = np.asarray([float(row["terms"].get(term_name, 0.0)) for row in rows], dtype=float)
        bottoms = np.where(values >= 0.0, positive_bottom, negative_bottom)
        ax.bar(x, values, bottom=bottoms, width=0.58, label=term_name)
        positive_bottom += np.where(values >= 0.0, values, 0.0)
        negative_bottom += np.where(values < 0.0, values, 0.0)

    totals = np.asarray([float(row["total"]) for row in rows], dtype=float)
    term_sums = np.asarray([sum(float(value) for value in row["terms"].values()) for row in rows], dtype=float)
    ax.scatter(x, totals, marker="x", color="black", label="Reported total", zorder=5)
    ax.set_title(title)
    ax.set_ylabel("Contribution value")
    ax.set_xticks(x, labels, rotation=35, ha="right")
    ax.axhline(0.0, color="0.82", linewidth=0.8)
    ax.legend(ncol=min(3, max(1, len(term_names))), fontsize="small")

    scale = np.maximum(np.abs(totals), 1.0e-30)
    err_ax.bar(x, (term_sums - totals) / scale, width=0.55, label="Closure error")
    err_ax.axhline(0.0, color="0.25", linewidth=0.8)
    err_ax.set_ylabel("Rel. closure")
    err_ax.set_xticks(x, labels, rotation=35, ha="right")

    output_path = plot_outputs.test_plot_path(__file__, filename)
    try:
        plot_outputs.save_plot_figure(fig, output_path, dpi=120)
    finally:
        plt.close(fig)
    _assert_plot_with_data(output_path)
    return output_path


def _append_payload_rows(rows: list[dict[str, object]], prefix: str, payload: dict[str, object]) -> None:
    total = np.asarray(payload["total"], dtype=float)
    term_arrays = {key: np.asarray(value, dtype=float) for key, value in payload["terms"].items()}
    if total.ndim == 0:
        rows.append(
            {
                "label": prefix,
                "terms": {key: float(value) for key, value in payload["terms"].items()},
                "total": float(total),
            }
        )
        return
    for index, total_value in enumerate(total.tolist()):
        rows.append(
            {
                "label": f"{prefix}[{index}]",
                "terms": {key: float(value[index]) for key, value in term_arrays.items()},
                "total": float(total_value),
            }
        )


def _finite_difference_gradient_values() -> tuple[np.ndarray, np.ndarray]:
    theta = {"m": 1.05, "s": 3.68, "e": 151.0}

    def objective_at(m: float, s: float, e: float) -> float:
        debug = _debug_native_pure_neutral_objective(
            _methane_like_records(),
            "Methane",
            assoc_scheme="",
            fixed_parameters=_minimal_neutral_metadata(16.043e-3),
            initial_guess=theta,
            x={"m": m, "s": s, "e": e},
        )
        return float(debug["objective"])

    debug = _debug_native_pure_neutral_objective(
        _methane_like_records(),
        "Methane",
        assoc_scheme="",
        fixed_parameters=_minimal_neutral_metadata(16.043e-3),
        initial_guess=theta,
        x=theta,
    )
    exact = np.asarray(debug["gradient"], dtype=float)
    eps = np.asarray([1.0e-6, 1.0e-6, 1.0e-5], dtype=float)
    fd = np.empty(3, dtype=float)
    base = np.asarray([theta["m"], theta["s"], theta["e"]], dtype=float)
    for i in range(3):
        forward = base.copy()
        backward = base.copy()
        forward[i] += eps[i]
        backward[i] -= eps[i]
        fd[i] = (objective_at(*forward) - objective_at(*backward)) / (2.0 * eps[i])
    return exact, fd


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

    _save_comparison_plot(
        "regression_hydrocarbon_basis_parameter_comparison.png",
        "Hydrocarbon fitted parameters vs workbook/literature references",
        labels,
        np.asarray(actual, dtype=float),
        np.asarray(expected, dtype=float),
    )


def test_runtime_method_surface_neutral_reference_comparison_plot() -> None:
    state, _species = _neutral_state()
    mures = state.residual_chemical_potential()
    lnfug = state.fugacity_coefficient()
    fugcoef = state.fugacity_coefficient(natural_log=False)
    expected_values = {
        "rho": 14330.417110,
        "P": 1276374.1152948933,
        "Z": 0.04594621208078564,
        "ares": -3.54988545131505,
        "dadt": 0.03077401856781036,
        "hres": -15758.229958475444,
        "sres": -55.751451436621096,
        "gres": -2759.779056027235,
        "mures A": -1.1478687523834008,
        "mures B": -3.6543804288405415,
        "mures C": -5.488063725572939,
        "lnphi A": 1.9324151168689134,
        "lnphi B": -0.5740965595882255,
        "lnphi C": -2.407779856320623,
        "phi A": 6.906169322700795,
        "phi B": 0.5632134688356544,
        "phi C": 0.09001491894620331,
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
        "mures A": mures[0],
        "mures B": mures[1],
        "mures C": mures[2],
        "lnphi A": lnfug[0],
        "lnphi B": lnfug[1],
        "lnphi C": lnfug[2],
        "phi A": fugcoef[0],
        "phi B": fugcoef[1],
        "phi C": fugcoef[2],
    }
    labels = list(expected_values)
    _save_comparison_plot(
        "runtime_method_surface_neutral_reference_comparison.png",
        "Neutral public method surface vs pinned expected values",
        labels,
        np.asarray([actual_values[label] for label in labels], dtype=float),
        np.asarray([expected_values[label] for label in labels], dtype=float),
    )


def test_runtime_method_surface_ionic_reference_comparison_plot() -> None:
    state, species = _ionic_state()
    component_activity = state.activity_coefficient(species=species)
    solvation = state.solvation_free_energy(species=species)
    mures = state.residual_chemical_potential()
    lnfug = state.fugacity_coefficient()
    expected_values = {
        "P": 100000.0,
        "rho": 55344.274540081075,
        "rho mass": 997.1665703121223,
        "Z": 0.000728884077611683,
        "ares": -9.7214027218058,
        "dadt": 0.032388021640507005,
        "hres": -26415.160790583413,
        "sres": -59.523895812302186,
        "gres": -8668.111254145517,
        "mures water": -10.682420304620588,
        "mures Na": -199.10395742942775,
        "mures Cl": -204.79630395556683,
        "lnphi water": -3.458424439279275,
        "lnphi Na": -191.8799615776576,
        "lnphi Cl": -197.57230810636238,
        "gamma water": 1.0000051724037697,
        "gamma Na": 0.9222113778654043,
        "gamma Cl": 0.9222258090371313,
        "epsr": 78.075982,
        "osmotic": 0.9739566103279091,
        "gsolv Na": -475461.4260703414,
        "gsolv Cl": -489572.50284416083,
    }
    actual_values = {
        "P": state.pressure(),
        "rho": state.density(),
        "rho mass": state.mass_density(),
        "Z": state.compressibility_factor(),
        "ares": state.residual_helmholtz(),
        "dadt": state.temperature_derivative_residual_helmholtz(),
        "hres": state.residual_enthalpy(),
        "sres": state.residual_entropy(),
        "gres": state.residual_gibbs(),
        "mures water": mures[0],
        "mures Na": mures[1],
        "mures Cl": mures[2],
        "lnphi water": lnfug[0],
        "lnphi Na": lnfug[1],
        "lnphi Cl": lnfug[2],
        "gamma water": component_activity["water"],
        "gamma Na": component_activity["Na+"],
        "gamma Cl": component_activity["Cl-"],
        "epsr": state.relative_permittivity()[0],
        "osmotic": state.osmotic_coefficient()[0],
        "gsolv Na": solvation["Na+"],
        "gsolv Cl": solvation["Cl-"],
    }
    labels = list(expected_values)
    _save_comparison_plot(
        "runtime_method_surface_ionic_reference_comparison.png",
        "Ionic public method surface vs pinned expected values",
        labels,
        np.asarray([actual_values[label] for label in labels], dtype=float),
        np.asarray([expected_values[label] for label in labels], dtype=float),
    )


def test_runtime_alias_canonical_method_parity_plot() -> None:
    state, species = _ionic_state()
    labels = [
        "p",
        "rho",
        "rho_molar",
        "rho_mass",
        "z",
        "ares",
        "dadt",
        "hres",
        "sres",
        "gres",
        "mures water",
        "mures Na",
        "mures Cl",
        "fugcoef water",
        "fugcoef Na",
        "fugcoef Cl",
        "epsr mixture",
        "osmotic_coef",
        "gamma water",
        "gamma Na",
        "gamma Cl",
        "diag pressure",
        "gsolv Na",
        "gsolv Cl",
    ]
    alias_values = [
        state.p(),
        state.rho(),
        state.rho_molar(),
        state.rho_mass(),
        state.z(),
        state.ares(),
        state.dadt(),
        state.hres(),
        state.sres(),
        state.gres(),
        *state.mures().tolist(),
        *state.fugcoef().tolist(),
        state.epsr()[0],
        state.osmotic_coef()[0],
        *[state.gamma(species=species)[label] for label in species],
        state.diag(species=species)["pressure"],
        *[state.gsolv(species=species)[label] for label in ("Na+", "Cl-")],
    ]
    canonical_values = [
        state.pressure(),
        state.density(),
        state.molar_density(),
        state.mass_density(),
        state.compressibility_factor(),
        state.residual_helmholtz(),
        state.temperature_derivative_residual_helmholtz(),
        state.residual_enthalpy(),
        state.residual_entropy(),
        state.residual_gibbs(),
        *state.residual_chemical_potential().tolist(),
        *state.fugacity_coefficient().tolist(),
        state.relative_permittivity()[0],
        state.osmotic_coefficient()[0],
        *[state.activity_coefficient(species=species)[label] for label in species],
        state.state_diagnostics(species=species)["pressure"],
        *[state.solvation_free_energy(species=species)[label] for label in ("Na+", "Cl-")],
    ]
    _save_parity_plot(
        "runtime_alias_canonical_method_parity.png",
        "State aliases vs canonical public methods",
        labels,
        np.asarray(alias_values, dtype=float),
        np.asarray(canonical_values, dtype=float),
    )


def test_runtime_diagnostics_public_method_parity_plot() -> None:
    neutral_state, neutral_species = _neutral_state()
    ionic_state, ionic_species = _ionic_state()
    neutral_diag = neutral_state.state_diagnostics(species=neutral_species)
    ionic_diag = ionic_state.state_diagnostics(species=ionic_species)
    labels = [
        "neutral P",
        "neutral rho",
        "neutral Z",
        "neutral ares",
        "neutral mures A",
        "neutral mures B",
        "neutral mures C",
        "neutral phi A",
        "neutral phi B",
        "neutral phi C",
        "ionic P",
        "ionic rho",
        "ionic rho mass",
        "ionic Z",
        "ionic ares",
        "ionic gamma water",
        "ionic gamma Na",
        "ionic gamma Cl",
        "ionic mean gamma mole",
        "ionic mean gamma molality",
        "ionic osmotic",
        "ionic gsolv Na",
        "ionic gsolv Cl",
    ]
    diag_values = [
        neutral_diag["pressure"],
        neutral_diag["density"],
        neutral_diag["compressibility_factor"],
        neutral_diag["residual_helmholtz"],
        *np.asarray(neutral_diag["residual_chemical_potential"], dtype=float).tolist(),
        *np.asarray(neutral_diag["fugacity_coefficient"], dtype=float).tolist(),
        ionic_diag["pressure"],
        ionic_diag["density"],
        ionic_diag["mass_density"],
        ionic_diag["compressibility_factor"],
        ionic_diag["residual_helmholtz"],
        *[ionic_diag["activity_coefficient"][label] for label in ionic_species],
        ionic_diag["mean_ionic_activity_coefficient_mole"]["Na+Cl-"],
        ionic_diag["mean_ionic_activity_coefficient_molality"]["Na+Cl-"],
        np.asarray(ionic_diag["osmotic_coefficient"], dtype=float)[0],
        *[ionic_diag["solvation_free_energy"][label] for label in ("Na+", "Cl-")],
    ]
    public_values = [
        neutral_state.pressure(),
        neutral_state.density(),
        neutral_state.compressibility_factor(),
        neutral_state.residual_helmholtz(),
        *neutral_state.residual_chemical_potential().tolist(),
        *neutral_state.fugacity_coefficient(natural_log=False).tolist(),
        ionic_state.pressure(),
        ionic_state.density(),
        ionic_state.mass_density(),
        ionic_state.compressibility_factor(),
        ionic_state.residual_helmholtz(),
        *[ionic_state.activity_coefficient(species=ionic_species)[label] for label in ionic_species],
        ionic_state.activity_coefficient(species=ionic_species, mean_ionic_form=True, basis="mole")["Na+Cl-"],
        ionic_state.activity_coefficient(species=ionic_species, mean_ionic_form=True, basis="molality")["Na+Cl-"],
        ionic_state.osmotic_coefficient()[0],
        *[ionic_state.solvation_free_energy(species=ionic_species)[label] for label in ("Na+", "Cl-")],
    ]
    _save_parity_plot(
        "runtime_diagnostics_public_method_parity.png",
        "State diagnostics payload vs public methods",
        labels,
        np.asarray(diag_values, dtype=float),
        np.asarray(public_values, dtype=float),
    )


def test_runtime_neutral_contribution_closure_plot() -> None:
    state, _species = _neutral_state()
    rows: list[dict[str, object]] = []
    _append_payload_rows(rows, "ares", state.ares(return_contribution_terms=True))
    _append_payload_rows(rows, "Z", state.z(return_contribution_terms=True))
    _append_payload_rows(rows, "dadt", state.dadt(return_contribution_terms=True))
    _append_payload_rows(rows, "mures", state.mures(return_contribution_terms=True))
    _append_payload_rows(rows, "lnphi", state.fugcoef(return_contribution_terms=True))
    _save_contribution_closure_plot(
        "runtime_neutral_contribution_closure.png",
        "Neutral contribution closure: term sums vs reported totals",
        rows,
    )


def test_runtime_ionic_contribution_closure_plot() -> None:
    state, _species = _ionic_state()
    rows: list[dict[str, object]] = []
    _append_payload_rows(rows, "ares", state.ares(return_contribution_terms=True))
    _append_payload_rows(rows, "Z", state.z(return_contribution_terms=True))
    _append_payload_rows(rows, "dadt", state.dadt(return_contribution_terms=True))
    _append_payload_rows(rows, "mures", state.mures(return_contribution_terms=True))
    _append_payload_rows(rows, "lnphi", state.fugcoef(return_contribution_terms=True))
    _save_contribution_closure_plot(
        "runtime_ionic_contribution_closure.png",
        "Ionic contribution closure: term sums vs reported totals",
        rows,
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

    _save_parity_plot(
        "native_temperature_derivative_finite_difference_parity.png",
        "Native temperature derivative vs finite difference",
        labels,
        np.asarray(actual, dtype=float),
        np.asarray(expected, dtype=float),
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

    _save_parity_plot(
        "native_composition_derivative_finite_difference_parity.png",
        "Native composition derivative vs constrained finite difference",
        labels,
        np.asarray(actual, dtype=float),
        np.asarray(expected, dtype=float),
    )


def test_regression_gradient_finite_difference_parity_plot() -> None:
    actual, expected = _finite_difference_gradient_values()
    _save_parity_plot(
        "regression_gradient_finite_difference_parity.png",
        "Native regression gradient vs finite difference",
        ["m", "sigma", "epsilon"],
        actual,
        expected,
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

    _save_parity_plot(
        "runtime_pressure_density_constructor_parity.png",
        "Pressure-constructed vs density-constructed state parity",
        labels,
        np.asarray(pressure_values, dtype=float),
        np.asarray(density_values, dtype=float),
        xlabel="Density constructor",
        ylabel="Pressure constructor",
    )
