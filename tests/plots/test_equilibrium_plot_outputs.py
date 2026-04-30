from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go

import epcsaft
from epcsaft.equilibrium_core.thermo_diagnostics import evaluate_khudaida_solver_gate
from scripts import plot_outputs
from tests.plots.plot_helpers import assert_plot_with_data
from tests.plots.plot_helpers import hydrocarbon_basis_mixture
from tests.plots.plot_helpers import methanol_cyclohexane_mixture
from tests.plots.plot_helpers import save_comparison_plot
from tests.plots.plot_helpers import save_plotly_html


def _ascani_electrolyte_lle_result() -> tuple[epcsaft.ePCSAFTMixture, epcsaft.EquilibriumResult]:
    mix = epcsaft.ePCSAFTMixture.from_dataset(
        "2022_Ascani",
        ["H2O", "Butanol", "Na+", "Cl-"],
        [0.55, 0.40, 0.025, 0.025],
        298.15,
    )
    result = mix.equilibrium(
        kind="electrolyte_lle",
        T=298.15,
        P=1.013e5,
        solvent_feed={"H2O": 0.58, "Butanol": 0.42},
        salt_molality={"NaCl": 1.0},
        options=epcsaft.EquilibriumOptions(include_phase_diagnostics=True),
    )
    return mix, result


def test_equilibrium_vle_composition_plot_is_written_to_gallery() -> None:
    mix = hydrocarbon_basis_mixture()
    result = mix.equilibrium(kind="tp_flash", T=220.0, P=1.0e5, z=[0.1, 0.3, 0.6])
    liquid, vapor = result.phases
    species = np.asarray(mix.species)
    x = np.arange(species.size)

    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.bar(x - 0.18, liquid.composition, width=0.36, label="Liquid")
    ax.bar(x + 0.18, vapor.composition, width=0.36, label="Vapor")
    ax.set_xticks(x, species)
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel(r"Mole fraction, $x_i$ or $y_i$")
    ax.set_title("Hydrocarbon TP flash phase compositions")
    ax.legend()

    output_path = plot_outputs.test_plot_path(
        __file__,
        "equilibrium_vle_compositions.png",
        category=("equilibrium", "vle"),
    )
    try:
        plot_outputs.save_plot_figure(fig, output_path, dpi=120, svg_companion=True)
        interactive = go.Figure()
        interactive.add_trace(go.Bar(x=species.tolist(), y=liquid.composition.tolist(), name="Liquid"))
        interactive.add_trace(go.Bar(x=species.tolist(), y=vapor.composition.tolist(), name="Vapor"))
        interactive.update_layout(
            title="Hydrocarbon TP flash phase compositions",
            barmode="group",
            xaxis_title="Species",
            yaxis_title=r"Mole fraction, $x_i$ or $y_i$",
            yaxis_range=[0.0, 1.0],
        )
        save_plotly_html(interactive, Path(output_path))
    finally:
        plt.close(fig)

    assert_plot_with_data(Path(output_path))


def test_equilibrium_vle_reference_comparison_plot() -> None:
    mix = hydrocarbon_basis_mixture()
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
    save_comparison_plot(
        "equilibrium_vle_reference_comparison.png",
        "Hydrocarbon VLE benchmark vs pinned expected values",
        labels,
        np.asarray([actual_values[label] for label in labels], dtype=float),
        np.asarray([expected_values[label] for label in labels], dtype=float),
        category=("equilibrium", "vle"),
        ylabel="Composition, fraction, or residual",
    )


def test_equilibrium_vle_residual_closure_plot() -> None:
    mix = hydrocarbon_basis_mixture()
    result = mix.equilibrium(kind="tp_flash", T=220.0, P=1.0e5, z=[0.1, 0.3, 0.6])
    labels = ["material balance", "fugacity residual"]
    save_comparison_plot(
        "equilibrium_vle_residual_closure.png",
        "Hydrocarbon VLE residual closure vs test tolerances",
        labels,
        np.asarray(
            [
                result.diagnostics["material_balance_error"],
                result.diagnostics["fugacity_residual_norm"],
            ],
            dtype=float,
        ),
        np.asarray([1.0e-10, 1.0e-6], dtype=float),
        category=("equilibrium", "vle"),
        ylabel="Residual norm",
        relative_error=False,
    )


def test_equilibrium_lle_tie_line_plot_is_written_to_gallery() -> None:
    mix = methanol_cyclohexane_mixture()
    feed = np.asarray([0.45, 0.55], dtype=float)
    result = mix.equilibrium(
        kind="lle_flash",
        T=298.15,
        P=1.013e5,
        z=feed,
        options=epcsaft.EquilibriumOptions(max_iterations=240, tolerance=1.0e-10, damping=0.5),
    )
    liq1, liq2 = result.phases

    fig, ax = plt.subplots(figsize=(6.6, 4.0))
    ax.plot([liq1.composition[0], liq2.composition[0]], [0.0, 0.0], label="Tie line")
    ax.scatter([feed[0]], [0.0], label="Feed")
    ax.scatter([liq1.composition[0]], [0.0], label="Liquid 1")
    ax.scatter([liq2.composition[0]], [0.0], label="Liquid 2")
    ax.set_xlim(0.0, 1.0)
    ax.set_yticks([])
    ax.set_xlabel(r"Methanol mole fraction, $x_{\mathrm{MeOH}}$")
    ax.set_title("Methanol/cyclohexane LLE tie line")
    ax.legend(loc="upper center", ncol=4)

    output_path = plot_outputs.test_plot_path(
        __file__,
        "equilibrium_lle_tie_line.png",
        category=("equilibrium", "lle"),
    )
    try:
        plot_outputs.save_plot_figure(fig, output_path, dpi=120, svg_companion=True)
        interactive = go.Figure()
        interactive.add_trace(
            go.Scatter(
                x=[liq1.composition[0], liq2.composition[0]],
                y=[0.0, 0.0],
                mode="lines+markers",
                name="Tie line",
                text=["Liquid 1", "Liquid 2"],
                hovertemplate=r"%{text}<br>$x_{\mathrm{MeOH}}$=%{x:.6g}<extra></extra>",
            )
        )
        interactive.add_trace(go.Scatter(x=[feed[0]], y=[0.0], mode="markers", name="Feed"))
        interactive.update_layout(
            title="Methanol/cyclohexane LLE tie line",
            xaxis_title=r"Methanol mole fraction, $x_{\mathrm{MeOH}}$",
            xaxis_range=[0.0, 1.0],
            yaxis_visible=False,
        )
        save_plotly_html(interactive, Path(output_path))
    finally:
        plt.close(fig)

    assert_plot_with_data(Path(output_path))


def test_equilibrium_lle_reference_comparison_plot() -> None:
    mix = methanol_cyclohexane_mixture()
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
    save_comparison_plot(
        "equilibrium_lle_reference_comparison.png",
        "Methanol/cyclohexane LLE benchmark vs pinned expected values",
        labels,
        np.asarray([actual_values[label] for label in labels], dtype=float),
        np.asarray([expected_values[label] for label in labels], dtype=float),
        category=("equilibrium", "lle"),
        ylabel="Composition, fraction, or residual",
    )


def test_equilibrium_lle_residual_closure_plot() -> None:
    mix = methanol_cyclohexane_mixture()
    result = mix.equilibrium(
        kind="lle_flash",
        T=298.15,
        P=1.013e5,
        z=[0.45, 0.55],
        options=epcsaft.EquilibriumOptions(max_iterations=240, tolerance=1.0e-10, damping=0.5),
    )
    labels = ["material balance", "fugacity residual"]
    save_comparison_plot(
        "equilibrium_lle_residual_closure.png",
        "Methanol/cyclohexane LLE residual closure vs test tolerances",
        labels,
        np.asarray(
            [
                result.diagnostics["material_balance_error"],
                result.diagnostics["fugacity_residual_norm"],
            ],
            dtype=float,
        ),
        np.asarray([1.0e-10, 1.0e-9], dtype=float),
        category=("equilibrium", "lle"),
        ylabel="Residual norm",
        relative_error=False,
    )


def test_equilibrium_electrolyte_lle_phase_composition_plot() -> None:
    mix, result = _ascani_electrolyte_lle_result()
    diagnostics = result.diagnostics
    phases = {phase.label: phase for phase in result.phases}
    species = np.asarray(mix.species)
    x = np.arange(species.size)

    assert result.split_detected is True
    assert diagnostics["acceptance_gate"] == "predictive_nonlinear_solve"
    assert diagnostics["phase_equilibrium_model"] == "electrolyte_lle_v4_charge_constrained_solve"
    assert phases["aq"].composition[0] > phases["org"].composition[0]
    assert phases["org"].composition[1] > phases["aq"].composition[1]

    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    ax.bar(x - 0.24, diagnostics["feed_composition"], width=0.24, label="Feed")
    ax.bar(x, phases["aq"].composition, width=0.24, label="Aqueous phase")
    ax.bar(x + 0.24, phases["org"].composition, width=0.24, label="Organic phase")
    ax.set_xticks(x, species)
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel(r"Mole fraction, $x_i$")
    ax.set_title("V4 electrolyte LLE phase compositions")
    ax.legend(loc="upper right")

    output_path = plot_outputs.test_plot_path(
        __file__,
        "equilibrium_electrolyte_lle_phase_compositions.png",
        category=("equilibrium", "electrolyte_lle"),
    )
    try:
        plot_outputs.save_plot_figure(fig, output_path, dpi=120, svg_companion=True)
        interactive = go.Figure()
        interactive.add_trace(go.Bar(x=species.tolist(), y=diagnostics["feed_composition"], name="Feed"))
        interactive.add_trace(go.Bar(x=species.tolist(), y=phases["aq"].composition.tolist(), name="Aqueous phase"))
        interactive.add_trace(go.Bar(x=species.tolist(), y=phases["org"].composition.tolist(), name="Organic phase"))
        interactive.update_layout(
            title="V4 electrolyte LLE phase compositions",
            barmode="group",
            xaxis_title="Species",
            yaxis_title=r"Mole fraction, $x_i$",
            yaxis_range=[0.0, 1.0],
        )
        save_plotly_html(interactive, Path(output_path))
    finally:
        plt.close(fig)

    assert_plot_with_data(Path(output_path))


def test_equilibrium_electrolyte_lle_residual_closure_plot() -> None:
    _mix, result = _ascani_electrolyte_lle_result()
    diagnostics = result.diagnostics

    assert diagnostics["acceptance_gate"] == "predictive_nonlinear_solve"
    labels = ["solver residual", "material balance", "charge balance", "Gibbs favored"]
    actual = np.asarray(
        [
            diagnostics["solver_residual_norm"],
            diagnostics["material_balance_error"],
            diagnostics["charge_balance_error"],
            max(diagnostics["gibbs_delta"], 0.0),
        ],
        dtype=float,
    )
    tolerances = np.asarray([1.0e-6, 1.0e-10, 1.0e-8, 1.0e-15], dtype=float)
    save_comparison_plot(
        "equilibrium_electrolyte_lle_residual_closure.png",
        "V4 electrolyte LLE residual closure vs acceptance tolerances",
        labels,
        actual,
        tolerances,
        category=("equilibrium", "electrolyte_lle"),
        ylabel="Residual norm or clipped Gibbs delta",
        relative_error=False,
    )


def test_equilibrium_electrolyte_khudaida_seeded_solver_gate_plot() -> None:
    diagnostics = evaluate_khudaida_solver_gate(figure=2, tie_line=1, source="package", seeded=True)
    solver = diagnostics["solver_diagnostics"]

    assert diagnostics["solver_outcome"] == "accepted"
    assert diagnostics["decision"] == "solver_accepts_package_fixed_tieline_feed"
    assert solver["solver_seed_name"] == "initial_phases"
    labels = ["fixed tie-line residual", "solver residual", "material balance", "charge balance", "Gibbs favored"]
    actual = np.asarray(
        [
            diagnostics["fixed_phase_residual_norm"],
            solver["solver_residual_norm"],
            solver["material_balance_error"],
            solver["charge_balance_error"],
            max(solver["gibbs_delta"], 0.0),
        ],
        dtype=float,
    )
    tolerances = np.asarray([1.0e-6, 1.0e-6, 1.0e-10, 1.0e-8, 1.0e-15], dtype=float)
    save_comparison_plot(
        "equilibrium_electrolyte_khudaida_seeded_solver_gate.png",
        "V4 Khudaida seeded electrolyte LLE solver gate",
        labels,
        actual,
        tolerances,
        category=("equilibrium", "electrolyte_lle"),
        ylabel="Residual norm or clipped Gibbs delta",
        relative_error=False,
    )


def test_equilibrium_stability_reference_case_plot() -> None:
    lle_mix = methanol_cyclohexane_mixture()
    unstable = lle_mix.equilibrium(
        kind="stability",
        T=298.15,
        P=1.013e5,
        z=[0.45, 0.55],
        parent_phase="liq",
        trial_phases=("liq",),
    )
    vle_mix = hydrocarbon_basis_mixture()
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
    save_comparison_plot(
        "equilibrium_stability_reference_cases.png",
        "Neutral TPD stability cases vs decision thresholds",
        labels,
        actual,
        expected,
        category=("equilibrium", "stability"),
        ylabel="TPD or mole fraction",
        relative_error=False,
    )
