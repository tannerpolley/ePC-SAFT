from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import epcsaft
from epcsaft import ePCSAFTMixture
from scripts import plot_outputs


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


def test_equilibrium_vle_composition_plot_is_written_to_gallery() -> None:
    mix = _hydrocarbon_basis_mixture()
    result = mix.equilibrium(kind="tp_flash", T=220.0, P=1.0e5, z=[0.1, 0.3, 0.6])
    liquid, vapor = result.phases
    species = np.asarray(mix.species)
    x = np.arange(species.size)

    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.bar(x - 0.18, liquid.composition, width=0.36, label="Liquid")
    ax.bar(x + 0.18, vapor.composition, width=0.36, label="Vapor")
    ax.set_xticks(x, species)
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel("Mole fraction")
    ax.set_title("Hydrocarbon TP flash phase compositions")
    ax.legend()

    output_path = plot_outputs.test_plot_path(__file__, "equilibrium_vle_compositions.png")
    try:
        plot_outputs.save_plot_figure(fig, output_path, dpi=120)
    finally:
        plt.close(fig)

    _assert_plot_with_data(output_path)


def test_equilibrium_lle_tie_line_plot_is_written_to_gallery() -> None:
    mix = _methanol_cyclohexane_mixture()
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
    ax.set_xlabel("Methanol mole fraction")
    ax.set_title("Methanol/cyclohexane LLE tie line")
    ax.legend(loc="upper center", ncol=4)

    output_path = plot_outputs.test_plot_path(__file__, "equilibrium_lle_tie_line.png")
    try:
        plot_outputs.save_plot_figure(fig, output_path, dpi=120)
    finally:
        plt.close(fig)

    _assert_plot_with_data(output_path)
