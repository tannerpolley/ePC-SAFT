from __future__ import annotations

from dataclasses import fields

import numpy as np
import pytest

import epcsaft
from epcsaft import ePCSAFTMixture


def _methanol_cyclohexane_mixture() -> ePCSAFTMixture:
    # Gross/Sadowski 2002 methanol parameters and kij; Gross/Sadowski 2001 cyclohexane parameters.
    params = {
        "MW": np.asarray([32.042e-3, 84.147e-3]),
        "m": np.asarray([1.5255, 2.5303]),
        "s": np.asarray([3.2300, 3.8499]),
        "e": np.asarray([188.90, 278.11]),
        "e_assoc": np.asarray([2899.5, 0.0]),
        "vol_a": np.asarray([0.035176, 0.0]),
        "assoc_scheme": ["2B", None],
        "k_ij": np.asarray(
            [
                [0.0, 0.051],
                [0.051, 0.0],
            ]
        ),
        "z": np.asarray([0.0, 0.0]),
        "dielc": np.asarray([33.05, 2.02]),
    }
    return ePCSAFTMixture.from_params(params, species=["Methanol", "Cyclohexane"])


def _methanol_cyclohexane_lle_benchmark() -> tuple[np.ndarray, dict[str, object]]:
    methanol_poor = np.asarray([0.05, 0.95], dtype=float)
    methanol_rich = np.asarray([0.85, 0.15], dtype=float)
    feed = 0.5 * methanol_poor + 0.5 * methanol_rich
    initial_phases = {
        "liq1": methanol_poor,
        "liq2": methanol_rich,
        "phase_fraction": 0.5,
    }
    return feed, initial_phases


def _assert_neutral_lle_route_pending(excinfo: pytest.ExceptionInfo[epcsaft.InputError]) -> None:
    message = str(excinfo.value)
    assert "lle_flash requires a native Ipopt equilibrium NLP route" in message
    assert "previous Ceres residual LLE route is disabled" in message


def test_methanol_cyclohexane_lle_flash_requires_native_ipopt_with_seed() -> None:
    mix = _methanol_cyclohexane_mixture()
    feed, initial_phases = _methanol_cyclohexane_lle_benchmark()

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="lle_flash",
            T=298.15,
            P=1.013e5,
            z=feed,
            backend="neutral_lle",
            initial_phases=initial_phases,
            options=epcsaft.EquilibriumOptions(max_iterations=240, tolerance=1.0e-10),
        )

    _assert_neutral_lle_route_pending(excinfo)


def test_lle_flash_without_initial_phases_requires_native_ipopt_after_validation() -> None:
    mix = _methanol_cyclohexane_mixture()
    feed, _initial_phases = _methanol_cyclohexane_lle_benchmark()

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="lle_flash",
            T=298.15,
            P=1.013e5,
            z=feed,
            options=epcsaft.EquilibriumOptions(max_iterations=240, tolerance=1.0e-10),
        )

    _assert_neutral_lle_route_pending(excinfo)


def test_equilibrium_options_expose_explicit_solver_backend_controls() -> None:
    option_fields = {field.name for field in fields(epcsaft.EquilibriumOptions)}

    assert "solver_backend" in option_fields
    assert "timeout_seconds" in option_fields
    assert "max_seed_attempts" in option_fields
    assert "max_density_failures" in option_fields
    assert "max_total_objective_evaluations" in option_fields
    assert epcsaft.EquilibriumOptions().solver_backend == "auto"
    assert epcsaft.EquilibriumOptions().timeout_seconds is None
    assert epcsaft.EquilibriumOptions().max_seed_attempts is None


@pytest.mark.parametrize(
    ("options", "match"),
    [
        (epcsaft.EquilibriumOptions(max_iterations=1.5), "max_iterations"),
        (epcsaft.EquilibriumOptions(max_iterations=True), "max_iterations"),
        (epcsaft.EquilibriumOptions(tolerance=float("nan")), "tolerance"),
        (epcsaft.EquilibriumOptions(min_composition=float("nan")), "min_composition"),
        (epcsaft.EquilibriumOptions(include_phase_diagnostics="yes"), "include_phase_diagnostics"),
        (epcsaft.EquilibriumOptions(stability_precheck="yes"), "stability_precheck"),
        (epcsaft.EquilibriumOptions(solver_backend="python_ipopt"), "solver_backend"),
        (epcsaft.EquilibriumOptions(solver_backend="newton"), "solver_backend"),
        (epcsaft.EquilibriumOptions(timeout_seconds=0.0), "timeout_seconds"),
        (epcsaft.EquilibriumOptions(timeout_seconds=float("nan")), "timeout_seconds"),
        (epcsaft.EquilibriumOptions(max_seed_attempts=0), "max_seed_attempts"),
        (epcsaft.EquilibriumOptions(max_seed_attempts=1.5), "max_seed_attempts"),
        (epcsaft.EquilibriumOptions(max_density_failures=0), "max_density_failures"),
        (epcsaft.EquilibriumOptions(max_total_objective_evaluations=0), "max_total_objective_evaluations"),
    ],
)
def test_lle_flash_rejects_invalid_options_through_public_api(options, match) -> None:
    mix = _methanol_cyclohexane_mixture()
    feed, _initial_phases = _methanol_cyclohexane_lle_benchmark()

    with pytest.raises(epcsaft.InputError, match=match):
        mix.equilibrium(
            kind="lle_flash",
            T=298.15,
            P=1.013e5,
            z=feed,
            initial_phases={"liq1": feed, "liq2": feed, "phase_fraction": 0.5},
            options=options,
        )


def test_lle_flash_requested_ipopt_requires_native_ipopt_route() -> None:
    mix = _methanol_cyclohexane_mixture()
    feed, _initial_phases = _methanol_cyclohexane_lle_benchmark()

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="lle_flash",
            T=298.15,
            P=1.013e5,
            z=feed,
            options=epcsaft.EquilibriumOptions(solver_backend="ipopt"),
        )

    _assert_neutral_lle_route_pending(excinfo)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"kind": "lle_flash", "P": 1.0e5, "z": [0.5, 0.5]}, "T"),
        ({"kind": "lle_flash", "T": 298.15, "z": [0.5, 0.5]}, "P"),
        ({"kind": "lle_flash", "T": 298.15, "P": 1.0e5}, "z"),
        ({"kind": "lle_flash", "T": 298.15, "P": 1.0e5, "z": [1.0]}, "length"),
        ({"kind": "lle_flash", "T": 298.15, "P": 1.0e5, "z": [0.5, -0.5]}, "non-negative"),
        (
            {
                "kind": "lle_flash",
                "T": 298.15,
                "P": 1.0e5,
                "z": [0.5, 0.5],
                "initial_phases": {"liq1": [0.5, 0.5], "phase_fraction": 0.5},
            },
            "initial_phases",
        ),
        (
            {
                "kind": "lle_flash",
                "T": 298.15,
                "P": 1.0e5,
                "z": [0.5, 0.5],
                "initial_phases": {"liq1": [0.5, 0.5], "liq2": [0.2, 0.8], "phase_fraction": 1.2},
            },
            "phase_fraction",
        ),
    ],
)
def test_lle_flash_rejects_invalid_public_inputs(kwargs, match) -> None:
    mix = _methanol_cyclohexane_mixture()

    with pytest.raises(epcsaft.InputError, match=match):
        mix.equilibrium(**kwargs)


def test_lle_flash_rejects_ionic_mixtures_for_v2() -> None:
    params = {
        "m": np.asarray([1.2047, 1.0, 1.0]),
        "s": np.asarray([2.7927, 2.8232, 2.7560]),
        "e": np.asarray([353.95, 230.0, 170.0]),
        "z": np.asarray([0.0, 1.0, -1.0]),
        "dielc": np.asarray([78.09, 8.0, 8.0]),
    }
    mix = ePCSAFTMixture.from_params(params, species=["water", "Na+", "Cl-"])

    with pytest.raises(epcsaft.InputError, match="ion-containing"):
        mix.equilibrium(kind="lle_flash", T=298.15, P=1.0e5, z=[0.9998, 1.0e-4, 1.0e-4])
