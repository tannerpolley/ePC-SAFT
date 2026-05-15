from __future__ import annotations

import json
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


def _assert_json_like(value):
    if isinstance(value, dict):
        for item in value.values():
            _assert_json_like(item)
    elif isinstance(value, list):
        for item in value:
            _assert_json_like(item)
    else:
        assert not isinstance(value, np.ndarray)


def _assert_no_backend_or_numerical_derivative_payload(value) -> None:
    payload = json.dumps(value, default=str).lower()
    assert "not_available" not in payload
    assert "numerical_derivative" not in payload


def _assert_methanol_cyclohexane_split(result: epcsaft.EquilibriumResult, feed: np.ndarray) -> None:
    assert result.split_detected is True
    assert result.stable is False
    assert result.phase_labels == ["liq1", "liq2"]
    assert len(result.phases) == 2
    beta = float(result.phases[1].phase_fraction)
    material = (1.0 - beta) * result.phases[0].composition + beta * result.phases[1].composition
    np.testing.assert_allclose(material, feed, atol=1.0e-10)
    assert result.diagnostics["fugacity_residual_norm"] < 1.0e-8
    assert result.diagnostics["material_balance_error"] < 1.0e-8
    assert result.diagnostics["phase_distance"] > 0.1
    assert result.diagnostics["nonlinear_solver"] == "native_derivative_free_nelder_mead"
    assert result.diagnostics["stability_analysis"] == "neutral_tpd"
    assert result.diagnostics["anti_trivial_solution_strategy"] == "phase_fraction_and_phase_distance_gate"
    assert result.diagnostics["derivative_backend"] == "not_applicable"
    assert result.diagnostics["derivative_status"] == "not_required"
    assert result.diagnostics["derivative_available"] is False
    assert result.diagnostics["jacobian_available"] is False
    _assert_json_like(result.to_dict())
    _assert_no_backend_or_numerical_derivative_payload(result.to_dict())


def test_methanol_cyclohexane_lle_flash_solves_seeded_phase_split() -> None:
    mix = _methanol_cyclohexane_mixture()
    feed, initial_phases = _methanol_cyclohexane_lle_benchmark()

    result = mix.equilibrium(
        kind="lle_flash",
        T=298.15,
        P=1.013e5,
        z=feed,
        backend="neutral_lle",
        initial_phases=initial_phases,
        options=epcsaft.EquilibriumOptions(max_iterations=240, tolerance=1.0e-10, damping=0.5),
    )

    _assert_methanol_cyclohexane_split(result, feed)
    assert result.diagnostics["seed_name"] == "user"
    assert result.diagnostics["attempt_count"] == 1


def test_lle_flash_without_initial_phases_solves_from_stability_seed() -> None:
    mix = _methanol_cyclohexane_mixture()
    feed, _initial_phases = _methanol_cyclohexane_lle_benchmark()

    result = mix.equilibrium(
        kind="lle_flash",
        T=298.15,
        P=1.013e5,
        z=feed,
        options=epcsaft.EquilibriumOptions(max_iterations=240, tolerance=1.0e-10, damping=0.5),
    )

    _assert_methanol_cyclohexane_split(result, feed)
    assert result.diagnostics["seed_name"].startswith("tpd_")
    assert result.diagnostics["attempt_count"] >= 1


def test_lle_flash_reports_no_split_for_identical_initial_phases() -> None:
    mix = _methanol_cyclohexane_mixture()
    feed, _initial_phases = _methanol_cyclohexane_lle_benchmark()

    result = mix.equilibrium(
        kind="lle_flash",
        T=298.15,
        P=1.013e5,
        z=feed,
        initial_phases={"liq1": feed, "liq2": feed, "phase_fraction": 0.5},
    )

    assert result.split_detected is False
    assert result.stable is False
    assert result.phase_labels == ["liq"]
    assert "no V2 LLE split" in result.diagnostics["message"]
    assert result.diagnostics["point_solver_split_detected"] is False
    assert "no V2 LLE split" in result.diagnostics["point_solver_message"]
    assert result.diagnostics["seed_name"] == "user"
    assert result.diagnostics["attempt_count"] == 1
    assert result.diagnostics["stability_analysis"] == "neutral_tpd"
    assert result.diagnostics["stability_stable"] is False
    assert np.isfinite(result.diagnostics["min_tpd"])


def test_lle_flash_can_skip_stability_precheck_for_debug_workflows() -> None:
    mix = _methanol_cyclohexane_mixture()
    feed, _initial_phases = _methanol_cyclohexane_lle_benchmark()

    result = mix.equilibrium(
        kind="lle_flash",
        T=298.15,
        P=1.013e5,
        z=feed,
        initial_phases={"liq1": feed, "liq2": feed, "phase_fraction": 0.5},
        options=epcsaft.EquilibriumOptions(stability_precheck=False),
    )

    assert result.split_detected is False
    assert result.stable is False
    assert result.diagnostics["stability_analysis"] == "not_run"
    assert result.diagnostics["stability_checked"] is False
    assert result.diagnostics["stability_stable"] is None
    assert "stability precheck skipped" in result.diagnostics["stability_message"]
    assert result.diagnostics["point_solver_split_detected"] is False
    assert "min_tpd" not in result.diagnostics


def test_lle_flash_phase_diagnostics_request_returns_clear_phase_details() -> None:
    mix = _methanol_cyclohexane_mixture()
    feed, initial_phases = _methanol_cyclohexane_lle_benchmark()

    result = mix.equilibrium(
        kind="lle_flash",
        T=298.15,
        P=1.013e5,
        z=feed,
        initial_phases=initial_phases,
        options=epcsaft.EquilibriumOptions(
            max_iterations=240,
            tolerance=1.0e-10,
            damping=0.5,
            include_phase_diagnostics=True,
        ),
    )

    _assert_methanol_cyclohexane_split(result, feed)
    for phase in result.phases:
        assert phase.diagnostics is not None
        assert phase.diagnostics["phase"] == phase.label
        assert phase.diagnostics["density"] > 0.0
        assert "fugacity_coefficient_terms" in phase.diagnostics


def test_equilibrium_options_expose_explicit_solver_backend_controls() -> None:
    option_fields = {field.name for field in fields(epcsaft.EquilibriumOptions)}

    assert "solver_backend" in option_fields
    assert "hessian_strategy" in option_fields
    assert "timeout_seconds" in option_fields
    assert "max_seed_attempts" in option_fields
    assert "max_density_failures" in option_fields
    assert "max_total_objective_evaluations" in option_fields
    assert epcsaft.EquilibriumOptions().solver_backend == "auto"
    assert epcsaft.EquilibriumOptions().hessian_strategy == "gauss_newton"
    assert epcsaft.EquilibriumOptions().timeout_seconds is None
    assert epcsaft.EquilibriumOptions().max_seed_attempts is None
    assert epcsaft.EquilibriumOptions().return_best_effort is False


def test_lle_flash_distinct_poor_seed_fails_loudly_without_forced_pass() -> None:
    mix = _methanol_cyclohexane_mixture()
    feed, _initial_phases = _methanol_cyclohexane_lle_benchmark()

    with pytest.raises(epcsaft.SolutionError, match="neutral LLE flash did not converge") as excinfo:
        mix.equilibrium(
            kind="lle_flash",
            T=298.15,
            P=1.013e5,
            z=feed,
            initial_phases={"liq1": [0.25, 0.75], "liq2": [0.65, 0.35], "phase_fraction": 0.5},
            options=epcsaft.EquilibriumOptions(max_iterations=80, tolerance=1.0e-6, damping=0.5),
        )

    message = str(excinfo.value)
    assert "best_seed=user" in message
    assert "maximum iterations" in message
    assert "numerical_derivative" not in message.lower()


@pytest.mark.parametrize(
    ("options", "match"),
    [
        (epcsaft.EquilibriumOptions(max_iterations=1.5), "max_iterations"),
        (epcsaft.EquilibriumOptions(max_iterations=True), "max_iterations"),
        (epcsaft.EquilibriumOptions(tolerance=float("nan")), "tolerance"),
        (epcsaft.EquilibriumOptions(damping=float("inf")), "damping"),
        (epcsaft.EquilibriumOptions(min_composition=float("nan")), "min_composition"),
        (epcsaft.EquilibriumOptions(include_phase_diagnostics="yes"), "include_phase_diagnostics"),
        (epcsaft.EquilibriumOptions(stability_precheck="yes"), "stability_precheck"),
        (epcsaft.EquilibriumOptions(solver_backend="cyipopt"), "solver_backend"),
        (epcsaft.EquilibriumOptions(hessian_strategy="exact"), "hessian_strategy"),
        (epcsaft.EquilibriumOptions(timeout_seconds=0.0), "timeout_seconds"),
        (epcsaft.EquilibriumOptions(timeout_seconds=float("nan")), "timeout_seconds"),
        (epcsaft.EquilibriumOptions(max_seed_attempts=0), "max_seed_attempts"),
        (epcsaft.EquilibriumOptions(max_seed_attempts=1.5), "max_seed_attempts"),
        (epcsaft.EquilibriumOptions(max_density_failures=0), "max_density_failures"),
        (epcsaft.EquilibriumOptions(max_total_objective_evaluations=0), "max_total_objective_evaluations"),
        (epcsaft.EquilibriumOptions(return_best_effort="yes"), "return_best_effort"),
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


def test_lle_flash_requested_ipopt_requires_cyipopt(monkeypatch) -> None:
    import epcsaft.ipopt_backend as ipopt_backend

    monkeypatch.setattr(ipopt_backend, "cyipopt_available", lambda: False)
    mix = _methanol_cyclohexane_mixture()
    feed, _initial_phases = _methanol_cyclohexane_lle_benchmark()

    with pytest.raises(epcsaft.InputError, match=r"cyipopt.*solver_backend='ipopt'"):
        mix.equilibrium(
            kind="lle_flash",
            T=298.15,
            P=1.013e5,
            z=feed,
            options=epcsaft.EquilibriumOptions(solver_backend="ipopt"),
        )


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
