from __future__ import annotations

from dataclasses import fields

import numpy as np
import pytest

import epcsaft
from epcsaft import _core, ePCSAFTMixture


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


def _methanol_cyclohexane_lle_feed() -> np.ndarray:
    methanol_poor = np.asarray([0.05, 0.95], dtype=float)
    methanol_rich = np.asarray([0.85, 0.15], dtype=float)
    return 0.5 * methanol_poor + 0.5 * methanol_rich


def _assert_neutral_lle_route_pending(excinfo: pytest.ExceptionInfo[epcsaft.InputError]) -> None:
    message = str(excinfo.value)
    assert "lle_flash requires a native Ipopt equilibrium NLP route" in message
    assert "No package-owned alternate LLE solver is available" in message


def test_lle_flash_builds_one_native_route_request_before_ipopt_gate(monkeypatch: pytest.MonkeyPatch) -> None:
    mix = _methanol_cyclohexane_mixture()
    feed = _methanol_cyclohexane_lle_feed()
    calls: list[dict[str, object]] = []

    def fake_route(
        _native,
        temperature,
        pressure,
        feed_amounts,
        max_iterations,
        tolerance,
        material_tolerance,
        pressure_tolerance,
        chemical_potential_tolerance,
        phase_distance_tolerance,
    ):
        calls.append(
            {
                "temperature": temperature,
                "pressure": pressure,
                "feed_amounts": feed_amounts,
                "max_iterations": max_iterations,
                "tolerance": tolerance,
                "material_tolerance": material_tolerance,
                "pressure_tolerance": pressure_tolerance,
                "chemical_potential_tolerance": chemical_potential_tolerance,
                "phase_distance_tolerance": phase_distance_tolerance,
            }
        )
        return {
            "backend": "ipopt",
            "compiled": False,
            "ran": False,
            "accepted": False,
            "status": "ipopt_dependency_required",
            "postsolve": {"accepted": False},
        }

    monkeypatch.setattr(_core, "_native_neutral_lle_eos_route_result", fake_route)

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="lle_flash",
            T=298.15,
            P=1.013e5,
            z=feed,
            options=epcsaft.EquilibriumOptions(max_iterations=19, tolerance=3.0e-8),
        )

    _assert_neutral_lle_route_pending(excinfo)
    assert len(calls) == 1
    call = calls[0]
    assert call["feed_amounts"] == pytest.approx(feed.tolist())
    assert call["max_iterations"] == 19
    assert call["tolerance"] == pytest.approx(3.0e-8)


def test_lle_flash_converts_accepted_native_route_payload(monkeypatch: pytest.MonkeyPatch) -> None:
    mix = _methanol_cyclohexane_mixture()
    feed = _methanol_cyclohexane_lle_feed()
    route_amounts = [[0.03, 0.47], [0.42, 0.08]]
    route_volumes = [0.001, 0.002]

    def fake_route(
        _native,
        temperature,
        pressure,
        feed_amounts,
        max_iterations,
        tolerance,
        material_tolerance,
        pressure_tolerance,
        chemical_potential_tolerance,
        phase_distance_tolerance,
    ):
        assert temperature == pytest.approx(298.15)
        assert pressure == pytest.approx(1.013e5)
        assert feed_amounts == pytest.approx(feed.tolist())
        assert max_iterations > 0
        assert tolerance > 0.0
        assert material_tolerance > 0.0
        assert pressure_tolerance > 0.0
        assert chemical_potential_tolerance > 0.0
        assert phase_distance_tolerance > 0.0
        return {
            "backend": "ipopt",
            "compiled": True,
            "ran": True,
            "accepted": True,
            "status": "accepted",
            "phase_amounts": route_amounts,
            "phase_volumes": route_volumes,
            "postsolve": {"accepted": True},
        }

    def fake_result(
        _native,
        temperature,
        pressure,
        phase_amounts,
        phase_volumes,
        feed_amounts,
        material_tolerance,
        pressure_tolerance,
        chemical_potential_tolerance,
        phase_distance_tolerance,
    ):
        assert temperature == pytest.approx(298.15)
        assert pressure == pytest.approx(1.013e5)
        assert phase_amounts == route_amounts
        assert phase_volumes == route_volumes
        assert feed_amounts == pytest.approx(feed.tolist())
        assert material_tolerance > 0.0
        assert pressure_tolerance > 0.0
        assert chemical_potential_tolerance > 0.0
        assert phase_distance_tolerance > 0.0
        return {
            "accepted": True,
            "backend": "native_equilibrium_nlp",
            "problem_kind": "neutral_two_phase_eos",
            "stable": False,
            "split_detected": True,
            "phases": [
                {
                    "label": "phase_0",
                    "composition": [0.06, 0.94],
                    "density": 800.0,
                    "temperature": 298.15,
                    "pressure": 1.013e5,
                    "phase_fraction": 0.5,
                    "ln_fugacity_coefficient": [0.0, 0.1],
                    "fugacity_coefficient": [1.0, float(np.exp(0.1))],
                },
                {
                    "label": "phase_1",
                    "composition": [0.84, 0.16],
                    "density": 900.0,
                    "temperature": 298.15,
                    "pressure": 1.013e5,
                    "phase_fraction": 0.5,
                    "ln_fugacity_coefficient": [0.2, 0.0],
                    "fugacity_coefficient": [float(np.exp(0.2)), 1.0],
                },
            ],
        }

    monkeypatch.setattr(_core, "_native_neutral_lle_eos_route_result", fake_route)
    monkeypatch.setattr(_core, "_native_neutral_two_phase_eos_result", fake_result)

    result = mix.equilibrium(kind="lle_flash", T=298.15, P=1.013e5, z=feed)

    assert result.backend == "native_equilibrium_nlp"
    assert result.problem_kind == "neutral_lle"
    assert result.split_detected is True
    assert [phase.label for phase in result.phases] == ["liq1", "liq2"]
    assert result.phases[1].fugacity_coefficient == pytest.approx(np.exp(result.phases[1].ln_fugacity_coefficient))


def test_equilibrium_options_public_surface_is_current_fields() -> None:
    option_fields = {field.name for field in fields(epcsaft.EquilibriumOptions)}

    assert option_fields == {
        "max_iterations",
        "tolerance",
        "min_composition",
        "include_phase_diagnostics",
        "stability_precheck",
        "jacobian_backend",
        "solver_backend",
        "timeout_seconds",
    }
    assert epcsaft.EquilibriumOptions().solver_backend == "auto"
    assert epcsaft.EquilibriumOptions().timeout_seconds is None


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
        (epcsaft.EquilibriumOptions(solver_backend="new" + "ton"), "solver_backend"),
        (epcsaft.EquilibriumOptions(timeout_seconds=0.0), "timeout_seconds"),
        (epcsaft.EquilibriumOptions(timeout_seconds=float("nan")), "timeout_seconds"),
    ],
)
def test_lle_flash_rejects_invalid_options_through_public_api(options, match) -> None:
    mix = _methanol_cyclohexane_mixture()
    feed = _methanol_cyclohexane_lle_feed()

    with pytest.raises(epcsaft.InputError, match=match):
        mix.equilibrium(
            kind="lle_flash",
            T=298.15,
            P=1.013e5,
            z=feed,
            options=options,
        )


def test_lle_flash_rejects_unknown_option_dict_keys() -> None:
    mix = _methanol_cyclohexane_mixture()
    feed = _methanol_cyclohexane_lle_feed()

    with pytest.raises(epcsaft.InputError, match="unknown_solver_control"):
        mix.equilibrium(
            kind="lle_flash",
            T=298.15,
            P=1.013e5,
            z=feed,
            options={"unknown_solver_control": 1},
        )


def test_lle_flash_requested_ipopt_requires_native_ipopt_route() -> None:
    mix = _methanol_cyclohexane_mixture()
    feed = _methanol_cyclohexane_lle_feed()

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
    ],
)
def test_lle_flash_rejects_invalid_public_inputs(kwargs, match) -> None:
    mix = _methanol_cyclohexane_mixture()

    with pytest.raises(epcsaft.InputError, match=match):
        mix.equilibrium(**kwargs)


def test_lle_problem_has_no_public_phase_seed_field() -> None:
    field_name = "initial" + "_phases"

    assert field_name not in epcsaft.LLEProblem.__dataclass_fields__


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
