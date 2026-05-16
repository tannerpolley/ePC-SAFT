"""Native runtime contracts for pressure-vs-density closure and contribution terms."""

from __future__ import annotations

import json

import numpy as np
import pytest

from epcsaft import InputError, SolutionError
from tests.helpers.native_cases import _ionic_state, _neutral_state


def _assert_close_terms(observed: dict[str, float], expected: dict[str, float]) -> None:
    assert set(observed) == set(expected)
    for key, value in expected.items():
        assert observed[key] == pytest.approx(value, rel=1e-10, abs=1e-12)

def _assert_finite_mapping_values(values: dict[str, float]) -> None:
    assert values
    assert all(np.isfinite(float(value)) for value in values.values())

def test_pressure_based_and_density_based_states_match_for_neutral_system() -> None:
    mix, _, pressure, density, temperature, composition = _neutral_state()
    from_rho = mix.state(T=temperature, x=composition, rho=density)
    from_p = mix.state(T=temperature, x=composition, P=pressure, phase="liq")

    assert from_p.density() == pytest.approx(from_rho.density())
    assert from_p.pressure() == pytest.approx(from_rho.pressure())
    assert from_p.z() == pytest.approx(from_rho.z())
    assert from_p.ares() == pytest.approx(from_rho.ares())

def test_pressure_based_and_density_based_states_match_for_ionic_system() -> None:
    mix, _, pressure, density, temperature, composition = _ionic_state()
    from_rho = mix.state(T=temperature, x=composition, rho=density)
    from_p = mix.state(T=temperature, x=composition, P=pressure, phase="liq")

    assert from_p.density() == pytest.approx(from_rho.density())
    assert from_p.pressure() == pytest.approx(from_rho.pressure())
    assert from_p.z() == pytest.approx(from_rho.z())
    assert from_p.ares() == pytest.approx(from_rho.ares())

def test_pressure_density_edge_cases_cover_vapor_and_liquid_extremes() -> None:
    mix, _, _, _, _, composition = _neutral_state()

    vapor = mix.state(T=600.0, x=composition, P=1.0, phase="vap")
    liquid = mix.state(T=220.0, x=composition, P=5.0e7, phase="liq")

    assert vapor.phase == 1
    assert liquid.phase == 0
    assert vapor.density() == pytest.approx(2.0045400150430712e-4, rel=1e-10)
    assert liquid.density() == pytest.approx(16076.977238412512, rel=1e-10)
    assert vapor.pressure() == pytest.approx(1.0)
    assert liquid.pressure() == pytest.approx(5.0e7)
    assert np.isfinite(vapor.z())
    assert np.isfinite(liquid.z())

def test_pressure_density_phase_branches_do_not_cross_at_two_phase_like_state() -> None:
    mix, _, _, _, _, composition = _neutral_state()
    temperature = 300.0
    pressure = 1.0e3

    vapor = mix.state(T=temperature, x=composition, P=pressure, phase="vap")
    liquid = mix.state(T=temperature, x=composition, P=pressure, phase="liq")

    assert vapor.phase == 1
    assert liquid.phase == 0
    assert vapor.pressure() == pytest.approx(pressure)
    assert liquid.pressure() == pytest.approx(pressure)
    assert vapor.density() == pytest.approx(0.4009505832238275, rel=1e-10)
    assert liquid.density() == pytest.approx(10700.137898056397, rel=1e-10)
    assert liquid.density() / vapor.density() > 1.0e4

def test_ionic_high_pressure_liquid_density_branch_remains_stable() -> None:
    mix, species, _, _, _, composition = _ionic_state()

    state = mix.state(T=320.0, x=composition, P=5.0e7, phase="liq")
    gamma = state.activity_coefficient(species=species, mean_ionic_form=True, basis="molality")

    assert state.phase == 0
    assert state.pressure() == pytest.approx(5.0e7)
    assert state.density() == pytest.approx(55177.833730914346, rel=1e-10)
    assert state.z() == pytest.approx(0.3405817367443399, rel=1e-10)
    assert gamma == {"Na+Cl-": pytest.approx(0.9295936238832602, rel=1e-10)}

def test_pressure_density_failure_reports_state_context_and_native_outcome() -> None:
    mix, _, _, _, temperature, composition = _ionic_state()

    with pytest.raises(SolutionError) as excinfo:
        mix.state(T=temperature, x=composition, P=1.0e12, phase="liq")

    message = str(excinfo.value)
    assert "pressure-based state solve failed" in message
    assert "T=298.15" in message
    assert "P=1000000000000.0" in message
    assert "phase=liq" in message
    assert "ncomp=3" in message
    assert "No continuous density root brackets were found" in message
    diagnostics = excinfo.value.diagnostics
    assert diagnostics["density_failure_count"] == 1
    assert diagnostics["density_validity_gate"] == "failed"
    assert diagnostics["density_best_candidate_refinement_used"] in {True, False}
    context = diagnostics["density_failure_contexts"][0]
    assert context["phase_label"] == "state"
    assert context["phase_kind"] == "liq"
    assert context["T"] == pytest.approx(temperature)
    assert context["P"] == pytest.approx(1.0e12)
    assert context["composition"] == pytest.approx(composition.tolist())
    assert context["scan_point_count"] > 0
    assert context["finite_point_count"] >= 0
    assert context["coarse_bracket_count"] >= 0
    assert context["refined_bracket_count"] >= 0
    assert context["candidate_root_count"] >= 0
    assert "best_near_root_pressure_error" in context
    assert "gres" in context
    assert context["rejection_reason"]
    json.dumps(diagnostics, allow_nan=False)

def test_pressure_density_invalid_inputs_fail_before_native_density_search() -> None:
    mix, _, pressure, _, temperature, composition = _neutral_state()

    with pytest.raises(InputError, match="State composition length"):
        mix.state(T=temperature, x=composition[:-1], P=pressure, phase="liq")
    with pytest.raises(InputError, match="phase must be"):
        mix.state(T=temperature, x=composition, P=pressure, phase="solid")
    with pytest.raises(InputError, match="Provide exactly one of P or rho"):
        mix.state(T=temperature, x=composition, P=pressure, rho=1.0, phase="liq")

def test_density_guess_uses_native_warm_start_diagnostics() -> None:
    mix, _, pressure, _, temperature, composition = _ionic_state()
    reference = mix.state(T=temperature, x=composition, P=pressure, phase="liq")
    mix.clear_runtime_caches()
    mix.reset_runtime_cache_stats()

    seeded = mix.state(T=temperature, x=composition, P=pressure, phase="liq", rho_guess=reference.density())

    assert seeded.density() == pytest.approx(reference.density())
    diagnostics = dict(mix._native.last_density_diagnostics())
    assert diagnostics["density_warm_start_source"] == "rho_guess"
    assert diagnostics["density_validity_gate"] == "passed"
    stats = mix.runtime_cache_stats()
    assert stats["density_warm_start_hits"] == 1
    assert stats["density_warm_start_rejections"] == 0
