from __future__ import annotations

import math

import numpy as np
import pytest

import epcsaft


def _salt_speciation_mixture() -> epcsaft.ePCSAFTMixture:
    params = {
        "m": np.asarray([1.2047, 1.0, 1.0, 1.0]),
        "s": np.asarray([2.7927, 3.0, 2.8232, 2.7560]),
        "e": np.asarray([353.95, 200.0, 230.0, 170.0]),
        "z": np.asarray([0.0, 0.0, 1.0, -1.0]),
        "dielc": np.asarray([78.09, 8.0, 8.0, 8.0]),
        "d_born": np.asarray([0.0, 0.0, 3.445, 4.1]),
        "MW": np.asarray([18.01528e-3, 58.44e-3, 22.989e-3, 35.45e-3]),
    }
    return epcsaft.ePCSAFTMixture.from_params(params, species=["H2O", "NaCl", "Na+", "Cl-"])


def _concentration_reaction_case() -> tuple[epcsaft.ePCSAFTMixture, float, list[float]]:
    mix = _salt_speciation_mixture()
    initial_x = [0.998, 0.001, 0.0005, 0.0005]
    density = mix.state(T=298.15, P=1.0e5, x=initial_x, phase="liq").molar_density()
    log_k = math.log(density * initial_x[2]) + math.log(density * initial_x[3])
    log_k -= math.log(density * initial_x[1])
    return mix, log_k, initial_x


def _solve_concentration_case(
    *,
    mix: epcsaft.ePCSAFTMixture,
    log_k: float,
    initial_x: list[float],
    options: epcsaft.ReactiveSpeciationOptions,
) -> epcsaft.ReactiveSpeciationResult:
    species = ["H2O", "NaCl", "Na+", "Cl-"]
    return epcsaft.solve_reactive_speciation(
        species=species,
        mixture_factory=lambda x, T, P: mix,
        T=298.15,
        P=1.0e5,
        balances={
            "water_total": {"H2O": 1.0},
            "sodium_total": {"NaCl": 1.0, "Na+": 1.0},
            "chloride_total": {"NaCl": 1.0, "Cl-": 1.0},
        },
        totals={"water_total": 0.998, "sodium_total": 0.0015, "chloride_total": 0.0015},
        reactions=[
            epcsaft.ReactionDefinition(
                stoichiometry={"NaCl": -1.0, "Na+": 1.0, "Cl-": 1.0},
                log_equilibrium_constant=log_k,
                standard_state="concentration",
            )
        ],
        initial_x=initial_x,
        options=options,
    )


def test_auto_unsupported_derivative_fallback_is_rejected_without_debug_gate() -> None:
    mix, log_k, initial_x = _concentration_reaction_case()

    with pytest.raises(Exception, match="unsupported_derivative"):
        _solve_concentration_case(
            mix=mix,
            log_k=log_k,
            initial_x=initial_x,
            options=epcsaft.ReactiveSpeciationOptions(max_iterations=20),
        )


def test_explicit_unsupported_derivative_requires_debug_gate() -> None:
    mix, log_k, initial_x = _concentration_reaction_case()

    with pytest.raises(Exception, match="debug-only"):
        _solve_concentration_case(
            mix=mix,
            log_k=log_k,
            initial_x=initial_x,
            options=epcsaft.ReactiveSpeciationOptions(
                max_iterations=20,
                jacobian_backend="unsupported_derivative",
            ),
        )


def test_explicit_unsupported_derivative_runs_when_debug_gate_is_enabled(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("EPCSAFT_ALLOW_DERIVATIVE_BACKEND_DEBUG", "1")
    mix, log_k, initial_x = _concentration_reaction_case()

    result = _solve_concentration_case(
        mix=mix,
        log_k=log_k,
        initial_x=initial_x,
        options=epcsaft.ReactiveSpeciationOptions(
            max_iterations=20,
            jacobian_backend="unsupported_derivative",
        ),
    )

    assert result.success is True
    assert result.diagnostics["jacobian_backend"] == "unsupported_derivative"
    assert result.diagnostics["unsupported_derivative_allowed"] is True
    assert result.diagnostics["explicit_unsupported_derivative"] is True



