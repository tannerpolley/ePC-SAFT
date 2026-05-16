from __future__ import annotations

import math

import numpy as np
import pytest

import epcsaft
from epcsaft import _core


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

def _mea_like_mixture() -> epcsaft.ePCSAFTMixture:
    params = {
        "m": np.asarray([1.2047, 2.5, 1.0, 1.0, 2.2, 2.7, 1.0]),
        "s": np.asarray([2.7927, 3.3, 3.0, 3.2, 3.5, 3.6, 2.0]),
        "e": np.asarray([353.95, 260.0, 190.0, 230.0, 250.0, 245.0, 120.0]),
        "z": np.asarray([0.0, 0.0, 0.0, -1.0, 1.0, -1.0, 1.0]),
        "dielc": np.asarray([78.09, 35.0, 12.0, 20.0, 25.0, 22.0, 8.0]),
        "d_born": np.asarray([0.0, 0.0, 0.0, 3.8, 3.4, 4.0, 2.0]),
        "MW": np.asarray([18.01528e-3, 61.08e-3, 44.01e-3, 61.02e-3, 62.09e-3, 104.1e-3, 1.008e-3]),
    }
    return epcsaft.ePCSAFTMixture.from_params(
        params,
        species=["H2O", "MEA", "CO2", "HCO3-", "MEAH+", "MEACOO-", "H+"],
    )

def _methanol_cyclohexane_mixture(kij: float = 0.051) -> epcsaft.ePCSAFTMixture:
    params = {
        "MW": np.asarray([32.042e-3, 84.147e-3]),
        "m": np.asarray([1.5255, 2.5303]),
        "s": np.asarray([3.2300, 3.8499]),
        "e": np.asarray([188.90, 278.11]),
        "e_assoc": np.asarray([2899.5, 0.0]),
        "vol_a": np.asarray([0.035176, 0.0]),
        "assoc_scheme": ["2B", None],
        "k_ij": np.asarray([[0.0, kij], [kij, 0.0]]),
        "z": np.asarray([0.0, 0.0]),
        "dielc": np.asarray([33.05, 2.02]),
    }
    return epcsaft.ePCSAFTMixture.from_params(params, species=["Methanol", "Cyclohexane"])

def _log_k_from_state(
    mix: epcsaft.ePCSAFTMixture,
    T: float,
    P: float,
    x: np.ndarray,
    stoichiometry: dict[str, float],
) -> float:
    state = mix.state(T=T, P=P, x=x, phase="liq")
    gamma = state.activity_coefficient(species=mix.species)
    return float(
        sum(
            nu * math.log(max(x[mix.species.index(label)] * gamma[label], 1.0e-300))
            for label, nu in stoichiometry.items()
        )
    )

def _neutral_log_k_from_fugacity_activity(
    mix: epcsaft.ePCSAFTMixture,
    T: float,
    P: float,
    x: np.ndarray,
    stoichiometry: dict[str, float],
) -> float:
    state = mix.state(T=T, P=P, x=x, phase="liq")
    ln_phi = state.fugacity_coefficient(natural_log=True)
    ln_gamma = []
    for idx in range(mix.ncomp):
        x_ref = np.full(mix.ncomp, 1.0e-14, dtype=float)
        x_ref[idx] = 1.0 - 1.0e-14 * float(mix.ncomp - 1)
        ref = mix.state(T=T, P=P, x=x_ref, phase="liq")
        ln_phi_ref = ref.fugacity_coefficient(natural_log=True)
        ln_gamma.append(float(ln_phi[idx] - ln_phi_ref[idx]))
    return float(
        sum(
            nu * (math.log(max(x[mix.species.index(label)], 1.0e-300)) + ln_gamma[mix.species.index(label)])
            for label, nu in stoichiometry.items()
        )
    )

def test_native_chemical_equilibrium_residual_evaluator_rejects_removed_backend() -> None:
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )
    request = {
        "T": 298.15,
        "P": 1.0e5,
        "initial_x": [0.5, 0.5],
        "balance_matrix": [1.0, 1.0],
        "balance_rows": 1,
        "total_vector": [1.0],
        "reaction_stoichiometry": [-1.0, 1.0],
        "reaction_rows": 1,
        "log_equilibrium_constants": [math.log(3.0)],
        "reaction_standard_states": [1],
        "options": {"jacobian_backend": "finite" + "_difference", "finite" + "_difference_step": 1.0e-7},
    }

    with pytest.raises(_core.NativeValueError, match="jacobian_backend"):
        _core._evaluate_chemical_equilibrium_residual_native(mix._native, request)

def test_mixture_equilibrium_rejects_non_native_chemical_equilibrium_backend() -> None:
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )

    with pytest.raises(epcsaft.InputError, match="native"):
        mix.equilibrium(
            kind="chemical_equilibrium",
            T=298.15,
            P=1.0e5,
            z=[0.5, 0.5],
            balances={"total": {"A": 1.0, "B": 1.0}},
            totals={"total": 1.0},
            reactions=[epcsaft.ReactionDefinition({"A": -1.0, "B": 1.0}, math.log(3.0))],
            backend="legacy",
        )

def test_native_chemical_equilibrium_skips_soft_start_when_no_reactions() -> None:
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )

    result = epcsaft.solve_reactive_speciation(
        species=["A", "B"],
        mixture_factory=lambda x, T, P: mix,
        T=298.15,
        P=1.0e5,
        balances={"total": {"A": 1.0, "B": 1.0}},
        totals={"total": 1.0},
        reactions=[],
        initial_x=[0.2, 0.8],
        options=epcsaft.ReactiveSpeciationOptions(tolerance=1.0e-10),
    )

    assert result.success is True
    assert result.diagnostics["soft_start_enabled"] is True
    assert result.diagnostics["soft_start_attempted"] is False
    assert result.diagnostics["soft_start_success"] is False
    assert result.diagnostics["soft_start_used"] is False
    assert result.diagnostics["soft_start_rejection_reason"] == "no_reactions"
    assert result.diagnostics["implicit_solve_results"] == {}

def test_native_chemical_equilibrium_handles_trace_species_seed_without_invalid_numbers() -> None:
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )
    trace_target = 1.0e-9
    log_k = math.log(trace_target / (1.0 - trace_target))

    result = epcsaft.solve_reactive_speciation(
        species=["A", "B"],
        mixture_factory=lambda x, T, P: mix,
        T=298.15,
        P=1.0e5,
        balances={"total": {"A": 1.0, "B": 1.0}},
        totals={"total": 1.0},
        reactions=[
            epcsaft.ReactionDefinition(
                {"A": -1.0, "B": 1.0},
                log_k,
                standard_state="ideal_mole_fraction",
            )
        ],
        initial_x=[1.0 - trace_target, trace_target],
        options=epcsaft.ReactiveSpeciationOptions(
            tolerance=1.0e-8,
            min_mole_fraction=1.0e-14,
        ),
    )

    x_values = list(result.x.values())
    activity_values = list(result.activity_coefficients.values())
    diagnostics = result.diagnostics

    assert result.success is True
    assert min(x_values) >= 0.0
    assert result.x["B"] == pytest.approx(trace_target, rel=1.0e-6)
    assert max(abs(value) for value in result.reaction_residuals) <= 1.0e-8
    assert all(math.isfinite(value) for value in x_values)
    assert all(math.isfinite(value) and value > 0.0 for value in activity_values)
    assert all(math.isfinite(value) for value in result.mass_balance_residuals.values())
    assert math.isfinite(result.charge_residual)
    assert all(math.isfinite(value) for value in result.reaction_residuals)
    assert all(math.isfinite(value) for value in diagnostics["history"])
    assert all(math.isfinite(value) for value in diagnostics["phase_equilibrium_handoff"]["composition"])
