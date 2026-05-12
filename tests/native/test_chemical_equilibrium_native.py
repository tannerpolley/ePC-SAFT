from __future__ import annotations

import math

import numpy as np
import pytest

import epcsaft
from epcsaft import _core


@pytest.fixture(autouse=True)
def _allow_finite_difference_debug(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("EPCSAFT_ALLOW_FINITE_DIFFERENCE_DEBUG", "1")


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


def test_native_chemical_equilibrium_entrypoint_is_exposed() -> None:
    assert hasattr(_core, "_solve_chemical_equilibrium_native")
    assert hasattr(_core, "_evaluate_chemical_equilibrium_residual_native")


def test_native_chemical_equilibrium_residual_evaluator_uses_analytic_jacobian_by_default() -> None:
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
        "options": {"tolerance": 1.0e-10},
    }

    payload = _core._evaluate_chemical_equilibrium_residual_native(mix._native, request)

    assert payload["variable_model"] == "log_species_amounts"
    assert payload["jacobian_backend"] == "analytic"
    assert payload["hessian_backend"] == "gauss_newton"
    diagnostics = payload["diagnostics"]
    assert diagnostics["derivative_backend_selected"] == "analytic"
    assert diagnostics["derivative_capability_path"] == "chemical_equilibrium:ideal_mole_fraction:log_amounts"
    assert diagnostics["finite_difference_allowed"] is False
    assert diagnostics["unsupported_derivative_reason"] == ""
    assert diagnostics["exact_hessian_available"] is False
    assert diagnostics["hessian_kind"] == "approximate_least_squares_gauss_newton"
    assert diagnostics["hessian_includes_second_residual_derivatives"] is False
    residual = np.asarray(payload["residual"], dtype=float)
    gradient = np.asarray(payload["gradient"], dtype=float)
    jacobian = np.asarray(payload["jacobian_row_major"], dtype=float).reshape(payload["jacobian_shape"])
    lower = np.asarray(payload["lower_bounds"], dtype=float)
    variables = np.asarray(payload["variables"], dtype=float)
    upper = np.asarray(payload["upper_bounds"], dtype=float)
    assert residual.shape == (3,)
    assert gradient.shape == (2,)
    assert jacobian.shape == (3, 2)
    assert np.isfinite(payload["objective"])
    assert np.all(np.isfinite(residual))
    assert np.all(np.isfinite(gradient))
    assert np.all(np.isfinite(jacobian))


def test_native_chemical_equilibrium_residual_evaluator_supports_explicit_autodiff() -> None:
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
        "options": {"tolerance": 1.0e-10, "jacobian_backend": "autodiff"},
    }

    cppad_enabled = bool(epcsaft.runtime_build_info()["native_dependencies"]["cppad"]["enabled"])
    if not cppad_enabled:
        with pytest.raises(ValueError, match="CppAD-enabled build"):
            _core._evaluate_chemical_equilibrium_residual_native(mix._native, request)
        return

    payload = _core._evaluate_chemical_equilibrium_residual_native(mix._native, request)

    assert payload["jacobian_backend"] == "autodiff"
    diagnostics = payload["diagnostics"]
    assert diagnostics["derivative_backend_selected"] == "autodiff"
    assert diagnostics["derivative_capability_path"] == "chemical_equilibrium:ideal_mole_fraction:log_amounts:cppad"
    jacobian = np.asarray(payload["jacobian_row_major"], dtype=float).reshape(payload["jacobian_shape"])
    residual = np.asarray(payload["residual"], dtype=float)
    gradient = np.asarray(payload["gradient"], dtype=float)
    lower = np.asarray(payload["lower_bounds"], dtype=float)
    variables = np.asarray(payload["variables"], dtype=float)
    upper = np.asarray(payload["upper_bounds"], dtype=float)
    assert np.all(np.isfinite(jacobian))
    assert np.all(np.isfinite(lower))
    assert np.all(np.isfinite(variables))
    assert np.all(np.isfinite(upper))
    assert np.all(lower < upper)
    assert np.all(variables >= lower)
    assert np.all(variables <= upper)
    np.testing.assert_allclose(gradient, jacobian.T @ residual, rtol=1.0e-10, atol=1.0e-10)
    assert payload["objective"] == pytest.approx(0.5 * float(residual @ residual))
    assert len(payload["lower_bounds"]) == len(payload["variables"]) == len(payload["upper_bounds"])


def test_native_chemical_equilibrium_residual_evaluator_keeps_explicit_finite_difference() -> None:
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
        "options": {"jacobian_backend": "finite_difference", "finite_difference_step": 1.0e-7},
    }

    payload = _core._evaluate_chemical_equilibrium_residual_native(mix._native, request)

    assert payload["jacobian_backend"] == "finite_difference"
    diagnostics = payload["diagnostics"]
    assert diagnostics["derivative_backend_selected"] == "finite_difference"
    assert diagnostics["finite_difference_allowed"] is True
    assert diagnostics["explicit_finite_difference"] is True
    assert diagnostics["finite_difference_scheme"] == "forward"


def test_native_chemical_equilibrium_concentration_autodiff_matches_finite_difference() -> None:
    mix = _salt_speciation_mixture()
    initial_x = np.asarray([0.998, 0.001, 0.0005, 0.0005], dtype=float)
    density = mix.state(T=298.15, P=1.0e5, x=initial_x, phase="liq").molar_density()
    log_k = math.log(density * initial_x[2]) + math.log(density * initial_x[3]) - math.log(density * initial_x[1])

    request = {
        "T": 298.15,
        "P": 1.0e5,
        "initial_x": initial_x.tolist(),
        "balance_matrix": [
            1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 1.0, 0.0,
            0.0, 1.0, 0.0, 1.0,
        ],
        "balance_rows": 3,
        "total_vector": [0.998, 0.0015, 0.0015],
        "reaction_stoichiometry": [0.0, -1.0, 1.0, 1.0],
        "reaction_rows": 1,
        "log_equilibrium_constants": [log_k],
        "reaction_standard_states": [2],
        "options": {"tolerance": 1.0e-9, "jacobian_backend": "autodiff"},
    }

    cppad_enabled = bool(epcsaft.runtime_build_info()["native_dependencies"]["cppad"]["enabled"])
    if not cppad_enabled:
        with pytest.raises(ValueError, match="CppAD-enabled build"):
            _core._evaluate_chemical_equilibrium_residual_native(mix._native, request)
        return

    autodiff_payload = _core._evaluate_chemical_equilibrium_residual_native(mix._native, request)
    request["options"] = {"tolerance": 1.0e-9, "jacobian_backend": "finite_difference", "finite_difference_step": 1.0e-7}
    fd_payload = _core._evaluate_chemical_equilibrium_residual_native(mix._native, request)

    assert autodiff_payload["jacobian_backend"] == "autodiff"
    assert autodiff_payload["diagnostics"]["derivative_capability_path"] == (
        "chemical_equilibrium:concentration:log_amounts:pressure_closure_cppad"
    )
    autodiff_jac = np.asarray(autodiff_payload["jacobian_row_major"], dtype=float).reshape(autodiff_payload["jacobian_shape"])
    fd_jac = np.asarray(fd_payload["jacobian_row_major"], dtype=float).reshape(fd_payload["jacobian_shape"])
    np.testing.assert_allclose(autodiff_jac, fd_jac, rtol=3.0e-4, atol=3.0e-4)


def test_native_chemical_equilibrium_activity_autodiff_matches_finite_difference() -> None:
    mix = _salt_speciation_mixture()
    initial_x = np.asarray([0.998, 0.001, 0.0005, 0.0005], dtype=float)
    stoich = {"NaCl": -1.0, "Na+": 1.0, "Cl-": 1.0}
    log_k = _log_k_from_state(mix, 298.15, 1.0e5, initial_x, stoich)

    request = {
        "T": 298.15,
        "P": 1.0e5,
        "initial_x": initial_x.tolist(),
        "balance_matrix": [
            1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 1.0, 0.0,
            0.0, 1.0, 0.0, 1.0,
        ],
        "balance_rows": 3,
        "total_vector": [initial_x[0], initial_x[1] + initial_x[2], initial_x[1] + initial_x[3]],
        "reaction_stoichiometry": [
            0.0, -1.0, 1.0, 1.0,
        ],
        "reaction_rows": 1,
        "log_equilibrium_constants": [log_k],
        "reaction_standard_states": [0],
        "options": {"tolerance": 1.0e-9, "jacobian_backend": "autodiff"},
    }

    cppad_enabled = bool(epcsaft.runtime_build_info()["native_dependencies"]["cppad"]["enabled"])
    if not cppad_enabled:
        with pytest.raises(ValueError, match="CppAD-enabled build"):
            _core._evaluate_chemical_equilibrium_residual_native(mix._native, request)
        return

    autodiff_payload = _core._evaluate_chemical_equilibrium_residual_native(mix._native, request)
    request["options"] = {"tolerance": 1.0e-9, "jacobian_backend": "finite_difference", "finite_difference_step": 1.0e-7}
    fd_payload = _core._evaluate_chemical_equilibrium_residual_native(mix._native, request)

    assert autodiff_payload["jacobian_backend"] == "autodiff"
    assert autodiff_payload["diagnostics"]["derivative_capability_path"] == (
        "chemical_equilibrium:mole_fraction_activity:log_amounts:component_activity_cppad"
    )
    autodiff_jac = np.asarray(autodiff_payload["jacobian_row_major"], dtype=float).reshape(autodiff_payload["jacobian_shape"])
    fd_jac = np.asarray(fd_payload["jacobian_row_major"], dtype=float).reshape(fd_payload["jacobian_shape"])
    np.testing.assert_allclose(autodiff_jac, fd_jac, rtol=5.0e-4, atol=5.0e-4)


def test_native_chemical_equilibrium_solves_easy_ideal_reaction() -> None:
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
        reactions=[
            epcsaft.ReactionDefinition(
                {"A": -1.0, "B": 1.0},
                log_equilibrium_constant=math.log(3.0),
                standard_state="ideal_mole_fraction",
            )
        ],
        initial_x=[0.5, 0.5],
        options=epcsaft.ReactiveSpeciationOptions(tolerance=1.0e-10),
    )

    assert result.success is True
    assert result.x["B"] / result.x["A"] == pytest.approx(3.0, rel=1.0e-8)
    assert result.diagnostics["solver_language"] == "c++"
    assert result.diagnostics["native_entrypoint"] == "_solve_chemical_equilibrium_native"
    assert result.diagnostics["activity_model"] == "epcsaft_neutral_fugacity_activity"
    assert result.diagnostics["requested_jacobian_backend"] == "auto"
    assert result.diagnostics["jacobian_backend"] == "analytic"
    assert result.diagnostics["jacobian_available"] is True
    assert result.diagnostics["jacobian_fallback_used"] is False
    assert result.diagnostics["finite_difference_fallback_used"] is False
    assert result.diagnostics["derivative_backend_selected"] == "analytic"
    assert result.diagnostics["finite_difference_allowed"] is False
    assert result.diagnostics["hessian_available"] is False
    assert result.diagnostics["hessian_backend"] == "not_implemented"
    assert result.diagnostics["hessian_fallback_used"] is False
    assert "IPOPT-compatible optimizer integration" in result.diagnostics["hessian_fallback_reason"]


def test_mixture_equilibrium_routes_chemical_equilibrium_to_native_speciation() -> None:
    mix = epcsaft.ePCSAFTMixture.from_params(
        {
            "m": np.asarray([1.0, 1.0]),
            "s": np.asarray([3.0, 3.0]),
            "e": np.asarray([200.0, 200.0]),
        },
        species=["A", "B"],
    )

    result = mix.equilibrium(
        kind="chemical_equilibrium",
        T=298.15,
        P=1.0e5,
        z=[0.5, 0.5],
        balances={"total": {"A": 1.0, "B": 1.0}},
        totals={"total": 1.0},
        reactions=[
            epcsaft.ReactionDefinition(
                {"A": -1.0, "B": 1.0},
                math.log(3.0),
                standard_state="ideal_mole_fraction",
            )
        ],
        options=epcsaft.ReactiveSpeciationOptions(tolerance=1.0e-10),
    )

    assert isinstance(result, epcsaft.ReactiveSpeciationResult)
    assert result.x["B"] / result.x["A"] == pytest.approx(3.0, rel=1.0e-8)
    assert result.diagnostics["native_entrypoint"] == "_solve_chemical_equilibrium_native"


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


def test_native_chemical_equilibrium_matches_activity_coupled_salt_speciation() -> None:
    mix = _salt_speciation_mixture()
    initial_x = np.asarray([0.998, 0.001, 0.0005, 0.0005], dtype=float)
    stoich = {"NaCl": -1.0, "Na+": 1.0, "Cl-": 1.0}
    log_k = _log_k_from_state(mix, 298.15, 1.0e5, initial_x, stoich)

    result = epcsaft.solve_reactive_speciation(
        species=mix.species,
        mixture_factory=lambda x, T, P: mix,
        T=298.15,
        P=1.0e5,
        balances={
            "water_total": {"H2O": 1.0},
            "sodium_total": {"NaCl": 1.0, "Na+": 1.0},
            "chloride_total": {"NaCl": 1.0, "Cl-": 1.0},
        },
        totals={"water_total": 0.998, "sodium_total": 0.0015, "chloride_total": 0.0015},
        reactions=[epcsaft.ReactionDefinition(stoich, log_k)],
        initial_x=initial_x,
        options=epcsaft.ReactiveSpeciationOptions(tolerance=1.0e-8, jacobian_backend="finite_difference"),
    )

    assert result.success is True
    assert max(abs(value) for value in result.mass_balance_residuals.values()) <= 1.0e-8
    assert abs(result.charge_residual) <= 1.0e-10
    assert max(abs(value) for value in result.reaction_residuals) <= 1.0e-8
    assert result.diagnostics["activity_model"] == "epcsaft_component_activity"


def test_native_chemical_equilibrium_solves_hard_mea_like_speciation_and_returns_phase_handoff() -> None:
    mix = _mea_like_mixture()
    target_x = np.asarray([0.865, 0.075, 0.02, 0.008, 0.012, 0.012, 0.008], dtype=float)
    initial_x = np.asarray([0.84, 0.09, 0.03, 0.004, 0.02, 0.006, 0.01], dtype=float)
    reactions = [
        {"CO2": -1.0, "H2O": -1.0, "HCO3-": 1.0, "H+": 1.0},
        {"MEA": -1.0, "H+": -1.0, "MEAH+": 1.0},
        {"MEA": -1.0, "CO2": -1.0, "MEACOO-": 1.0, "H+": 1.0},
    ]
    log_k = [_log_k_from_state(mix, 313.15, 1.0e5, target_x, reaction) for reaction in reactions]

    result = epcsaft.solve_reactive_speciation(
        species=mix.species,
        mixture_factory=lambda x, T, P: mix,
        T=313.15,
        P=1.0e5,
        balances={
            "water_total": {"H2O": 1.0},
            "amine_total": {"MEA": 1.0, "MEAH+": 1.0, "MEACOO-": 1.0},
            "carbon_total": {"CO2": 1.0, "HCO3-": 1.0, "MEACOO-": 1.0},
        },
        totals={
            "water_total": float(target_x[0]),
            "amine_total": float(target_x[1] + target_x[4] + target_x[5]),
            "carbon_total": float(target_x[2] + target_x[3] + target_x[5]),
        },
        reactions=[epcsaft.ReactionDefinition(reaction, value) for reaction, value in zip(reactions, log_k)],
        initial_x=initial_x,
        options=epcsaft.ReactiveSpeciationOptions(
            tolerance=1.0e-7,
            max_iterations=40,
            damping=0.7,
            jacobian_backend="finite_difference",
        ),
    )

    assert result.success is True
    assert max(abs(value) for value in result.mass_balance_residuals.values()) <= 1.0e-7
    assert abs(result.charge_residual) <= 1.0e-8
    assert max(abs(value) for value in result.reaction_residuals) <= 1.0e-7
    assert result.diagnostics["problem_class"] == "homogeneous_chemical_equilibrium"
    assert result.diagnostics["phase_equilibrium_handoff"]["composition"] == pytest.approx(list(result.x.values()))
    assert result.diagnostics["phase_equilibrium_handoff"]["activity_basis"] == "mole_fraction"


def test_native_chemical_equilibrium_uses_convex_soft_start_for_bad_neutral_seed() -> None:
    mix = _methanol_cyclohexane_mixture()
    target_x = np.asarray([0.35, 0.65], dtype=float)
    stoich = {"Methanol": -1.0, "Cyclohexane": 1.0}
    log_k = _neutral_log_k_from_fugacity_activity(mix, 298.15, 1.013e5, target_x, stoich)

    result = epcsaft.solve_reactive_speciation(
        species=mix.species,
        mixture_factory=lambda x, T, P: mix,
        T=298.15,
        P=1.013e5,
        balances={"total": {"Methanol": 1.0, "Cyclohexane": 1.0}},
        totals={"total": 1.0},
        reactions=[epcsaft.ReactionDefinition(stoich, log_k)],
        initial_x=[0.999, 0.001],
        options=epcsaft.ReactiveSpeciationOptions(
            tolerance=1.0e-9,
            max_iterations=35,
            jacobian_backend="finite_difference",
        ),
    )

    assert result.success is True
    assert result.x["Methanol"] == pytest.approx(target_x[0], rel=5.0e-6)
    assert result.x["Cyclohexane"] == pytest.approx(target_x[1], rel=5.0e-6)
    assert result.diagnostics["soft_start_enabled"] is True
    assert result.diagnostics["soft_start_attempted"] is True
    assert result.diagnostics["soft_start_success"] is True
    assert result.diagnostics["soft_start_used"] is True
    assert result.diagnostics["soft_start_rejection_reason"] == ""
    assert result.diagnostics["soft_start_iterations"] >= 0
    assert len(result.diagnostics["soft_start_composition"]) == len(mix.species)
    assert min(result.diagnostics["soft_start_composition"]) > 0.0
    assert result.diagnostics["soft_start_composition"] != pytest.approx([0.999, 0.001])
    assert result.diagnostics["soft_start_initial_residual_norm"] > result.diagnostics["soft_start_residual_norm"]
    assert len(result.diagnostics["soft_start_history"]) >= 1


def test_native_chemical_equilibrium_soft_start_reports_for_hard_ionic_speciation() -> None:
    mix = _mea_like_mixture()
    target_x = np.asarray([0.865, 0.075, 0.02, 0.008, 0.012, 0.012, 0.008], dtype=float)
    hard_initial_x = np.asarray([0.93, 0.035, 0.026, 0.001, 0.002, 0.001, 0.005], dtype=float)
    reactions = [
        {"CO2": -1.0, "H2O": -1.0, "HCO3-": 1.0, "H+": 1.0},
        {"MEA": -1.0, "H+": -1.0, "MEAH+": 1.0},
        {"MEA": -1.0, "CO2": -1.0, "MEACOO-": 1.0, "H+": 1.0},
    ]
    log_k = [_log_k_from_state(mix, 313.15, 1.0e5, target_x, reaction) for reaction in reactions]

    result = epcsaft.solve_reactive_speciation(
        species=mix.species,
        mixture_factory=lambda x, T, P: mix,
        T=313.15,
        P=1.0e5,
        balances={
            "water_total": {"H2O": 1.0},
            "amine_total": {"MEA": 1.0, "MEAH+": 1.0, "MEACOO-": 1.0},
            "carbon_total": {"CO2": 1.0, "HCO3-": 1.0, "MEACOO-": 1.0},
        },
        totals={
            "water_total": float(target_x[0]),
            "amine_total": float(target_x[1] + target_x[4] + target_x[5]),
            "carbon_total": float(target_x[2] + target_x[3] + target_x[5]),
        },
        reactions=[epcsaft.ReactionDefinition(reaction, value) for reaction, value in zip(reactions, log_k)],
        initial_x=hard_initial_x,
        options=epcsaft.ReactiveSpeciationOptions(
            tolerance=1.0e-7,
            max_iterations=45,
            damping=0.7,
            jacobian_backend="finite_difference",
        ),
    )

    assert result.success is True
    assert result.diagnostics["activity_model"] == "epcsaft_component_activity"
    assert result.diagnostics["soft_start_enabled"] is True
    assert result.diagnostics["soft_start_attempted"] is True
    assert result.diagnostics["soft_start_success"] is True
    assert len(result.diagnostics["soft_start_composition"]) == len(mix.species)
    assert min(result.diagnostics["soft_start_composition"]) > 0.0
    assert len(result.diagnostics["soft_start_history"]) >= 1


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
        initial_x=[1.0 - 1.0e-14, 1.0e-14],
        options=epcsaft.ReactiveSpeciationOptions(
            tolerance=1.0e-8,
            min_mole_fraction=1.0e-14,
            jacobian_backend="finite_difference",
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


def test_reactive_stability_chemical_equilibrates_feed_before_native_tpd() -> None:
    mix = _methanol_cyclohexane_mixture()
    target_x = np.asarray([0.45, 0.55], dtype=float)
    stoich = {"Methanol": -1.0, "Cyclohexane": 1.0}
    log_k = _neutral_log_k_from_fugacity_activity(mix, 298.15, 1.013e5, target_x, stoich)

    result = mix.equilibrium(
        kind="reactive_stability",
        T=298.15,
        P=1.013e5,
        z=[0.5, 0.5],
        balances={"total": {"Methanol": 1.0, "Cyclohexane": 1.0}},
        totals={"total": 1.0},
        reactions=[
            epcsaft.ReactionDefinition(
                stoich,
                log_k,
            )
        ],
        parent_phase="liq",
        trial_phases=("liq",),
        options=epcsaft.ReactiveSpeciationOptions(tolerance=1.0e-10, jacobian_backend="finite_difference"),
    )

    assert isinstance(result, epcsaft.StabilityResult)
    assert result.backend == "neutral_tpd"
    assert result.stable is False
    assert result.diagnostics["reactive_phase_method"] == "chemical_equilibrium_then_native_stability"
    assert result.diagnostics["reactive_chemical_equilibrium"]["x"]["Methanol"] == pytest.approx(0.45, rel=1.0e-8)
    assert result.diagnostics["reactive_chemical_equilibrium"]["x"]["Cyclohexane"] == pytest.approx(0.55, rel=1.0e-8)
    assert result.diagnostics["reactive_feed_composition"] == pytest.approx([0.45, 0.55])
    assert result.diagnostics["reactive_chemical_equilibrium"]["diagnostics"]["soft_start_attempted"] is True
    assert result.diagnostics["reactive_chemical_equilibrium"]["diagnostics"]["soft_start_used"] is True


def test_native_chemical_equilibrium_uses_epcsaft_activities_for_neutral_reaction() -> None:
    mix = _methanol_cyclohexane_mixture()
    target_x = np.asarray([0.35, 0.65], dtype=float)
    stoich = {"Methanol": -1.0, "Cyclohexane": 1.0}
    log_k = _neutral_log_k_from_fugacity_activity(mix, 298.15, 1.013e5, target_x, stoich)

    result = epcsaft.solve_reactive_speciation(
        species=mix.species,
        mixture_factory=lambda x, T, P: mix,
        T=298.15,
        P=1.013e5,
        balances={"total": {"Methanol": 1.0, "Cyclohexane": 1.0}},
        totals={"total": 1.0},
        reactions=[epcsaft.ReactionDefinition(stoich, log_k)],
        initial_x=[0.5, 0.5],
        options=epcsaft.ReactiveSpeciationOptions(tolerance=1.0e-9, jacobian_backend="finite_difference"),
    )

    assert result.success is True
    assert result.x["Methanol"] == pytest.approx(target_x[0], rel=5.0e-6)
    assert result.x["Cyclohexane"] == pytest.approx(target_x[1], rel=5.0e-6)
    assert result.diagnostics["activity_model"] == "epcsaft_neutral_fugacity_activity"
    assert result.activity_coefficients["Methanol"] != pytest.approx(1.0, abs=1.0e-3)


def test_native_chemical_equilibrium_solution_shifts_when_fugacity_model_changes() -> None:
    base_mix = _methanol_cyclohexane_mixture(kij=0.051)
    target_x = np.asarray([0.35, 0.65], dtype=float)
    stoich = {"Methanol": -1.0, "Cyclohexane": 1.0}
    log_k = _neutral_log_k_from_fugacity_activity(base_mix, 298.15, 1.013e5, target_x, stoich)

    base_result = epcsaft.solve_reactive_speciation(
        species=base_mix.species,
        mixture_factory=lambda x, T, P: base_mix,
        T=298.15,
        P=1.013e5,
        balances={"total": {"Methanol": 1.0, "Cyclohexane": 1.0}},
        totals={"total": 1.0},
        reactions=[epcsaft.ReactionDefinition(stoich, log_k)],
        initial_x=[0.5, 0.5],
        options=epcsaft.ReactiveSpeciationOptions(tolerance=1.0e-9, jacobian_backend="finite_difference"),
    )

    perturbed_mix = _methanol_cyclohexane_mixture(kij=0.0)
    perturbed_result = epcsaft.solve_reactive_speciation(
        species=perturbed_mix.species,
        mixture_factory=lambda x, T, P: perturbed_mix,
        T=298.15,
        P=1.013e5,
        balances={"total": {"Methanol": 1.0, "Cyclohexane": 1.0}},
        totals={"total": 1.0},
        reactions=[epcsaft.ReactionDefinition(stoich, log_k)],
        initial_x=[0.5, 0.5],
        options=epcsaft.ReactiveSpeciationOptions(tolerance=1.0e-9, jacobian_backend="finite_difference"),
    )

    assert base_result.success is True
    assert base_result.x["Methanol"] == pytest.approx(target_x[0], rel=5.0e-6)
    assert base_result.x["Cyclohexane"] == pytest.approx(target_x[1], rel=5.0e-6)
    assert base_result.diagnostics["activity_model"] == "epcsaft_neutral_fugacity_activity"

    assert perturbed_result.success is True
    assert perturbed_result.diagnostics["activity_model"] == "epcsaft_neutral_fugacity_activity"
    assert max(abs(value) for value in perturbed_result.reaction_residuals) <= 1.0e-8
    assert perturbed_result.x["Methanol"] - base_result.x["Methanol"] > 0.1
