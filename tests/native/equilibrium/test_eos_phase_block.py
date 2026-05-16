from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from epcsaft import _core


def _neutral_binary_mixture() -> epcsaft.ePCSAFTMixture:
    params = {
        "m": np.asarray([1.0, 1.6069]),
        "s": np.asarray([3.7039, 3.5206]),
        "e": np.asarray([150.03, 191.42]),
        "k_ij": np.asarray([[0.0, 3.0e-4], [3.0e-4, 0.0]]),
    }
    return epcsaft.ePCSAFTMixture.from_params(params, species=["Methane", "Ethane"])


def _two_phase_binary_case() -> tuple[
    epcsaft.ePCSAFTMixture,
    float,
    list[np.ndarray],
    list[float],
    np.ndarray,
    float,
]:
    mix = _neutral_binary_mixture()
    temperature = 300.0
    phase_amounts = [
        np.asarray([0.7, 0.3], dtype=float),
        np.asarray([0.1, 0.9], dtype=float),
    ]
    volumes = [float(phase_amounts[0].sum() / 80.0), float(phase_amounts[1].sum() / 140.0)]
    feed_amounts = phase_amounts[0] + phase_amounts[1]
    target_pressure = mix.state(
        T=temperature,
        rho=phase_amounts[0].sum() / volumes[0],
        x=phase_amounts[0] / phase_amounts[0].sum(),
        phase="liquid",
    ).pressure()
    return mix, temperature, phase_amounts, volumes, feed_amounts, target_pressure


def test_eos_phase_block_builds_amount_volume_variables_and_pressure_gate() -> None:
    mix = _neutral_binary_mixture()
    amounts = np.asarray([0.8, 1.2], dtype=float)
    temperature = 300.0
    density = 120.0
    volume = float(amounts.sum() / density)
    composition = amounts / amounts.sum()
    state = mix.state(T=temperature, rho=density, x=composition, phase="vapor")
    pressure = state.pressure()

    payload = _core._native_eos_phase_block(mix._native, temperature, pressure, amounts.tolist(), volume)

    assert payload["block"] == "eos_phase"
    assert payload["variable_names"] == ["n_0", "n_1", "volume"]
    assert payload["constraint_names"] == ["pressure_consistency"]
    assert payload["derivative_backend"] == "analytic"
    assert payload["composition"] == pytest.approx(composition.tolist(), abs=1.0e-14)
    assert payload["density"] == pytest.approx(density, rel=1.0e-14)
    assert payload["eos_pressure"] == pytest.approx(pressure, rel=1.0e-14)
    assert abs(payload["pressure_consistency_residual"]) <= 1.0e-8
    assert abs(payload["gradient"][-1]) <= 1.0e-10
    assert payload["objective_terms"]["pressure_work"] == pytest.approx(
        pressure * volume / payload["gas_constant_temperature"]
    )


def test_eos_phase_block_gradient_matches_chemical_potential_and_pressure_identities() -> None:
    mix = _neutral_binary_mixture()
    amounts = np.asarray([0.8, 1.2], dtype=float)
    temperature = 300.0
    density = 120.0
    volume = float(amounts.sum() / density)
    composition = amounts / amounts.sum()
    state = mix.state(T=temperature, rho=density, x=composition, phase="vapor")
    eos_pressure = state.pressure()
    target_pressure = eos_pressure + 2500.0

    payload = _core._native_eos_phase_block(mix._native, temperature, target_pressure, amounts.tolist(), volume)

    residual_mu = np.asarray(state.residual_chemical_potential(), dtype=float)
    expected_amount_gradient = np.log(amounts / volume) + residual_mu
    expected_volume_gradient = (target_pressure - eos_pressure) / payload["gas_constant_temperature"]

    assert payload["gradient"][:-1] == pytest.approx(expected_amount_gradient.tolist(), rel=1.0e-12, abs=1.0e-12)
    assert payload["gradient"][-1] == pytest.approx(expected_volume_gradient, rel=1.0e-12, abs=1.0e-12)
    assert payload["pressure_consistency_residual"] == pytest.approx(eos_pressure - target_pressure, abs=1.0e-8)
    assert np.isfinite(payload["objective"])


def test_eos_phase_block_pressure_jacobian_uses_exact_cppad_homogeneity_identity() -> None:
    mix = _neutral_binary_mixture()
    amounts = np.asarray([0.8, 1.2], dtype=float)
    temperature = 300.0
    density = 120.0
    volume = float(amounts.sum() / density)
    composition = amounts / amounts.sum()
    state = mix.state(T=temperature, rho=density, x=composition, phase="vapor")
    target_pressure = state.pressure() + 2500.0

    payload = _core._native_eos_phase_block(mix._native, temperature, target_pressure, amounts.tolist(), volume)

    jacobian = np.asarray(payload["pressure_jacobian"], dtype=float)
    assert payload["pressure_jacobian_backend"] == "cppad"
    assert payload["pressure_jacobian_shape"] == (1, 3)
    assert jacobian.size == 3
    assert jacobian[-1] == pytest.approx(
        -payload["pressure_density_derivative"] * payload["density"] / payload["volume"],
        rel=1.0e-12,
        abs=1.0e-8,
    )
    assert float(np.dot(amounts, jacobian[:-1]) + volume * jacobian[-1]) == pytest.approx(0.0, abs=1.0e-8)


def test_eos_phase_block_reports_pressure_constraint_jacobian_from_exact_curvature_identity() -> None:
    mix = _neutral_binary_mixture()
    amounts = np.asarray([0.8, 1.2], dtype=float)
    temperature = 300.0
    density = 120.0
    volume = float(amounts.sum() / density)
    composition = amounts / amounts.sum()
    state = mix.state(T=temperature, rho=density, x=composition, phase="vapor")

    payload = _core._native_eos_phase_block(mix._native, temperature, state.pressure(), amounts.tolist(), volume)

    assert payload["objective_curvature_backend"] == "cppad"
    assert payload["constraint_jacobian_backend"] == "cppad"
    assert payload["objective_curvature_shape"] == (3, 3)
    assert payload["constraint_jacobian_shape"] == (1, 3)
    objective_curvature = np.asarray(payload["objective_curvature_row_major"], dtype=float).reshape((3, 3))
    constraint_jacobian = np.asarray(payload["constraint_jacobian_row_major"], dtype=float).reshape((1, 3))
    expected_pressure_jacobian = -payload["gas_constant_temperature"] * objective_curvature[-1, :]

    assert np.all(np.isfinite(objective_curvature))
    assert np.all(np.isfinite(constraint_jacobian))
    assert constraint_jacobian[0, :] == pytest.approx(expected_pressure_jacobian, rel=1.0e-11, abs=1.0e-8)


def test_eos_phase_system_assembles_two_phase_material_balance_and_objective() -> None:
    mix, temperature, phase_amounts, volumes, feed_amounts, target_pressure = _two_phase_binary_case()

    payload = _core._native_eos_phase_system(
        mix._native,
        temperature,
        target_pressure,
        [phase.tolist() for phase in phase_amounts],
        volumes,
        feed_amounts.tolist(),
    )
    phase_blocks = [
        _core._native_eos_phase_block(mix._native, temperature, target_pressure, phase.tolist(), volume)
        for phase, volume in zip(phase_amounts, volumes, strict=True)
    ]

    assert payload["block"] == "eos_phase_system"
    assert payload["phase_count"] == 2
    assert payload["species_count"] == 2
    assert payload["variable_names"] == [
        "phase_0.n_0",
        "phase_0.n_1",
        "phase_0.volume",
        "phase_1.n_0",
        "phase_1.n_1",
        "phase_1.volume",
    ]
    assert payload["constraint_names"] == [
        "material_balance_0",
        "material_balance_1",
        "phase_0.pressure_consistency",
        "phase_1.pressure_consistency",
    ]
    assert payload["constraints"][:2] == pytest.approx([0.0, 0.0], abs=1.0e-14)
    assert payload["objective"] == pytest.approx(sum(block["objective"] for block in phase_blocks))
    assert payload["gradient"] == pytest.approx(
        phase_blocks[0]["gradient"] + phase_blocks[1]["gradient"],
        rel=1.0e-12,
        abs=1.0e-12,
    )


def test_eos_phase_system_reports_exact_material_and_pressure_jacobian_rows() -> None:
    mix, temperature, phase_amounts, volumes, feed_amounts, target_pressure = _two_phase_binary_case()

    payload = _core._native_eos_phase_system(
        mix._native,
        temperature,
        target_pressure,
        [phase.tolist() for phase in phase_amounts],
        volumes,
        feed_amounts.tolist(),
    )
    phase_blocks = [
        _core._native_eos_phase_block(mix._native, temperature, target_pressure, phase.tolist(), volume)
        for phase, volume in zip(phase_amounts, volumes, strict=True)
    ]

    assert payload["constraint_jacobian_backend"] == "analytic_cppad"
    assert payload["constraint_jacobian_shape"] == (4, 6)
    jacobian = np.asarray(payload["constraint_jacobian_row_major"], dtype=float).reshape((4, 6))
    assert np.array_equal(
        jacobian[:2, :],
        np.asarray(
            [
                [1.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0, 1.0, 0.0],
            ],
            dtype=float,
        ),
    )
    assert jacobian[2, :3] == pytest.approx(phase_blocks[0]["pressure_jacobian"], rel=1.0e-12, abs=1.0e-8)
    assert jacobian[2, 3:] == pytest.approx([0.0, 0.0, 0.0], abs=0.0)
    assert jacobian[3, :3] == pytest.approx([0.0, 0.0, 0.0], abs=0.0)
    assert jacobian[3, 3:] == pytest.approx(phase_blocks[1]["pressure_jacobian"], rel=1.0e-12, abs=1.0e-8)


def test_eos_phase_system_can_append_phase_charge_balance_rows() -> None:
    mix, temperature, phase_amounts, volumes, feed_amounts, target_pressure = _two_phase_binary_case()
    charges = np.asarray([1.0, -1.0], dtype=float)

    payload = _core._native_eos_phase_system(
        mix._native,
        temperature,
        target_pressure,
        [phase.tolist() for phase in phase_amounts],
        volumes,
        feed_amounts.tolist(),
        charges.tolist(),
    )

    assert payload["phase_charge_residuals"] == pytest.approx(
        [float(phase @ charges) for phase in phase_amounts],
        abs=1.0e-14,
    )
    assert payload["constraint_names"][-2:] == ["phase_0.charge_balance", "phase_1.charge_balance"]
    assert payload["constraint_jacobian_shape"] == (6, 6)
    jacobian = np.asarray(payload["constraint_jacobian_row_major"], dtype=float).reshape((6, 6))
    assert jacobian[4, :] == pytest.approx([1.0, -1.0, 0.0, 0.0, 0.0, 0.0], abs=0.0)
    assert jacobian[5, :] == pytest.approx([0.0, 0.0, 0.0, 1.0, -1.0, 0.0], abs=0.0)


def test_eos_phase_system_can_append_association_mass_action_rows() -> None:
    mix, temperature, phase_amounts, volumes, feed_amounts, target_pressure = _two_phase_binary_case()
    site_fractions = [
        np.asarray([0.8, 0.6], dtype=float),
        np.asarray([0.9, 0.7], dtype=float),
    ]
    delta = np.asarray([[0.0, 2.0], [3.0, 0.0]], dtype=float)

    payload = _core._native_eos_phase_system(
        mix._native,
        temperature,
        target_pressure,
        [phase.tolist() for phase in phase_amounts],
        volumes,
        feed_amounts.tolist(),
        [],
        [values.tolist() for values in site_fractions],
        delta.ravel().tolist(),
    )

    assert payload["variable_names"] == [
        "phase_0.n_0",
        "phase_0.n_1",
        "phase_0.volume",
        "phase_0.association_site_0",
        "phase_0.association_site_1",
        "phase_1.n_0",
        "phase_1.n_1",
        "phase_1.volume",
        "phase_1.association_site_0",
        "phase_1.association_site_1",
    ]
    assert payload["constraint_names"][-4:] == [
        "phase_0.association_site_0",
        "phase_0.association_site_1",
        "phase_1.association_site_0",
        "phase_1.association_site_1",
    ]
    assert payload["constraint_jacobian_shape"] == (8, 10)

    jacobian = np.asarray(payload["constraint_jacobian_row_major"], dtype=float).reshape((8, 10))
    association_rows = jacobian[4:, :]
    assert association_rows[:2, 5:] == pytest.approx(np.zeros((2, 5)), abs=0.0)
    assert association_rows[2:, :5] == pytest.approx(np.zeros((2, 5)), abs=0.0)

    for phase_index, (amounts, volume, fractions) in enumerate(
        zip(phase_amounts, volumes, site_fractions, strict=True)
    ):
        composition = amounts / amounts.sum()
        density = amounts.sum() / volume
        block = _core._native_association_mass_action_block(
            density,
            fractions.tolist(),
            composition.tolist(),
            delta.ravel().tolist(),
        )
        row_offset = phase_index * 2
        col_offset = phase_index * 5
        assert payload["phase_association_residuals"][row_offset : row_offset + 2] == pytest.approx(
            block["residuals"],
            rel=1.0e-14,
            abs=1.0e-14,
        )
        assert payload["constraints"][4 + row_offset : 6 + row_offset] == pytest.approx(
            block["residuals"],
            rel=1.0e-14,
            abs=1.0e-14,
        )

        expected_site_jacobian = np.asarray(block["site_fraction_jacobian_row_major"], dtype=float).reshape((2, 2))
        assert association_rows[row_offset : row_offset + 2, col_offset + 3 : col_offset + 5] == pytest.approx(
            expected_site_jacobian,
            rel=1.0e-14,
            abs=1.0e-14,
        )
        for site in range(2):
            expected_amount_derivative = fractions[site] * fractions * delta[site, :] / volume
            amount_sum = float(np.dot(amounts, fractions * delta[site, :]))
            expected_volume_derivative = -fractions[site] * amount_sum / (volume * volume)
            assert association_rows[row_offset + site, col_offset : col_offset + 2] == pytest.approx(
                expected_amount_derivative,
                rel=1.0e-14,
                abs=1.0e-14,
            )
            assert association_rows[row_offset + site, col_offset + 2] == pytest.approx(
                expected_volume_derivative,
                rel=1.0e-14,
                abs=1.0e-14,
            )
