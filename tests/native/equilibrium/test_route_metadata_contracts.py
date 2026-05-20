from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from epcsaft import _core
from tests.native.equilibrium.test_route_builders import (
    _ascani_electrolyte_mixture,
    _neutral_binary_mixture,
)

pytestmark = pytest.mark.native_contract

REPO_ROOT = Path(__file__).resolve().parents[3]


def test_neutral_two_phase_contract_declares_route_metadata_without_solving() -> None:
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

    payload = _core._native_neutral_two_phase_eos_nlp_contract(
        mix._native,
        temperature,
        target_pressure,
        [phase.tolist() for phase in phase_amounts],
        volumes,
        feed_amounts.tolist(),
    )

    assert payload["problem_name"] == "neutral_two_phase_eos"
    assert payload["variable_model"] == "phase_species_amounts_plus_phase_volume"
    assert payload["density_backend"] == "explicit_phase_volume_pressure_constraint"
    assert payload["residual_families"] == [
        "material_balance",
        "phase_pressure_consistency",
        "phase_equilibrium",
        "phase_distance",
    ]
    assert payload["constraint_families"] == ["material_balance", "phase_pressure_consistency"]


def test_seeded_lle_contract_declares_phase_distance_constraint_without_solving() -> None:
    mix = _neutral_binary_mixture()
    payload = _core._native_neutral_lle_eos_nlp_contract(
        mix._native,
        298.15,
        1.013e5,
        [0.45, 0.55],
    )

    assert payload["variable_model"] == "phase_species_amounts_plus_phase_volume"
    assert payload["constraint_families"] == [
        "material_balance",
        "phase_pressure_consistency",
        "phase_distance",
    ]


def test_fixed_composition_route_contracts_are_metadata_only_checks() -> None:
    mix = _neutral_binary_mixture()

    bubble = _core._native_neutral_bubble_p_eos_nlp_contract(mix._native, 300.0, [0.35, 0.65])
    dew = _core._native_neutral_dew_t_eos_nlp_contract(mix._native, 1.0e5, [0.35, 0.65])

    assert bubble["variable_model"] == "phase_species_amounts_plus_phase_volume_plus_pressure"
    assert bubble["constraint_families"] == [
        "fixed_composition",
        "phase_amount_total",
        "phase_pressure_consistency",
        "phase_equilibrium",
        "phase_volume_gap",
    ]
    assert dew["variable_model"] == "phase_species_amounts_plus_phase_volume_plus_temperature"
    assert dew["constraint_families"] == bubble["constraint_families"]


def test_liquid_root_contracts_declare_residual_families_without_running_ipopt() -> None:
    electrolyte_mix, electrolyte_feed = _ascani_electrolyte_mixture()
    electrolyte = _core._native_electrolyte_lle_eos_nlp_contract(
        electrolyte_mix._native,
        298.15,
        1.0e5,
        electrolyte_feed,
    )

    reactive_mix = _neutral_binary_mixture()
    reactive = _core._native_reactive_lle_eos_nlp_contract(
        reactive_mix._native,
        300.0,
        1.0e5,
        [0.3, 0.7],
        1,
        [1.0, 1.0],
        [1.0],
        1,
        [-1.0, 1.0],
        [float(np.log(3.0))],
    )

    assert electrolyte["variable_model"] == "ascani_transformed_salt_pairs_explicit_density"
    assert electrolyte["residual_families"] == ["phase_equilibrium", "material_balance"]
    assert reactive["variable_model"] == "log_phase_species_amounts_plus_log_density"
    assert reactive["residual_families"] == [
        "conserved_balance",
        "reaction_stationarity",
        "phase_equilibrium",
        "phase_charge",
    ]


def test_liquid_root_and_nlp_routes_share_metadata_type() -> None:
    header = (REPO_ROOT / "src" / "epcsaft" / "native" / "epcsaft_equilibrium.h").read_text(encoding="utf-8")

    assert "struct NativeRouteMetadata" not in header
    assert "using NativeRouteMetadata = epcsaft::native::equilibrium_nlp::RouteMetadata;" in header


def test_liquid_root_route_metadata_factories_live_with_nlp_metadata() -> None:
    metadata_header = (
        REPO_ROOT / "src" / "epcsaft" / "native" / "equilibrium_nlp" / "route_metadata.h"
    ).read_text(encoding="utf-8")
    equilibrium_header = (REPO_ROOT / "src" / "epcsaft" / "native" / "epcsaft_equilibrium.h").read_text(
        encoding="utf-8"
    )

    assert "inline RouteMetadata electrolyte_liquid_root_route_metadata(" in metadata_header
    assert "inline RouteMetadata reactive_liquid_root_route_metadata(" in metadata_header
    assert "NativeRouteMetadata electrolyte_liquid_root_route_metadata" not in equilibrium_header
    assert "NativeRouteMetadata reactive_liquid_root_route_metadata" not in equilibrium_header


def test_phase_tagged_reactive_liquid_root_contract_declares_reaction_constraint() -> None:
    mix = _neutral_binary_mixture()
    contract = _core._native_reactive_lle_eos_nlp_contract(
        mix._native,
        300.0,
        1.0e5,
        [0.3, 0.7],
        1,
        [1.0, 1.0],
        [1.0],
        1,
        [-1.0, 1.0],
        [float(np.log(3.0))],
        [0],
        [-1.0, 0.0, 0.0, 1.0],
    )

    assert contract["constraint_families"] == [
        "conserved_balance",
        "reaction_stationarity",
        "phase_pressure_consistency",
        "phase_distance",
    ]
