from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from tests.api.reactive.test_reactive_speciation_options import _assert_reactive_speciation_route_pending


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


def test_solve_reactive_speciation_auto_requires_native_ipopt_route_before_result_payload() -> None:
    species = ["H2O", "NaCl", "Na+", "Cl-"]
    mix = _salt_speciation_mixture()

    with pytest.raises(epcsaft.InputError) as excinfo:
        epcsaft.solve_reactive_speciation(
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
                    log_equilibrium_constant=100.0,
                    name="salt_dissociation",
                    standard_state="mole_fraction_activity",
                )
            ],
            initial_x=[0.998, 0.001, 0.0005, 0.0005],
            options=epcsaft.ReactiveSpeciationOptions(
                max_iterations=0,
                tolerance=1.0e-12,
            ),
        )

    _assert_reactive_speciation_route_pending(excinfo)


def test_solve_reactive_speciation_result_mode_does_not_mask_route_gate() -> None:
    species = ["H2O", "NaCl", "Na+", "Cl-"]
    mix = _salt_speciation_mixture()

    with pytest.raises(epcsaft.InputError) as excinfo:
        epcsaft.solve_reactive_speciation(
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
                    log_equilibrium_constant=100.0,
                    name="salt_dissociation",
                    standard_state="mole_fraction_activity",
                )
            ],
            initial_x=[0.998, 0.001, 0.0005, 0.0005],
            options=epcsaft.ReactiveSpeciationOptions(
                max_iterations=0,
                tolerance=1.0e-12,
                error_mode="result",
            ),
        )

    _assert_reactive_speciation_route_pending(excinfo)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        (
            {"reactions": [epcsaft.ReactionDefinition({"Missing": -1.0, "Na+": 1.0}, 0.0)]},
            "Unknown species 'Missing' in reaction stoichiometry",
        ),
        (
            {"totals": {"water_total": 0.998, "sodium_total": 0.0015}},
            "Missing total for balance 'chloride_total'",
        ),
        (
            {"initial_x": [0.998, 0.001, 0.0005]},
            "initial_x length must match species length",
        ),
    ],
)
def test_solve_reactive_speciation_rejects_invalid_chemistry_inputs(kwargs: dict[str, object], message: str) -> None:
    species = ["H2O", "NaCl", "Na+", "Cl-"]
    mix = _salt_speciation_mixture()
    request = {
        "species": species,
        "mixture_factory": lambda x, T, P: mix,
        "T": 298.15,
        "P": 1.0e5,
        "balances": {
            "water_total": {"H2O": 1.0},
            "sodium_total": {"NaCl": 1.0, "Na+": 1.0},
            "chloride_total": {"NaCl": 1.0, "Cl-": 1.0},
        },
        "totals": {"water_total": 0.998, "sodium_total": 0.0015, "chloride_total": 0.0015},
        "reactions": [
            epcsaft.ReactionDefinition(
                stoichiometry={"NaCl": -1.0, "Na+": 1.0, "Cl-": 1.0},
                log_equilibrium_constant=0.0,
            )
        ],
        "initial_x": [0.998, 0.001, 0.0005, 0.0005],
    }
    request.update(kwargs)

    with pytest.raises(epcsaft.InputError, match=message):
        epcsaft.solve_reactive_speciation(**request)
