from __future__ import annotations

import json

import numpy as np
import pytest

import epcsaft
from epcsaft import ePCSAFTMixture


def test_parameter_set_round_trips_to_legacy_payload_and_mixture() -> None:
    params = epcsaft.ParameterSet.from_records(
        [
            epcsaft.PureRecord(
                "Methanol",
                molar_mass=32.04,
                m=1.5255,
                sigma=3.23,
                epsilon_k=188.9,
                epsilon_k_ab=2899.5,
                kappa_ab=0.035176,
                association_scheme="2B",
                association_sites=(epcsaft.AssociationSite("H"), epcsaft.AssociationSite("O")),
            ),
            epcsaft.PureRecord("Cyclohexane", molar_mass=84.16, m=2.5303, sigma=3.8499, epsilon_k=278.11),
        ],
        [epcsaft.BinaryRecord(("Methanol", "Cyclohexane"), k_ij=0.051)],
        metadata={"source": "unit-test"},
    )

    legacy = params.to_legacy_dict()
    assert legacy["m"].tolist() == pytest.approx([1.5255, 2.5303])
    assert legacy["k_ij"][0, 1] == pytest.approx(0.051)
    assert legacy["assoc_scheme"] == ["2B", None]

    mix = ePCSAFTMixture.from_params(params)
    state = mix.state(T=298.15, P=101325.0, x=[0.5, 0.5])

    assert mix.species == ["Methanol", "Cyclohexane"]
    assert state.pressure() == pytest.approx(101325.0, rel=1.0e-5)


def test_parameter_set_from_legacy_dict_preserves_binary_records() -> None:
    payload = {
        "m": np.asarray([1.0, 2.0]),
        "s": np.asarray([3.0, 4.0]),
        "e": np.asarray([200.0, 250.0]),
        "MW": np.asarray([18.0, 46.0]),
        "k_ij": np.asarray([[0.0, 0.12], [0.12, 0.0]]),
    }

    params = epcsaft.ParameterSet.from_dict(payload, species=["water", "ethanol"])

    assert params.components == ("water", "ethanol")
    assert params.binary_records == (epcsaft.BinaryRecord(("water", "ethanol"), k_ij=0.12),)
    json.loads(params.to_json())


def test_parameter_set_validation_rejects_missing_records() -> None:
    with pytest.raises(epcsaft.InputError, match="Missing pure records"):
        epcsaft.ParameterSet(
            components=("A", "B"),
            pure_records=(epcsaft.PureRecord("A", molar_mass=10.0, m=1.0, sigma=3.0, epsilon_k=200.0),),
        )
