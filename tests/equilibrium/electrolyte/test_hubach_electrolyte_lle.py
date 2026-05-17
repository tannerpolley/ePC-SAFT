from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

import epcsaft
from epcsaft import ePCSAFTMixture
from tests.helpers.numeric import assert_allclose

REPO_ROOT = Path(__file__).resolve().parents[3]
FIXTURE_DIR = REPO_ROOT / "data" / "reference" / "equilibrium_benchmarks" / "electrolyte_lle" / "hubach_2024"
SPECIES = ["H2O", "TBP", "[emim][tcb]", "Li+", "Cl-"]
T_K = 294.15
P_PA = 1.013e5


def _row0_feed() -> np.ndarray:
    with (FIXTURE_DIR / "feed_compositions.csv").open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle))
    return np.asarray([float(row[name]) for name in SPECIES], dtype=float)


def _hubach_mixture(feed: np.ndarray) -> ePCSAFTMixture:
    return ePCSAFTMixture.from_dataset("2024_Hubach", SPECIES, feed, T_K)


def test_hubach_fixture_loads() -> None:
    feed = _row0_feed()
    mix = _hubach_mixture(feed)

    assert list(mix.species) == SPECIES
    assert_allclose(mix.parameters["z"], [0.0, 0.0, 0.0, 1.0, -1.0])
    assert float(np.sum(feed)) == pytest.approx(1.0)
    assert abs(float(np.dot(feed, mix.parameters["z"]))) <= 1.0e-12


def test_hubach_fixture_matches_lithium_canonical_option_surface() -> None:
    feed = _row0_feed()
    mix = _hubach_mixture(feed)
    model = mix.parameters["elec_model"]

    assert model["rel_perm"] == {"rule": 3, "differential_mode": 0}
    assert model["hc_model"] == {"dadx_differential_mode": 0}
    assert model["disp_model"] == {"dadx_differential_mode": 0}
    assert model["assoc_model"] == {"dadx_differential_mode": 0}
    assert model["DH_model"]["bjeruum_treatment"] is False
    assert model["DH_model"]["mu_DH_model"]["differential_mode"] == 0
    assert model["include_born_model"] is True
    assert model["born_model"]["d_Born_mode"] == 3
    assert model["born_model"]["solvation_shell_model"] is True
    assert model["born_model"]["dielectric_saturation"] is True
    assert model["born_model"]["bulk_mode"] == "mix"
    assert model["born_model"]["mu_born_model"]["differential_mode"] == 0
    assert model["born_model"]["mu_born_model"]["comp_dep_delta_d"] is True


def test_hubach_cold_start_rejects_removed_option_dict_keys() -> None:
    feed = _row0_feed()
    mix = _hubach_mixture(feed)

    with pytest.raises(epcsaft.InputError, match="force_seed_solve"):
        mix.equilibrium(
            kind="electrolyte_lle",
            T=T_K,
            P=P_PA,
            z=feed,
            options={
                "force_seed_solve": True,
            },
        )


def test_hubach_rejects_removed_density_diagnostics_option() -> None:
    feed = _row0_feed()
    mix = _hubach_mixture(feed)

    with pytest.raises(epcsaft.InputError, match="density_diagnostics"):
        mix.equilibrium(
            kind="electrolyte_lle",
            T=T_K,
            P=P_PA,
            z=feed,
            options={"density_diagnostics": "full"},
        )
