from __future__ import annotations

import csv
import os
from pathlib import Path

import numpy as np
import pytest

import epcsaft
from epcsaft import ePCSAFTMixture
from epcsaft.equilibrium_core.electrolyte_seeds import charge_neutral_lle_seed_from_org_phase
from tests.equilibrium.electrolyte.test_electrolyte_lle_smokes import _assert_electrolyte_lle_route_pending
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


def _row0_initial_phases() -> dict[str, object]:
    with (FIXTURE_DIR / "initial_phase_guesses.csv").open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle))
    aq = np.asarray([float(row["aq_" + name]) for name in SPECIES], dtype=float)
    org = np.asarray([float(row["org_" + name]) for name in SPECIES], dtype=float)
    return {"aq": aq, "org": org, "phase_fraction": float(row["beta_org"])}


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


def test_hubach_seed_helper_constructs_charge_neutral_material_balanced_guess() -> None:
    feed = _row0_feed()
    mix = _hubach_mixture(feed)
    org = np.asarray([0.55, 0.30, 0.10, 0.025, 0.025], dtype=float)

    seed = charge_neutral_lle_seed_from_org_phase(feed, org, 0.05, mix.parameters["z"])
    payload = seed.to_initial_phases()
    aq = payload["aq"]
    org_out = payload["org"]
    beta = payload["phase_fraction"]

    assert_allclose((1.0 - beta) * aq + beta * org_out, feed, atol=1.0e-12)
    assert abs(float(np.dot(aq, mix.parameters["z"]))) <= 1.0e-8
    assert abs(float(np.dot(org_out, mix.parameters["z"]))) <= 1.0e-8
    assert org_out[1] + org_out[2] > aq[1] + aq[2]
    assert aq[0] > org_out[0]


@pytest.mark.skipif(
    os.environ.get("EPCSAFT_RUN_HUBACH_LLE") != "1", reason="Hubach native LLE solve is an opt-in hard-case regression."
)
def test_hubach_row0_explicit_seed_requires_native_ipopt_route() -> None:
    feed = _row0_feed()
    mix = _hubach_mixture(feed)

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="electrolyte_lle",
            T=T_K,
            P=P_PA,
            z=feed,
            initial_phases=_row0_initial_phases(),
            options=epcsaft.EquilibriumOptions(max_iterations=180, tolerance=1.0e-8),
        )

    _assert_electrolyte_lle_route_pending(excinfo)


@pytest.mark.skipif(
    os.environ.get("EPCSAFT_RUN_HUBACH_LLE") != "1", reason="Hubach native LLE solve is an opt-in hard-case regression."
)
def test_hubach_cold_start_requires_native_ipopt_route() -> None:
    feed = _row0_feed()
    mix = _hubach_mixture(feed)

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="electrolyte_lle",
            T=T_K,
            P=P_PA,
            z=feed,
            options=epcsaft.EquilibriumOptions(max_iterations=2, tolerance=1.0e-12),
        )

    _assert_electrolyte_lle_route_pending(excinfo)


@pytest.mark.skipif(
    os.environ.get("EPCSAFT_RUN_HUBACH_LLE") != "1",
    reason="Hubach cold-start failure is an opt-in hard-case regression.",
)
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


@pytest.mark.skipif(
    os.environ.get("EPCSAFT_RUN_HUBACH_LLE") != "1",
    reason="Hubach density diagnostics are an opt-in hard-case regression.",
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
