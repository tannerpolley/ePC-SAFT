from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import numpy as np
import pytest

import epcsaft
from epcsaft import ePCSAFTMixture
from epcsaft.equilibrium_core.electrolyte_seeds import charge_neutral_lle_seed_from_org_phase

REPO_ROOT = Path(__file__).resolve().parents[2]
FIXTURE_DIR = REPO_ROOT / "data" / "equilibrium_benchmarks" / "electrolyte_lle" / "hubach_2024"
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
    np.testing.assert_allclose(mix.parameters["z"], [0.0, 0.0, 0.0, 1.0, -1.0])
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

    np.testing.assert_allclose((1.0 - beta) * aq + beta * org_out, feed, atol=1.0e-12)
    assert abs(float(np.dot(aq, mix.parameters["z"]))) <= 1.0e-8
    assert abs(float(np.dot(org_out, mix.parameters["z"]))) <= 1.0e-8
    assert org_out[1] + org_out[2] > aq[1] + aq[2]
    assert aq[0] > org_out[0]


@pytest.mark.skipif(
    os.environ.get("EPCSAFT_RUN_HUBACH_LLE") != "1", reason="Hubach native LLE solve is an opt-in hard-case regression."
)
def test_hubach_row0_explicit_seed_converges_to_distinct_fixed_species_lle() -> None:
    feed = _row0_feed()
    mix = _hubach_mixture(feed)

    result = mix.equilibrium(
        kind="electrolyte_lle",
        T=T_K,
        P=P_PA,
        z=feed,
        initial_phases=_row0_initial_phases(),
        options=epcsaft.EquilibriumOptions(max_iterations=180, tolerance=1.0e-8, damping=0.5),
    )

    diagnostics = result.diagnostics
    phases = {phase.label: phase for phase in result.phases}
    reconstructed = sum(phase.phase_fraction * phase.composition for phase in result.phases)

    assert result.split_detected is True
    assert len(result.phases) == 2
    assert diagnostics["phase_equilibrium_model"] == "electrolyte_lle_v5_native_charge_constrained_solve"
    assert diagnostics["acceptance_gate"] == "predictive_nonlinear_solve"
    assert diagnostics["solver_residual_norm"] <= 1.0e-6
    assert diagnostics["material_balance_error"] <= 1.0e-10
    assert diagnostics["charge_balance_error"] <= 1.0e-8
    assert diagnostics["phase_distance"] > 1.0e-4
    assert diagnostics["seed_attempt_count"] >= 1
    assert diagnostics["seed_attempts"]
    assert phases["aq"].composition[0] > phases["org"].composition[0]
    assert (
        phases["org"].composition[1] + phases["org"].composition[2]
        > phases["aq"].composition[1] + phases["aq"].composition[2]
    )
    np.testing.assert_allclose(reconstructed, feed, atol=1.0e-10)
    json.dumps(result.to_dict(), allow_nan=False)


@pytest.mark.skipif(
    os.environ.get("EPCSAFT_RUN_HUBACH_LLE") != "1", reason="Hubach native LLE solve is an opt-in hard-case regression."
)
def test_hubach_cold_start_failure_reports_seed_attempts_before_error() -> None:
    feed = _row0_feed()
    mix = _hubach_mixture(feed)

    try:
        mix.equilibrium(
            kind="electrolyte_lle",
            T=T_K,
            P=P_PA,
            z=feed,
            options=epcsaft.EquilibriumOptions(max_iterations=2, tolerance=1.0e-12),
        )
    except epcsaft.SolutionError as exc:
        diagnostics = exc.args[1]
    else:
        return

    seed_names = [attempt["seed_name"] for attempt in diagnostics["seed_attempts"]]
    assert any(name.startswith("native_tpd_trial") or name.startswith("native_gibbs_tpd_trial") for name in seed_names)
    assert any("org_endpoint" in name for name in seed_names)
    assert any("_b5" in name or "_b10" in name for name in seed_names)
    assert diagnostics["seed_attempt_count"] == len(seed_names)
    assert diagnostics["seed_attempt_count"] > 4


@pytest.mark.skipif(
    os.environ.get("EPCSAFT_RUN_HUBACH_LLE") != "1",
    reason="Hubach legacy-candidate fallback is an opt-in hard-case regression.",
)
def test_hubach_cold_start_preserves_distinct_legacy_candidate_on_strict_failure() -> None:
    feed = _row0_feed()
    mix = _hubach_mixture(feed)

    try:
        mix.equilibrium(
            kind="electrolyte_lle",
            T=T_K,
            P=P_PA,
            z=feed,
            options={
                "max_nfev": 20,
                "solver_tol": 1.0e-8,
                "split_tol": 1.0e-4,
                "solver_accept_norm": 0.5,
                "tpdf_global_trials": 1200,
                "tpdf_local_trials": 600,
                "charge_weight": 1000.0,
                "force_seed_solve": True,
            },
        )
    except epcsaft.SolutionError as exc:
        diagnostics = exc.args[1]
    else:
        return

    assert diagnostics["acceptance_gate"] == "predictive_solve_failed"
    assert diagnostics["legacy_candidate_found"] is True
    assert diagnostics["legacy_candidate_phase_distance"] > diagnostics["legacy_candidate_split_tolerance"]
    assert diagnostics["legacy_candidate_material_balance_error"] <= 1.0e-8
    assert diagnostics["legacy_candidate_charge_balance_error"] <= 1.0e-6
    assert diagnostics["legacy_candidate_rejected_by_strict_gate"] is True
    assert len(diagnostics["legacy_candidate_aq_composition"]) == mix.ncomp
    assert len(diagnostics["legacy_candidate_org_composition"]) == mix.ncomp
    assert diagnostics["density_failure_count"] >= 0
    assert diagnostics["density_diagnostics_mode"] == "auto"
    assert diagnostics["density_validity_gate"] in {"passed", "failed", "not_evaluated"}
    assert "tpdf_global_trials" in diagnostics["ignored_legacy_options"]
    assert "force_seed_solve" in diagnostics["ignored_legacy_options"]
    json.dumps(diagnostics, allow_nan=False)


@pytest.mark.skipif(
    os.environ.get("EPCSAFT_RUN_HUBACH_LLE") != "1",
    reason="Hubach density diagnostics are an opt-in hard-case regression.",
)
def test_hubach_cold_start_density_failure_payload_is_json_safe() -> None:
    feed = _row0_feed()
    mix = _hubach_mixture(feed)

    try:
        mix.equilibrium(
            kind="electrolyte_lle",
            T=T_K,
            P=P_PA,
            z=feed,
            options=epcsaft.EquilibriumOptions(max_iterations=2, tolerance=1.0e-12, density_diagnostics="full"),
        )
    except epcsaft.SolutionError as exc:
        diagnostics = exc.args[1]
    else:
        return

    assert diagnostics["density_diagnostics_mode"] == "full"
    assert diagnostics["density_failure_count"] >= 0
    assert isinstance(diagnostics["density_failure_contexts"], list)
    for context in diagnostics["density_failure_contexts"]:
        assert {"phase_label", "phase_kind", "T", "P", "composition", "scan_point_count", "rejection_reason"} <= set(
            context
        )
        assert len(context["composition"]) == mix.ncomp
    json.dumps(diagnostics, allow_nan=False)
