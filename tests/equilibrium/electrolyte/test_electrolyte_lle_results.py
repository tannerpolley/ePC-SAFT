from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest

import epcsaft
import epcsaft.ipopt_backend as ipopt_backend
from epcsaft import ePCSAFTMixture
from epcsaft.equilibrium import _explicit_to_formula_composition, _formula_to_explicit_composition
from epcsaft.equilibrium_core.electrolyte_basis import build_electrolyte_basis

REPO_ROOT = Path(__file__).resolve().parents[3]

ASCANI_CASE2 = REPO_ROOT / "data" / "multiphase" / "ascani_case2_model_comparison.csv"

def _ascani_water_butanol_nacl_mixture(feed=None) -> ePCSAFTMixture:
    species = ["H2O", "Butanol", "Na+", "Cl-"]
    if feed is None:
        feed = np.asarray([0.55, 0.40, 0.025, 0.025], dtype=float)
    return ePCSAFTMixture.from_dataset("2022_Ascani", species, feed, 298.15)

def _ascani_case2_reference() -> dict[str, float]:
    rows: dict[str, float] = {}
    with ASCANI_CASE2.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            value = row["model_2020"]
            if value and value.lower() not in {"true", "false"}:
                rows[row["quantity"]] = float(value)
            paper = row["paper"]
            if paper and paper.lower() not in {"true", "false"}:
                rows[row["quantity"] + " paper"] = float(paper)
    return rows

def _expand_formula_phase(water: float, butanol: float, nacl: float, kcl: float) -> np.ndarray:
    formula = np.asarray([water, butanol, nacl, kcl], dtype=float)
    denom = float(formula[0] + formula[1] + 2.0 * formula[2] + 2.0 * formula[3])
    return np.asarray([formula[0], formula[1], formula[2], formula[3], formula[2] + formula[3]], dtype=float) / denom

def _case2_feed() -> np.ndarray:
    return np.asarray(
        [
            0.940373242284748,
            0.04879624542603625,
            0.0019339313461782701,
            0.003481324798429627,
            0.005415256144607897,
        ],
        dtype=float,
    )

def _case2_mixture(feed=None) -> ePCSAFTMixture:
    if feed is None:
        feed = _case2_feed()
    return ePCSAFTMixture.from_dataset("2022_Ascani", ["H2O", "Butanol", "Na+", "K+", "Cl-"], feed, 298.15)

def test_ascani_counterion_basis_has_expected_rank_and_preserves_charge() -> None:
    mix = _case2_mixture()
    basis = build_electrolyte_basis(mix.species, mix.parameters["z"], _case2_feed())

    assert basis.variable_model == "ascani_transformed_salt_pairs"
    assert basis.rank == 2
    assert basis.e_matrix.shape == (2, 3)
    assert [pair["label"] for pair in basis.salt_pairs] == ["NaCl", "KCl"]

    xi = np.asarray([0.001, 0.002], dtype=float)
    charged_delta = basis.e_matrix.T @ xi
    charged_charges = np.asarray([mix.parameters["z"][i] for i in basis.charged_indices], dtype=float)
    assert abs(float(np.dot(charged_delta, charged_charges))) <= 1.0e-12

def test_divalent_two_to_one_salt_basis_reconstructs_charge_neutral_formula() -> None:
    species = ["H2O", "TBP", "Mg2+", "Cl-"]
    charges = np.asarray([0.0, 0.0, 2.0, -1.0], dtype=float)
    feed = np.asarray([0.80, 0.10, 0.02, 0.04], dtype=float)
    feed = feed / float(np.sum(feed))

    basis = build_electrolyte_basis(species, charges, feed)
    payload = basis.to_dict()
    basis_dict = {
        "neutral_indices": tuple(payload["neutral_indices"]),
        "salt_pairs": tuple(payload["salt_pairs"]),
    }
    formula = _explicit_to_formula_composition(feed, basis_dict)
    explicit, _scale = _formula_to_explicit_composition(formula, basis_dict, len(species))

    assert payload["salt_pairs"][0]["label"] == "MgCl2"
    assert payload["salt_pairs"][0]["cation_stoich"] == 1
    assert payload["salt_pairs"][0]["anion_stoich"] == 2
    np.testing.assert_allclose(explicit, feed, atol=1.0e-12)
    assert abs(float(np.dot(explicit, charges))) <= 1.0e-12

def test_one_to_two_salt_basis_reconstructs_charge_neutral_formula() -> None:
    species = ["H2O", "TBP", "Na+", "SO4--"]
    charges = np.asarray([0.0, 0.0, 1.0, -2.0], dtype=float)
    feed = np.asarray([0.80, 0.10, 0.04, 0.02], dtype=float)
    feed = feed / float(np.sum(feed))

    basis = build_electrolyte_basis(species, charges, feed)
    payload = basis.to_dict()
    basis_dict = {
        "neutral_indices": tuple(payload["neutral_indices"]),
        "salt_pairs": tuple(payload["salt_pairs"]),
    }
    formula = _explicit_to_formula_composition(feed, basis_dict)
    explicit, _scale = _formula_to_explicit_composition(formula, basis_dict, len(species))

    assert payload["salt_pairs"][0]["label"] == "Na2SO4"
    assert payload["salt_pairs"][0]["cation_stoich"] == 2
    assert payload["salt_pairs"][0]["anion_stoich"] == 1
    np.testing.assert_allclose(explicit, feed, atol=1.0e-12)
    assert abs(float(np.dot(explicit, charges))) <= 1.0e-12

def test_mixed_monovalent_divalent_shared_anion_basis_builds() -> None:
    species = ["H2O", "TBP", "Li+", "Mg2+", "Cl-"]
    charges = np.asarray([0.0, 0.0, 1.0, 2.0, -1.0], dtype=float)
    feed = np.asarray([0.80, 0.10, 0.02, 0.01, 0.04], dtype=float)
    feed = feed / float(np.sum(feed))

    basis = build_electrolyte_basis(species, charges, feed)
    labels = [pair["label"] for pair in basis.salt_pairs]

    assert labels == ["LiCl", "MgCl2"]
    assert basis.rank == 2

def test_electrolyte_stability_payload_is_json_serializable() -> None:
    mix = _case2_mixture()

    result = mix.equilibrium(
        kind="electrolyte_stability",
        T=298.15,
        P=1.0e5,
        z=_case2_feed(),
        options=epcsaft.EquilibriumOptions(max_iterations=80, tolerance=1.0e-8),
    )

    assert result.backend == "electrolyte_tpd"
    assert result.diagnostics["variable_model"] == "ascani_transformed_salt_pairs"
    assert result.diagnostics["basis_rank"] == 2
    assert result.diagnostics["tpd_method"] == "native_tpd_global_search"
    assert result.diagnostics["solver_language"] == "c++"
    assert result.diagnostics["tpd_objective_value"] == pytest.approx(result.min_tpd)
    assert result.diagnostics["phase_charge_balance"]["trial"] == pytest.approx(0.0, abs=1.0e-8)
    json.dumps(result.to_dict(), allow_nan=False)

def test_mixed_electrolyte_lle_reports_phase_charge_balance_for_distributed_ions() -> None:
    mix = _case2_mixture()
    result = mix.equilibrium(
        kind="electrolyte_lle",
        T=298.15,
        P=1.0e5,
        z=_case2_feed(),
        options=epcsaft.EquilibriumOptions(max_iterations=180, tolerance=1.0e-8, damping=0.5),
    )

    charge_balance = result.diagnostics["phase_charge_balance"]
    assert set(charge_balance) == {"feed", "aq", "org", "max_abs"}
    assert charge_balance["feed"] == pytest.approx(0.0, abs=1.0e-12)
    assert charge_balance["aq"] == pytest.approx(0.0, abs=1.0e-8)
    assert charge_balance["org"] == pytest.approx(0.0, abs=1.0e-8)
    assert charge_balance["max_abs"] == pytest.approx(result.diagnostics["charge_balance_error"], abs=1.0e-12)

    salt_pairs = result.diagnostics["salt_pairs"]
    assert [pair["label"] for pair in salt_pairs] == ["NaCl", "KCl"]
    assert [(pair["cation"], pair["anion"]) for pair in salt_pairs] == [(2, 4), (3, 4)]
    aq, org = result.phases
    for ion_index in (2, 3, 4):
        assert aq.composition[ion_index] > 0.0
        assert org.composition[ion_index] > 0.0
