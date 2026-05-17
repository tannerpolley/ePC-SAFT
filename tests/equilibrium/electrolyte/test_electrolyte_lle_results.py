from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest

import epcsaft
from epcsaft import ePCSAFTMixture
from epcsaft.equilibrium_core.electrolyte_basis import (
    build_electrolyte_basis,
    explicit_to_formula_composition,
    formula_to_explicit_composition,
)
from tests.equilibrium.core.test_stability import _assert_stability_route_pending
from tests.helpers.numeric import assert_allclose

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
    formula = explicit_to_formula_composition(feed, basis_dict)
    explicit, _scale = formula_to_explicit_composition(formula, basis_dict, len(species))

    assert payload["salt_pairs"][0]["label"] == "MgCl2"
    assert payload["salt_pairs"][0]["cation_stoich"] == 1
    assert payload["salt_pairs"][0]["anion_stoich"] == 2
    assert_allclose(explicit, feed, atol=1.0e-12)
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
    formula = explicit_to_formula_composition(feed, basis_dict)
    explicit, _scale = formula_to_explicit_composition(formula, basis_dict, len(species))

    assert payload["salt_pairs"][0]["label"] == "Na2SO4"
    assert payload["salt_pairs"][0]["cation_stoich"] == 2
    assert payload["salt_pairs"][0]["anion_stoich"] == 1
    assert_allclose(explicit, feed, atol=1.0e-12)
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

def test_electrolyte_stability_requires_native_ipopt_route_after_validation() -> None:
    mix = _case2_mixture()

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="electrolyte_stability",
            T=298.15,
            P=1.0e5,
            z=_case2_feed(),
            options=epcsaft.EquilibriumOptions(max_iterations=80, tolerance=1.0e-8),
        )

    _assert_stability_route_pending(excinfo, route="electrolyte_stability")
