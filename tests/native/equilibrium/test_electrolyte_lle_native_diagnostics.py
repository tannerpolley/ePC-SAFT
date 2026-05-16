from __future__ import annotations

import numpy as np
import pytest

from epcsaft import _core, ePCSAFTMixture
from tests.helpers.numeric import assert_allclose


def _case2_feed() -> list[float]:
    return [
        0.940373242284748,
        0.04879624542603625,
        0.0019339313461782701,
        0.003481324798429627,
        0.005415256144607897,
    ]


def _case2_mixture() -> ePCSAFTMixture:
    feed = _case2_feed()
    return ePCSAFTMixture.from_dataset("2022_Ascani", ["H2O", "Butanol", "Na+", "K+", "Cl-"], feed, 298.15)


def test_native_electrolyte_lle_reports_basis_and_transform_diagnostics() -> None:
    mix = _case2_mixture()
    payload = _core._solve_equilibrium_native(
        mix._native,
        {
            "kind": "electrolyte_lle",
            "T": 298.15,
            "P": 1.0e5,
            "z": _case2_feed(),
            "species": mix.species,
            "options": {
                "max_iterations": 180,
                "tolerance": 1.0e-8,
                "min_composition": 1.0e-12,
                "include_phase_diagnostics": False,
                "stability_precheck": True,
            },
        },
    )

    diagnostics = payload["diagnostics"]

    assert payload["backend"] == "electrolyte_lle"
    assert payload["split_detected"] is True
    assert diagnostics["basis_model"] == "charge_neutral_salt_pair_coordinates"
    assert diagnostics["basis_vector_model"] == "salt_pair_stoichiometry_rows_by_public_species"
    assert diagnostics["variable_model"] == "ascani_transformed_salt_pairs"
    assert diagnostics["basis_rank"] == 2
    assert diagnostics["neutral_species_indices"] == [0.0, 1.0]
    assert diagnostics["cation_species_indices"] == [2.0, 3.0]
    assert diagnostics["anion_species_indices"] == [4.0]
    assert diagnostics["charged_species_indices"] == [2.0, 3.0, 4.0]
    assert diagnostics["salt_pair_cation_indices"] == [2.0, 3.0]
    assert diagnostics["salt_pair_anion_indices"] == [4.0, 4.0]
    assert diagnostics["salt_pair_cation_stoich"] == [1.0, 1.0]
    assert diagnostics["salt_pair_anion_stoich"] == [1.0, 1.0]
    assert diagnostics["explicit_species_count"] == len(mix.species)
    assert diagnostics["formula_variable_count"] == 4
    assert diagnostics["transformed_variable_count"] == 4
    assert diagnostics["phase_charge_enforced_by_basis"] is True
    assert diagnostics["material_balance_enforced_by_formula_transform"] is True
    assert diagnostics["formula_phase_positivity_enforced_by_transform"] is True
    assert diagnostics["explicit_public_species_reported"] is True
    assert diagnostics["transformed_variable_candidate_available"] is True
    assert diagnostics["accepted_transformed_variables_feasible"] is True

    rows = diagnostics["basis_vector_rows"]
    cols = diagnostics["basis_vector_cols"]
    basis_vectors = np.asarray(diagnostics["basis_vectors_row_major"], dtype=float).reshape(rows, cols)
    species_charges = np.asarray(diagnostics["species_charge_vector"], dtype=float)
    assert_allclose(basis_vectors @ species_charges, np.zeros(rows), atol=1.0e-12)

    assert diagnostics["phase_charge_balance_feed"] == pytest.approx(0.0, abs=1.0e-12)
    assert diagnostics["phase_charge_balance_aq"] == pytest.approx(0.0, abs=1.0e-8)
    assert diagnostics["phase_charge_balance_org"] == pytest.approx(0.0, abs=1.0e-8)
    assert diagnostics["phase_charge_balance_max_abs"] == pytest.approx(
        diagnostics["charge_balance_error"],
        abs=1.0e-12,
    )
    assert len(diagnostics["accepted_aq_formula"]) == diagnostics["formula_variable_count"]
    assert len(diagnostics["accepted_org_formula"]) == diagnostics["formula_variable_count"]
    assert max(abs(value) for value in diagnostics["material_balance_residual"]) <= diagnostics["material_balance_error"]
