from __future__ import annotations

from scripts.multiphase_model_analysis import ascani_case2_ghat_audit as audit


def test_paper_case2_reported_phase_data_behave_like_mass_fractions():
    payload = audit.build_audit()

    beta = payload["componentwise_mass_split_beta"]
    assert abs(beta["water"] - beta["butanol"]) < 5.0e-5
    assert abs(beta["water"] - beta["KCl"]) < 5.0e-4
    assert abs(beta["water"] - beta["NaCl"]) < 5.0e-3

    molefit = payload["formula_molefraction_fit_if_reported_values_taken_literally"]
    assert molefit["max_abs_residual"] > 1.0e-2


def test_paper_case2_literal_ghat_does_not_match_printed_value():
    payload = audit.build_audit()

    literal = payload["reconstructions"]["literal_stated_formula"]["ghat_eq_j_per_mol"]
    nearest = payload["reconstructions"]["mass_fraction_sum_double_count_salts"]["ghat_eq_j_per_mol"]

    assert abs(literal - audit.PAPER_TARGETS["ghat_eq_j_per_mol"]) > 1.0e4
    assert abs(nearest - audit.PAPER_TARGETS["ghat_eq_j_per_mol"]) < 1.0e3


def test_general_ionic_reconstruction_is_unique_and_not_equal_to_paper_value():
    payload = audit.build_audit()

    ionic = payload["general_ionic_reconstruction"]

    assert abs(ionic["ghat_eq_j_per_mol"] - audit.PAPER_TARGETS["ghat_eq_j_per_mol"]) > 1.0e4
    assert ionic["ghat_delta_eq_minus_feed_j_per_mol"] > 0.0
