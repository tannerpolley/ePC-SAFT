from __future__ import annotations

import numpy as np

from scripts.multiphase_model_analysis.external_package_case2_comparison import (
    WORKED_EXAMPLE,
    charge_residual,
    map_species_to_pseudo_salt_basis,
    mean_ionic_lnfugacities,
)


def test_charge_residual_balanced_species_basis_is_zero() -> None:
    x = np.array([0.90, 0.05, 0.01, 0.015, 0.025], dtype=float)
    assert abs(charge_residual(x, WORKED_EXAMPLE.charges)) < 1.0e-12


def test_map_species_to_pseudo_salt_basis_drops_chloride_and_renormalizes() -> None:
    x = np.array([0.90, 0.05, 0.01, 0.015, 0.025], dtype=float)
    pseudo = map_species_to_pseudo_salt_basis(x, WORKED_EXAMPLE)
    assert set(pseudo) == {"water", "alcohol", "NaCl", "KCl"}
    assert abs(sum(pseudo.values()) - 1.0) < 1.0e-12
    assert abs(pseudo["water"] - 0.90 / 0.975) < 1.0e-12
    assert abs(pseudo["NaCl"] - 0.01 / 0.975) < 1.0e-12
    assert abs(pseudo["KCl"] - 0.015 / 0.975) < 1.0e-12


def test_mean_ionic_lnfugacities_uses_cation_chloride_average() -> None:
    lnf = np.array([-3.5, -5.1, -225.0, -207.0, -224.0], dtype=float)
    mean_values = mean_ionic_lnfugacities(lnf, WORKED_EXAMPLE)
    assert abs(mean_values["NaCl"] - (-224.5)) < 1.0e-12
    assert abs(mean_values["KCl"] - (-215.5)) < 1.0e-12
