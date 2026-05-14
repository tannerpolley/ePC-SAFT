from __future__ import annotations

import numpy as np
import pytest

from analyses.miac_fits.scripts import validate_miac_fits as vmf
from scripts._epcsaft_oop import as_mixture


def _figiel_nabr_water_case():
    combo = next(
        combo
        for combo in vmf.discover_combos(solvent_scope="water", salt_scope="NaBr")
        if not combo.get("comp_signature")
    )
    species = vmf._species_for_combo("NaBr", "water")
    params = vmf.build_params_for_variant("2025_Figiel", combo, user_options={})
    x = vmf._molality_to_molefraction_combo(0.5, "NaBr", "water", dict(combo.get("comp", {})))
    return combo, species, params, x


def test_figiel_2025_born_parameter_fixture_values_are_preserved() -> None:
    _, species, params, _ = _figiel_nabr_water_case()

    d_born = dict(zip(species, np.asarray(params["d_born"], dtype=float)))
    f_solv = dict(zip(species, np.asarray(params["f_solv"], dtype=float)))

    assert d_born["Na+"] == pytest.approx(3.445)
    assert d_born["Br-"] == pytest.approx(4.48)
    assert d_born["H2O-2B-NaCl"] == pytest.approx(0.0)
    assert f_solv["Na+"] == pytest.approx(1.0)
    assert f_solv["Br-"] == pytest.approx(1.0)
    assert f_solv["H2O-2B-NaCl"] == pytest.approx(1.5)

    born_model = params["elec_model"]["born_model"]
    assert born_model["d_Born_mode"] == 3
    assert born_model["solvation_shell_model"] is True
    assert born_model["dielectric_saturation"] is True


def test_figiel_2025_miac_liquid_outputs_and_born_derivatives_remain_finite() -> None:
    _, species, params, x = _figiel_nabr_water_case()
    mixture = as_mixture(params, species=species)
    state = mixture.state(T=vmf.T_REF, x=x, P=vmf.P_REF, phase="liq")

    miac = state.activity_coefficient(species=species, mean_ionic_form=True, basis="molality")
    fugacity_terms = state.fugacity_coefficient(return_contribution_terms=True)
    derivatives = state.born_ssmds_liquid_derivatives()

    assert miac["Na+Br-"] == pytest.approx(0.7732309439080085)
    assert np.all(np.isfinite(fugacity_terms["terms"]["born"]))
    assert derivatives["supported"] is True
    assert derivatives["backend"] in {"analytic", "cppad"}
    assert np.all(np.isfinite(derivatives["a_born_d_d_born"]))
    assert np.all(np.isfinite(derivatives["a_born_d_f_solv"]))
