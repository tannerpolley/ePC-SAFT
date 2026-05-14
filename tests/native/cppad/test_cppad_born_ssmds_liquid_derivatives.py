from __future__ import annotations

import numpy as np
import pytest

from analyses.miac_fits.scripts import validate_miac_fits as vmf
from scripts._epcsaft_oop import as_mixture


def _figiel_nabr_water_state(phase: str = "liq"):
    combo = next(
        combo
        for combo in vmf.discover_combos(solvent_scope="water", salt_scope="NaBr")
        if not combo.get("comp_signature")
    )
    species = vmf._species_for_combo("NaBr", "water")
    params = vmf.build_params_for_variant("2025_Figiel", combo, user_options={})
    mixture = as_mixture(params, species=species)
    x = vmf._molality_to_molefraction_combo(0.5, "NaBr", "water", dict(combo.get("comp", {})))
    state = mixture.state(T=vmf.T_REF, x=x, P=vmf.P_REF, phase=phase)
    return species, state


def _matrix(payload: dict[str, object], key: str) -> np.ndarray:
    return np.asarray(payload[key], dtype=float)


def test_liquid_ssmds_born_derivatives_are_supported_for_figiel_d_born_and_f_solv() -> None:
    species, state = _figiel_nabr_water_state("liq")

    payload = state.born_ssmds_liquid_derivatives()

    assert payload["supported"] is True
    assert payload["backend"] in {"analytic", "cppad"}
    assert payload["phase_scope"] == "liquid_electrolyte_only"
    assert payload["parameters"] == ("d_born", "f_solv")
    assert payload["vapor_support"] is False
    assert "finite" not in str(payload["backend"]).lower()

    d_born = np.asarray(payload["a_born_d_d_born"], dtype=float)
    f_solv = np.asarray(payload["a_born_d_f_solv"], dtype=float)
    assert d_born.shape == (len(species),)
    assert f_solv.shape == (len(species),)
    assert np.all(np.isfinite(d_born))
    assert np.all(np.isfinite(f_solv))
    assert d_born[species.index("Na+")] != pytest.approx(0.0)
    assert d_born[species.index("Br-")] != pytest.approx(0.0)
    assert d_born[species.index("H2O-2B-NaCl")] == pytest.approx(0.0)
    assert f_solv[species.index("Na+")] == pytest.approx(0.0)
    assert f_solv[species.index("Br-")] == pytest.approx(0.0)
    assert f_solv[species.index("H2O-2B-NaCl")] != pytest.approx(0.0)

    for key in (
        "mu_res_d_d_born",
        "mu_res_d_f_solv",
        "lnfug_d_d_born",
        "lnfug_d_f_solv",
        "lngamma_d_d_born",
        "lngamma_d_f_solv",
    ):
        values = _matrix(payload, key)
        assert values.shape == (len(species), len(species))
        assert np.all(np.isfinite(values))

    np.testing.assert_allclose(payload["lnfug_d_d_born"], payload["mu_res_d_d_born"])
    np.testing.assert_allclose(payload["lnfug_d_f_solv"], payload["mu_res_d_f_solv"])


def test_vapor_ssmds_born_derivatives_report_backend_unavailable() -> None:
    _, state = _figiel_nabr_water_state("vap")

    payload = state.born_ssmds_liquid_derivatives()

    assert payload["supported"] is False
    assert payload["backend"] == "backend_unavailable"
    assert payload["message"] == "SSM+DS Born derivatives are liquid-electrolyte only"
