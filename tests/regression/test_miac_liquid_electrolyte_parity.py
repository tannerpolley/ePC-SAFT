from __future__ import annotations

import numpy as np

from analyses.miac_fits.scripts import validate_miac_fits as vmf
from scripts._epcsaft_oop import as_mixture


def test_figiel_miac_liquid_electrolyte_activity_paths_stay_finite() -> None:
    cases = [
        ("water", "NaBr", (), "Na+Br-"),
        ("methanol", "KCl", (), "K+Cl-"),
        ("water-ethanol", "NaBr", (("water", 0.0), ("ethanol", 1.0)), "Na+Br-"),
    ]

    for solvent_system, salt, comp_signature, pair_key in cases:
        combo = next(
            combo
            for combo in vmf.discover_combos(solvent_scope=solvent_system, salt_scope=salt)
            if combo.get("comp_signature") == comp_signature
        )
        species = vmf._species_for_combo(salt, solvent_system)
        params = vmf.build_params_for_variant("2025_Figiel", combo, user_options={})
        mixture = as_mixture(params, species=species)
        x = vmf._molality_to_molefraction_combo(0.5, salt, solvent_system, dict(combo.get("comp", {})))
        state = mixture.state(T=vmf.T_REF, x=x, P=vmf.P_REF, phase="liq")

        mean_molality = state.activity_coefficient(species=species, mean_ionic_form=True, basis="molality")
        mean_mole = state.activity_coefficient(species=species, mean_ionic_form=True, basis="mole")
        derivatives = state.born_ssmds_liquid_derivatives()

        assert pair_key in mean_molality
        assert pair_key in mean_mole
        assert np.isfinite(mean_molality[pair_key])
        assert np.isfinite(mean_mole[pair_key])
        assert derivatives["supported"] is True
