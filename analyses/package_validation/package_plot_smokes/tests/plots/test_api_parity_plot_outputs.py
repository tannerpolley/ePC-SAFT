from __future__ import annotations

import numpy as np

from tests.helpers.runtime_cases import _ionic_state
from tests.helpers.runtime_cases import _neutral_state
from analyses.package_validation.package_plot_smokes.tests.plots.plot_helpers import save_parity_plot


def test_runtime_alias_canonical_method_parity_plot() -> None:
    state, species = _ionic_state()
    labels = [
        "p",
        "rho",
        "rho_molar",
        "rho_mass",
        "z",
        "ares",
        "dadt",
        "hres",
        "sres",
        "gres",
        "mures water",
        "mures Na",
        "mures Cl",
        "fugcoef water",
        "fugcoef Na",
        "fugcoef Cl",
        "epsr mixture",
        "osmotic_coef",
        "gamma water",
        "gamma Na",
        "gamma Cl",
        "diag pressure",
        "gsolv Na",
        "gsolv Cl",
    ]
    alias_values = [
        state.p(),
        state.rho(),
        state.rho_molar(),
        state.rho_mass(),
        state.z(),
        state.ares(),
        state.dadt(),
        state.hres(),
        state.sres(),
        state.gres(),
        *state.mures().tolist(),
        *state.fugcoef().tolist(),
        state.epsr()[0],
        state.osmotic_coef()[0],
        *[state.gamma(species=species)[label] for label in species],
        state.diag(species=species)["pressure"],
        *[state.gsolv(species=species)[label] for label in ("Na+", "Cl-")],
    ]
    canonical_values = [
        state.pressure(),
        state.density(),
        state.molar_density(),
        state.mass_density(),
        state.compressibility_factor(),
        state.residual_helmholtz(),
        state.temperature_derivative_residual_helmholtz(),
        state.residual_enthalpy(),
        state.residual_entropy(),
        state.residual_gibbs(),
        *state.residual_chemical_potential().tolist(),
        *state.fugacity_coefficient().tolist(),
        state.relative_permittivity()[0],
        state.osmotic_coefficient()[0],
        *[state.activity_coefficient(species=species)[label] for label in species],
        state.state_diagnostics(species=species)["pressure"],
        *[state.solvation_free_energy(species=species)[label] for label in ("Na+", "Cl-")],
    ]
    save_parity_plot(
        "runtime_alias_canonical_method_parity.png",
        "State aliases vs canonical public methods",
        labels,
        np.asarray(alias_values, dtype=float),
        np.asarray(canonical_values, dtype=float),
        category=("api", "parity"),
    )


def test_runtime_diagnostics_public_method_parity_plot() -> None:
    neutral_state, neutral_species = _neutral_state()
    ionic_state, ionic_species = _ionic_state()
    neutral_diag = neutral_state.state_diagnostics(species=neutral_species)
    ionic_diag = ionic_state.state_diagnostics(species=ionic_species)
    labels = [
        "neutral P",
        "neutral rho",
        "neutral Z",
        "neutral ares",
        "neutral mures A",
        "neutral mures B",
        "neutral mures C",
        "neutral phi A",
        "neutral phi B",
        "neutral phi C",
        "ionic P",
        "ionic rho",
        "ionic rho mass",
        "ionic Z",
        "ionic ares",
        "ionic gamma water",
        "ionic gamma Na",
        "ionic gamma Cl",
        "ionic mean gamma mole",
        "ionic mean gamma molality",
        "ionic osmotic",
        "ionic gsolv Na",
        "ionic gsolv Cl",
    ]
    diag_values = [
        neutral_diag["pressure"],
        neutral_diag["density"],
        neutral_diag["compressibility_factor"],
        neutral_diag["residual_helmholtz"],
        *np.asarray(neutral_diag["residual_chemical_potential"], dtype=float).tolist(),
        *np.asarray(neutral_diag["fugacity_coefficient"], dtype=float).tolist(),
        ionic_diag["pressure"],
        ionic_diag["density"],
        ionic_diag["mass_density"],
        ionic_diag["compressibility_factor"],
        ionic_diag["residual_helmholtz"],
        *[ionic_diag["activity_coefficient"][label] for label in ionic_species],
        ionic_diag["mean_ionic_activity_coefficient_mole"]["Na+Cl-"],
        ionic_diag["mean_ionic_activity_coefficient_molality"]["Na+Cl-"],
        np.asarray(ionic_diag["osmotic_coefficient"], dtype=float)[0],
        *[ionic_diag["solvation_free_energy"][label] for label in ("Na+", "Cl-")],
    ]
    public_values = [
        neutral_state.pressure(),
        neutral_state.density(),
        neutral_state.compressibility_factor(),
        neutral_state.residual_helmholtz(),
        *neutral_state.residual_chemical_potential().tolist(),
        *neutral_state.fugacity_coefficient(natural_log=False).tolist(),
        ionic_state.pressure(),
        ionic_state.density(),
        ionic_state.mass_density(),
        ionic_state.compressibility_factor(),
        ionic_state.residual_helmholtz(),
        *[ionic_state.activity_coefficient(species=ionic_species)[label] for label in ionic_species],
        ionic_state.activity_coefficient(species=ionic_species, mean_ionic_form=True, basis="mole")["Na+Cl-"],
        ionic_state.activity_coefficient(species=ionic_species, mean_ionic_form=True, basis="molality")["Na+Cl-"],
        ionic_state.osmotic_coefficient()[0],
        *[ionic_state.solvation_free_energy(species=ionic_species)[label] for label in ("Na+", "Cl-")],
    ]
    save_parity_plot(
        "runtime_diagnostics_public_method_parity.png",
        "State diagnostics payload vs public methods",
        labels,
        np.asarray(diag_values, dtype=float),
        np.asarray(public_values, dtype=float),
        category=("api", "parity"),
    )
