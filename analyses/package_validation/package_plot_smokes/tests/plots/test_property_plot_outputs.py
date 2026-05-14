from __future__ import annotations

import numpy as np

from tests.helpers.runtime_cases import _ionic_state
from tests.helpers.runtime_cases import _neutral_state
from analyses.package_validation.package_plot_smokes.tests.plots.plot_helpers import save_comparison_plot


def test_runtime_neutral_scalar_reference_comparison_plot() -> None:
    state, _species = _neutral_state()
    expected_values = {
        "rho": 14330.417110,
        "P": 1276374.1152948933,
        "Z": 0.04594621208078564,
        "ares": -3.54988545131505,
        "dadt": 0.03077401856781036,
        "hres": -15758.229958475444,
        "sres": -55.751451436621096,
        "gres": -2759.779056027235,
    }
    actual_values = {
        "rho": state.density(),
        "P": state.pressure(),
        "Z": state.compressibility_factor(),
        "ares": state.residual_helmholtz(),
        "dadt": state.temperature_derivative_residual_helmholtz(),
        "hres": state.residual_enthalpy(),
        "sres": state.residual_entropy(),
        "gres": state.residual_gibbs(),
    }

    labels = list(expected_values)
    save_comparison_plot(
        "runtime_neutral_scalar_reference_comparison.png",
        "Neutral runtime scalar outputs vs pinned expected values",
        labels,
        np.asarray([actual_values[label] for label in labels], dtype=float),
        np.asarray([expected_values[label] for label in labels], dtype=float),
        category=("properties", "residual_energy"),
    )


def test_runtime_residual_energy_family_reference_comparison_plot() -> None:
    state, _species = _neutral_state()
    expected_values = {
        "Z": 0.04594621208078564,
        "ares": -3.54988545131505,
        "dadt": 0.03077401856781036,
        "hres": -15758.229958475444,
        "sres": -55.751451436621096,
        "gres": -2759.779056027235,
    }
    actual_values = {
        "Z": state.z(),
        "ares": state.ares(),
        "dadt": state.dadt(),
        "hres": state.hres(),
        "sres": state.sres(),
        "gres": state.gres(),
    }
    labels = list(expected_values)
    save_comparison_plot(
        "runtime_residual_energy_property_family.png",
        "Residual energy property family vs pinned expected values",
        labels,
        np.asarray([actual_values[label] for label in labels], dtype=float),
        np.asarray([expected_values[label] for label in labels], dtype=float),
        category=("properties", "residual_energy"),
    )


def test_runtime_ionic_residual_energy_family_reference_comparison_plot() -> None:
    state, _species = _ionic_state()
    expected_values = {
        "Z": 0.000728884077611683,
        "ares": -9.7214027218058,
        "dadt": 0.032388021640507005,
        "hres": -26415.160790583413,
        "sres": -59.523895812302186,
        "gres": -8668.111254145517,
    }
    actual_values = {
        "Z": state.z(),
        "ares": state.ares(),
        "dadt": state.dadt(),
        "hres": state.hres(),
        "sres": state.sres(),
        "gres": state.gres(),
    }
    labels = list(expected_values)
    save_comparison_plot(
        "runtime_ionic_residual_energy_property_family.png",
        "Ionic residual energy property family vs pinned expected values",
        labels,
        np.asarray([actual_values[label] for label in labels], dtype=float),
        np.asarray([expected_values[label] for label in labels], dtype=float),
        category=("properties", "residual_energy"),
    )


def test_runtime_ionic_diagnostics_reference_comparison_plot() -> None:
    state, species = _ionic_state()
    component_activity = state.activity_coefficient(species=species)
    solvation = state.solvation_free_energy(species=species)
    expected_values = {
        "water gamma": 1.0000051724037697,
        "Na gamma": 0.9222113778654043,
        "Cl gamma": 0.9222258090371313,
        "mean gamma x": 0.9222185934230398,
        "mean gamma m": 0.9220341497043553,
        "Na gsolv": -475461.4260703414,
        "Cl gsolv": -489572.50284416083,
        "rho": 55344.274540081075,
        "rho mass": 997.1665703121223,
        "Z": 0.000728884077611683,
    }
    actual_values = {
        "water gamma": component_activity["water"],
        "Na gamma": component_activity["Na+"],
        "Cl gamma": component_activity["Cl-"],
        "mean gamma x": state.activity_coefficient(species=species, mean_ionic_form=True, basis="mole")["Na+Cl-"],
        "mean gamma m": state.activity_coefficient(species=species, mean_ionic_form=True, basis="molality")["Na+Cl-"],
        "Na gsolv": solvation["Na+"],
        "Cl gsolv": solvation["Cl-"],
        "rho": state.density(),
        "rho mass": state.mass_density(),
        "Z": state.compressibility_factor(),
    }

    labels = list(expected_values)
    save_comparison_plot(
        "runtime_ionic_diagnostics_reference_comparison.png",
        "Ionic runtime diagnostics vs pinned expected values",
        labels,
        np.asarray([actual_values[label] for label in labels], dtype=float),
        np.asarray([expected_values[label] for label in labels], dtype=float),
        category=("properties", "activity_fugacity"),
    )


def test_runtime_neutral_fugacity_chemical_potential_comparison_plot() -> None:
    state, _species = _neutral_state()
    mures = state.mures()
    lnfug = state.fugcoef()
    expected_values = {
        "mures A": -1.1478687523834008,
        "mures B": -3.6543804288405415,
        "mures C": -5.488063725572939,
        "lnphi A": 1.9324151168689134,
        "lnphi B": -0.5740965595882255,
        "lnphi C": -2.407779856320623,
    }
    actual_values = {
        "mures A": mures[0],
        "mures B": mures[1],
        "mures C": mures[2],
        "lnphi A": lnfug[0],
        "lnphi B": lnfug[1],
        "lnphi C": lnfug[2],
    }
    labels = list(expected_values)
    save_comparison_plot(
        "runtime_neutral_fugacity_chemical_potential_comparison.png",
        "Neutral fugacity and chemical potential outputs vs pinned expected values",
        labels,
        np.asarray([actual_values[label] for label in labels], dtype=float),
        np.asarray([expected_values[label] for label in labels], dtype=float),
        category=("properties", "activity_fugacity"),
    )


def test_runtime_activity_fugacity_solution_property_comparison_plot() -> None:
    state, species = _ionic_state()
    component_activity = state.gamma(species=species)
    lnfug = state.fugacity_coefficient()
    solvation = state.gsolv(species=species)
    expected_values = {
        "lnphi water": -3.458424439279275,
        "lnphi Na": -191.8799615776576,
        "lnphi Cl": -197.57230810636238,
        "gamma water": 1.0000051724037697,
        "gamma Na": 0.9222113778654043,
        "gamma Cl": 0.9222258090371313,
        "osmotic": 0.9739566103279091,
        "epsr": 78.075982,
        "gsolv Na": -475461.4260703414,
        "gsolv Cl": -489572.50284416083,
    }
    actual_values = {
        "lnphi water": lnfug[0],
        "lnphi Na": lnfug[1],
        "lnphi Cl": lnfug[2],
        "gamma water": component_activity["water"],
        "gamma Na": component_activity["Na+"],
        "gamma Cl": component_activity["Cl-"],
        "osmotic": state.osmotic_coef()[0],
        "epsr": state.epsr()[0],
        "gsolv Na": solvation["Na+"],
        "gsolv Cl": solvation["Cl-"],
    }
    labels = list(expected_values)
    save_comparison_plot(
        "runtime_activity_fugacity_solution_property_comparison.png",
        "Activity, fugacity, and solution diagnostics vs pinned expected values",
        labels,
        np.asarray([actual_values[label] for label in labels], dtype=float),
        np.asarray([expected_values[label] for label in labels], dtype=float),
        category=("properties", "activity_fugacity"),
    )


def test_runtime_method_surface_neutral_reference_comparison_plot() -> None:
    state, _species = _neutral_state()
    mures = state.residual_chemical_potential()
    lnfug = state.fugacity_coefficient()
    fugcoef = state.fugacity_coefficient(natural_log=False)
    expected_values = {
        "rho": 14330.417110,
        "P": 1276374.1152948933,
        "Z": 0.04594621208078564,
        "ares": -3.54988545131505,
        "dadt": 0.03077401856781036,
        "hres": -15758.229958475444,
        "sres": -55.751451436621096,
        "gres": -2759.779056027235,
        "mures A": -1.1478687523834008,
        "mures B": -3.6543804288405415,
        "mures C": -5.488063725572939,
        "lnphi A": 1.9324151168689134,
        "lnphi B": -0.5740965595882255,
        "lnphi C": -2.407779856320623,
        "phi A": 6.906169322700795,
        "phi B": 0.5632134688356544,
        "phi C": 0.09001491894620331,
    }
    actual_values = {
        "rho": state.density(),
        "P": state.pressure(),
        "Z": state.compressibility_factor(),
        "ares": state.residual_helmholtz(),
        "dadt": state.temperature_derivative_residual_helmholtz(),
        "hres": state.residual_enthalpy(),
        "sres": state.residual_entropy(),
        "gres": state.residual_gibbs(),
        "mures A": mures[0],
        "mures B": mures[1],
        "mures C": mures[2],
        "lnphi A": lnfug[0],
        "lnphi B": lnfug[1],
        "lnphi C": lnfug[2],
        "phi A": fugcoef[0],
        "phi B": fugcoef[1],
        "phi C": fugcoef[2],
    }
    labels = list(expected_values)
    save_comparison_plot(
        "runtime_method_surface_neutral_reference_comparison.png",
        "Neutral public method surface vs pinned expected values",
        labels,
        np.asarray([actual_values[label] for label in labels], dtype=float),
        np.asarray([expected_values[label] for label in labels], dtype=float),
        category=("properties", "residual_energy"),
    )


def test_runtime_method_surface_ionic_reference_comparison_plot() -> None:
    state, species = _ionic_state()
    component_activity = state.activity_coefficient(species=species)
    solvation = state.solvation_free_energy(species=species)
    mures = state.residual_chemical_potential()
    lnfug = state.fugacity_coefficient()
    expected_values = {
        "P": 100000.0,
        "rho": 55344.274540081075,
        "rho mass": 997.1665703121223,
        "Z": 0.000728884077611683,
        "ares": -9.7214027218058,
        "dadt": 0.032388021640507005,
        "hres": -26415.160790583413,
        "sres": -59.523895812302186,
        "gres": -8668.111254145517,
        "mures water": -10.682420304620588,
        "mures Na": -199.10395742942775,
        "mures Cl": -204.79630395556683,
        "lnphi water": -3.458424439279275,
        "lnphi Na": -191.8799615776576,
        "lnphi Cl": -197.57230810636238,
        "gamma water": 1.0000051724037697,
        "gamma Na": 0.9222113778654043,
        "gamma Cl": 0.9222258090371313,
        "epsr": 78.075982,
        "osmotic": 0.9739566103279091,
        "gsolv Na": -475461.4260703414,
        "gsolv Cl": -489572.50284416083,
    }
    actual_values = {
        "P": state.pressure(),
        "rho": state.density(),
        "rho mass": state.mass_density(),
        "Z": state.compressibility_factor(),
        "ares": state.residual_helmholtz(),
        "dadt": state.temperature_derivative_residual_helmholtz(),
        "hres": state.residual_enthalpy(),
        "sres": state.residual_entropy(),
        "gres": state.residual_gibbs(),
        "mures water": mures[0],
        "mures Na": mures[1],
        "mures Cl": mures[2],
        "lnphi water": lnfug[0],
        "lnphi Na": lnfug[1],
        "lnphi Cl": lnfug[2],
        "gamma water": component_activity["water"],
        "gamma Na": component_activity["Na+"],
        "gamma Cl": component_activity["Cl-"],
        "epsr": state.relative_permittivity()[0],
        "osmotic": state.osmotic_coefficient()[0],
        "gsolv Na": solvation["Na+"],
        "gsolv Cl": solvation["Cl-"],
    }
    labels = list(expected_values)
    save_comparison_plot(
        "runtime_method_surface_ionic_reference_comparison.png",
        "Ionic public method surface vs pinned expected values",
        labels,
        np.asarray([actual_values[label] for label in labels], dtype=float),
        np.asarray([expected_values[label] for label in labels], dtype=float),
        category=("properties", "activity_fugacity"),
    )
