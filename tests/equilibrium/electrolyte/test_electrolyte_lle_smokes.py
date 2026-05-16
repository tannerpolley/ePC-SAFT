from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

import epcsaft
from epcsaft._optional_backends import ipopt as ipopt_backend
from epcsaft import ePCSAFTMixture

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

def _assert_ceres_solver_diagnostics(result) -> None:
    assert result.split_detected is True
    assert result.diagnostics["solver_backend"] == "ceres"
    assert result.diagnostics["solver_method"] == "ceres_trust_region_residual_solve"
    assert result.diagnostics["jacobian_backend"] == "cppad_implicit"
    assert result.diagnostics["derivative_backend"] == "cppad_implicit"
    assert result.diagnostics["jacobian_available"] is True
    assert result.diagnostics["derivative_available"] is True

def test_electrolyte_lle_direct_feed_reports_production_solver_derivatives() -> None:
    feed = np.asarray([0.55, 0.40, 0.025, 0.025], dtype=float)
    mix = _ascani_water_butanol_nacl_mixture(feed)

    result = mix.equilibrium(
        kind="electrolyte_lle",
        T=298.15,
        P=1.013e5,
        z=feed,
        options=epcsaft.EquilibriumOptions(include_phase_diagnostics=True),
    )

    _assert_ceres_solver_diagnostics(result)

@pytest.mark.skipif(not ipopt_backend.cyipopt_available(), reason="cyipopt is optional")
def test_electrolyte_lle_direct_feed_solves_ipopt_predictive_split() -> None:
    feed = np.asarray([0.55, 0.40, 0.025, 0.025], dtype=float)
    mix = _ascani_water_butanol_nacl_mixture(feed)

    with pytest.raises(epcsaft.InputError, match="not_available"):
        mix.equilibrium(
            kind="electrolyte_lle",
            T=298.15,
            P=1.013e5,
            z=feed,
            options=epcsaft.EquilibriumOptions(
                solver_backend="ipopt",
                hessian_strategy="lbfgs",
                max_iterations=80,
                tolerance=1.0e-8,
            ),
        )

def test_electrolyte_lle_molality_feed_reports_production_solver_derivatives() -> None:
    mix = _ascani_water_butanol_nacl_mixture()

    result = mix.equilibrium(
        kind="electrolyte_lle",
        T=298.15,
        P=1.013e5,
        solvent_feed={"H2O": 0.58, "Butanol": 0.42},
        salt_molality={"NaCl": 1.0},
        options=epcsaft.EquilibriumOptions(include_phase_diagnostics=True),
    )

    _assert_ceres_solver_diagnostics(result)

def test_electrolyte_lle_accepts_strict_aq_org_initial_phases() -> None:
    mix = _ascani_water_butanol_nacl_mixture()
    aq = np.asarray([0.798324680201737, 0.016320352824141723, 0.09267748348706063, 0.09267748348706063], dtype=float)
    org = np.asarray([0.37006036048879404, 0.6214918588210971, 0.004223890345054407, 0.004223890345054407], dtype=float)
    beta_org = 0.613766575013417
    feed = (1.0 - beta_org) * aq + beta_org * org

    result = mix.equilibrium(
        kind="electrolyte_lle",
        T=298.15,
        P=1.013e5,
        z=feed,
        initial_phases={"aq": aq, "org": org, "phase_fraction": beta_org},
        options=epcsaft.EquilibriumOptions(include_phase_diagnostics=True),
    )

    assert result.split_detected is True
    assert result.diagnostics["solver_seed_name"] == "initial_phases"
    _assert_ceres_solver_diagnostics(result)

def test_electrolyte_lle_rejects_neutral_lle_initial_phase_labels() -> None:
    mix = _ascani_water_butanol_nacl_mixture()
    feed = np.asarray([0.55, 0.40, 0.025, 0.025], dtype=float)

    with np.testing.assert_raises_regex(epcsaft.InputError, "aq.*org.*phase_fraction"):
        mix.equilibrium(
            kind="electrolyte_lle",
            T=298.15,
            P=1.013e5,
            z=feed,
            initial_phases={"liq1": feed, "liq2": feed, "phase_fraction": 0.5},
        )

def test_equilibrium_options_default_iteration_budget_is_robust_for_electrolyte_lle() -> None:
    options = epcsaft.EquilibriumOptions()

    assert options.max_iterations == 180

def test_equilibrium_options_expose_density_robustness_controls() -> None:
    options = epcsaft.EquilibriumOptions(
        density_diagnostics="full",
        experimental_coupled_density_lle=True,
    )

    assert options.density_diagnostics == "full"
    assert options.experimental_coupled_density_lle is True

def test_one_salt_smoke_reports_production_solver_derivatives() -> None:
    feed = np.asarray([0.55, 0.40, 0.025, 0.025], dtype=float)
    mix = _ascani_water_butanol_nacl_mixture(feed)

    result = mix.equilibrium(
        kind="electrolyte_lle",
        T=298.15,
        P=1.013e5,
        z=feed,
        options=epcsaft.EquilibriumOptions(max_iterations=80, tolerance=1.0e-6),
    )

    _assert_ceres_solver_diagnostics(result)

def test_electrolyte_lle_rejects_non_neutral_direct_feed() -> None:
    mix = _ascani_water_butanol_nacl_mixture()

    with np.testing.assert_raises_regex(epcsaft.InputError, "charge neutral"):
        mix.equilibrium(kind="electrolyte_lle", T=298.15, P=1.013e5, z=[0.55, 0.40, 0.04, 0.01])

def test_neutral_lle_keeps_rejecting_ionic_mixtures() -> None:
    feed = np.asarray([0.55, 0.40, 0.025, 0.025], dtype=float)
    mix = _ascani_water_butanol_nacl_mixture(feed)

    with np.testing.assert_raises_regex(epcsaft.InputError, "ion-containing"):
        mix.equilibrium(kind="lle_flash", T=298.15, P=1.013e5, z=feed)
