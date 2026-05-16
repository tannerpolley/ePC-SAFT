from __future__ import annotations

import numpy as np
import pytest

import epcsaft
from epcsaft import ePCSAFTMixture


def _ascani_water_butanol_nacl_mixture(feed=None) -> ePCSAFTMixture:
    species = ["H2O", "Butanol", "Na+", "Cl-"]
    if feed is None:
        feed = np.asarray([0.55, 0.40, 0.025, 0.025], dtype=float)
    return ePCSAFTMixture.from_dataset("2022_Ascani", species, feed, 298.15)


def _assert_electrolyte_lle_route_pending(excinfo: pytest.ExceptionInfo[epcsaft.InputError]) -> None:
    message = str(excinfo.value)
    assert "electrolyte_lle requires a native Ipopt equilibrium NLP route" in message
    assert "No package-owned alternate LLE solver is available" in message


def test_one_salt_electrolyte_lle_direct_feed_requires_native_ipopt_route() -> None:
    feed = np.asarray([0.55, 0.40, 0.025, 0.025], dtype=float)
    mix = _ascani_water_butanol_nacl_mixture(feed)

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="electrolyte_lle",
            T=298.15,
            P=1.013e5,
            z=feed,
            options=epcsaft.EquilibriumOptions(include_phase_diagnostics=True),
        )

    _assert_electrolyte_lle_route_pending(excinfo)


def test_electrolyte_lle_direct_feed_requested_ipopt_requires_native_ipopt_route() -> None:
    feed = np.asarray([0.55, 0.40, 0.025, 0.025], dtype=float)
    mix = _ascani_water_butanol_nacl_mixture(feed)

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="electrolyte_lle",
            T=298.15,
            P=1.013e5,
            z=feed,
            options=epcsaft.EquilibriumOptions(
                solver_backend="ipopt",
                max_iterations=80,
                tolerance=1.0e-8,
            ),
        )

    _assert_electrolyte_lle_route_pending(excinfo)


def test_electrolyte_lle_molality_feed_requires_native_ipopt_route() -> None:
    mix = _ascani_water_butanol_nacl_mixture()

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="electrolyte_lle",
            T=298.15,
            P=1.013e5,
            solvent_feed={"H2O": 0.58, "Butanol": 0.42},
            salt_molality={"NaCl": 1.0},
            options=epcsaft.EquilibriumOptions(include_phase_diagnostics=True),
        )

    _assert_electrolyte_lle_route_pending(excinfo)


def test_electrolyte_lle_validates_strict_aq_org_initial_phases_before_route_gate() -> None:
    mix = _ascani_water_butanol_nacl_mixture()
    aq = np.asarray([0.798324680201737, 0.016320352824141723, 0.09267748348706063, 0.09267748348706063])
    org = np.asarray([0.37006036048879404, 0.6214918588210971, 0.004223890345054407, 0.004223890345054407])
    beta_org = 0.613766575013417
    feed = (1.0 - beta_org) * aq + beta_org * org

    with pytest.raises(epcsaft.InputError) as excinfo:
        mix.equilibrium(
            kind="electrolyte_lle",
            T=298.15,
            P=1.013e5,
            z=feed,
            initial_phases={"aq": aq, "org": org, "phase_fraction": beta_org},
            options=epcsaft.EquilibriumOptions(include_phase_diagnostics=True),
        )

    _assert_electrolyte_lle_route_pending(excinfo)


def test_electrolyte_lle_rejects_neutral_lle_initial_phase_labels() -> None:
    mix = _ascani_water_butanol_nacl_mixture()
    feed = np.asarray([0.55, 0.40, 0.025, 0.025], dtype=float)

    with pytest.raises(epcsaft.InputError, match=r"aq.*org.*phase_fraction"):
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


def test_electrolyte_lle_rejects_non_neutral_direct_feed() -> None:
    mix = _ascani_water_butanol_nacl_mixture()

    with pytest.raises(epcsaft.InputError, match="charge neutral"):
        mix.equilibrium(kind="electrolyte_lle", T=298.15, P=1.013e5, z=[0.55, 0.40, 0.04, 0.01])


def test_neutral_lle_keeps_rejecting_ionic_mixtures() -> None:
    feed = np.asarray([0.55, 0.40, 0.025, 0.025], dtype=float)
    mix = _ascani_water_butanol_nacl_mixture(feed)

    with pytest.raises(epcsaft.InputError, match="ion-containing"):
        mix.equilibrium(kind="lle_flash", T=298.15, P=1.013e5, z=feed)
