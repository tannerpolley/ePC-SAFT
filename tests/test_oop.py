import numpy as np
import pytest

from pcsaft import (
    PCSAFTMixture,
    pcsaft_den,
    pcsaft_ares,
    pcsaft_gres,
    pcsaft_lnfugcoef,
    pcsaft_multiphase_lle,
    pcsaft_p,
    pcsaft_Hvap,
)
from pcsaft.parameters import get_prop_dict


def test_dataset_state_matches_legacy_density_and_lnfugcoef():
    t = 298.15
    p = 1.0e5
    species = ["Li+", "Br-", "Ethanol"]
    x = np.asarray([0.03, 0.03, 0.94], dtype=float)

    mixture = PCSAFTMixture.from_dataset("2020_Bulow", species)
    state = mixture.state(T=t, x=x, P=p, phase="liq")

    legacy_params = get_prop_dict("2020_Bulow", species, x, t, user_options={})
    rho_legacy = pcsaft_den(t, p, x, legacy_params, phase="liq")

    assert state.density() == pytest.approx(rho_legacy, rel=0.0, abs=1e-12)
    assert state.pressure() == pytest.approx(p, rel=0.0, abs=1e-12)
    assert np.allclose(state.lnfugcoef(), np.asarray(pcsaft_lnfugcoef(t, rho_legacy, x, legacy_params)), atol=1e-12)


def test_raw_state_matches_legacy_energy_properties():
    t = 325.0
    x = np.asarray([1.0], dtype=float)
    params = {"m": np.asarray([2.8149]), "s": np.asarray([3.7169]), "e": np.asarray([285.69])}

    mixture = PCSAFTMixture.from_params(["Toluene"], params)
    state = mixture.state(T=t, x=x, rho=9135.590853014008, phase="liq")

    assert state.ares() == pytest.approx(pcsaft_ares(t, 9135.590853014008, x, params), rel=0.0, abs=1e-12)
    assert state.gres() == pytest.approx(pcsaft_gres(t, 9135.590853014008, x, params), rel=0.0, abs=1e-12)
    assert state.pressure() == pytest.approx(pcsaft_p(t, 9135.590853014008, x, params), rel=0.0, abs=1e-12)


def test_legacy_wrappers_emit_warning_and_preserve_shapes():
    t = 325.0
    x = np.asarray([1.0], dtype=float)
    params = {"m": np.asarray([2.8149]), "s": np.asarray([3.7169]), "e": np.asarray([285.69])}

    with pytest.warns(DeprecationWarning):
        rho = pcsaft_den(t, 101325.0, x, params, phase="liq")
    assert rho > 0.0

    with pytest.warns(DeprecationWarning):
        result = pcsaft_Hvap(380.0, x, params)
    assert isinstance(result, list)
    assert len(result) == 2


def test_multiphase_lle_state_result_round_trips_to_legacy_dict():
    t = 298.15
    p = 1.0e5
    species = ["H2O-2B-Li", "Na+", "Cl-"]
    x = np.asarray([1.0 / 0.01801528, 1e-4, 1e-4], dtype=float)
    x = x / np.sum(x)

    mixture = PCSAFTMixture.from_dataset("2020_Bulow", species)
    state = mixture.state(T=t, x=x, P=p, phase="liq")
    result = state.multiphase_lle(options={"tpdf_global_trials": 20, "tpdf_local_trials": 10, "max_nfev": 25})

    legacy = result.to_legacy()
    assert legacy["n_phases"] in (1, 2)
    assert "phases" in legacy
    assert "e_matrix" in legacy

    with pytest.warns(DeprecationWarning):
        wrapper = pcsaft_multiphase_lle(t, p, x, get_prop_dict("2020_Bulow", species, x, t), species, options={"tpdf_global_trials": 20, "tpdf_local_trials": 10, "max_nfev": 25})
    assert wrapper["n_phases"] == legacy["n_phases"]
