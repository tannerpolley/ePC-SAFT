import os
import importlib

import pytest

from epcsaft.regression import _fit_pure_neutral_associating_native

baygi = importlib.import_module("analyses.2015_baygi.scripts._common")
MEA_MW = baygi.MEA_MW
TABLE2_MEA_PARAMETERS = baygi.TABLE2_MEA_PARAMETERS


def _baygi_mea_psat_pa(T: float) -> float:
    return baygi.baygi_mea_psat_pa(T)


def _baygi_mea_density_molar(T: float) -> float:
    return baygi.baygi_mea_density_molar(T)


def _baygi_mea_records():
    return [
        {
            "T": T,
            "P": _baygi_mea_psat_pa(T),
            "rho": _baygi_mea_density_molar(T),
        }
        for T in (303.15, 323.15, 343.15, 363.15)
    ]


def _perturbed(values):
    return {
        "m": values["m"] * 0.88,
        "s": values["s"] * 1.05,
        "e": values["e"] * 0.9,
        "e_assoc": values["e_assoc"] * 0.85,
        "vol_a": values["vol_a"] * 1.2,
    }


@pytest.mark.parametrize("assoc_scheme", ("2B", "3B", "4C"))
def test_baygi_table2_mea_parameters_score_better_than_perturbed_seed(assoc_scheme):
    table_values = TABLE2_MEA_PARAMETERS[assoc_scheme]
    table_score = _fit_pure_neutral_associating_native(
        _baygi_mea_records(),
        "MEA",
        assoc_scheme=assoc_scheme,
        fixed_parameters={"MW": MEA_MW},
        initial_guess=table_values,
        max_nfev=1,
    )
    perturbed_score = _fit_pure_neutral_associating_native(
        _baygi_mea_records(),
        "MEA",
        assoc_scheme=assoc_scheme,
        fixed_parameters={"MW": MEA_MW},
        initial_guess=_perturbed(table_values),
        max_nfev=1,
    )

    assert table_score.problem.fit_targets == ("m", "s", "e", "e_assoc", "vol_a")
    assert table_score.backend == "least_squares_native"
    assert table_score.fitted_values == pytest.approx(table_values, rel=0.0, abs=1.0e-12)
    assert table_score.residual_norm < perturbed_score.residual_norm


@pytest.mark.skipif(os.environ.get("EPCSAFT_RUN_MEA_TABLE2_REGRESSION") != "1", reason="opt-in MEA Table 2 optimizer")
@pytest.mark.parametrize("assoc_scheme", ("2B", "3B", "4C"))
def test_baygi_table2_mea_association_scheme_optimizer_does_not_worsen_seed(assoc_scheme):
    table_values = TABLE2_MEA_PARAMETERS[assoc_scheme]
    result = _fit_pure_neutral_associating_native(
        _baygi_mea_records(),
        "MEA",
        assoc_scheme=assoc_scheme,
        fixed_parameters={"MW": MEA_MW},
        initial_guess=_perturbed(table_values),
        bounds={
            "m": (0.5 * table_values["m"], 1.5 * table_values["m"]),
            "s": (0.8 * table_values["s"], 1.2 * table_values["s"]),
            "e": (0.5 * table_values["e"], 1.5 * table_values["e"]),
            "e_assoc": (0.5 * table_values["e_assoc"], 1.5 * table_values["e_assoc"]),
            "vol_a": (0.5 * table_values["vol_a"], 1.5 * table_values["vol_a"]),
        },
        multistart=1,
        max_nfev=8,
    )

    assert result.success, result.message
    assert result.backend == "least_squares_native"
    assert result.nfev > 1
    assert result.residual_norm <= result.metrics_by_term["initial_residual_norm"] + 1.0e-12
