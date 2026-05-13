from __future__ import annotations

from tests.helpers.native_cases import _ionic_state


def _state():
    mix, _species, _pressure, density, temperature, composition = _ionic_state()
    return mix.state(T=temperature, x=composition, rho=density)


def test_derivative_coverage_matrix_enumerates_required_quantities() -> None:
    state = _state()

    rows = state.derivative_coverage_matrix()
    quantities = {row["quantity"] for row in rows}

    assert {
        "hard_chain",
        "dispersion",
        "association",
        "debye_huckel / ion",
        "born_direct",
        "born_ssmds_liquid",
        "relative_permittivity",
        "pressure",
        "fugacity",
        "activity",
        "chemical_potential",
        "k_hb_ij",
        "density_root",
    }.issubset(quantities)
    parameter_rows = [row for row in rows if row["derivative"] == "parameter"]
    assert {"relative_permittivity", "pressure", "fugacity", "activity", "chemical_potential"}.issubset(
        {row["quantity"] for row in parameter_rows}
    )
    assert "finite_difference" not in str(rows).lower()
    for row in rows:
        assert set(
            (
                "quantity",
                "derivative",
                "backend",
                "supported",
                "not_applicable",
                "classification",
                "backend_unavailable_reason",
                "source_equation_ids",
                "parameter_family",
                "future_owner",
            )
        ).issubset(row)


def test_derivative_coverage_matrix_uses_explicit_backend_labels() -> None:
    state = _state()

    rows = state.derivative_coverage_matrix()
    backend_labels = {row["backend"] for row in rows}

    assert "autodiff" not in backend_labels
    assert backend_labels.issubset(
        {
            "analytic",
            "cppad",
            "analytic_implicit",
            "cppad_implicit",
            "legacy_eigen_forward",
            "backend_unavailable",
        }
    )


def test_derivative_coverage_matrix_classifies_supported_blocked_and_out_of_scope_rows() -> None:
    state = _state()

    rows = state.derivative_coverage_matrix()
    classifications = {row["classification"] for row in rows}

    assert classifications.issubset(
        {
            "production_supported",
            "blocker",
            "blocker_requires_implicit_association_sensitivity",
            "out_of_scope",
        }
    )
    for row in rows:
        if row["not_applicable"]:
            assert row["classification"] == "out_of_scope"
        elif row["supported"]:
            assert row["classification"] == "production_supported"
        elif row["parameter_family"] == "k_hb_ij":
            assert row["classification"] == "blocker_requires_implicit_association_sensitivity"
        else:
            assert row["classification"] == "blocker"


def test_derivative_coverage_matrix_tracks_khbij_without_overclaiming() -> None:
    state = _state()

    rows = state.derivative_coverage_matrix()
    khb_rows = [row for row in rows if row["parameter_family"] == "k_hb_ij"]

    assert len(khb_rows) == 1
    row = khb_rows[0]
    assert row["quantity"] == "k_hb_ij"
    assert row["derivative"] == "parameter"
    assert row["backend"] == "backend_unavailable"
    assert row["supported"] is False
    assert row["classification"] == "blocker_requires_implicit_association_sensitivity"
    assert row["future_owner"] == "Task C"
    assert "implicit association site-fraction sensitivities" in row["backend_unavailable_reason"]
