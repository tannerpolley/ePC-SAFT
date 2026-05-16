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
        "relative_permittivity",
        "pressure",
    }.issubset(quantities)
    parameter_rows = [row for row in rows if row["derivative"] == "parameter"]
    assert {"relative_permittivity"}.issubset({row["quantity"] for row in parameter_rows})
    assert "numerical_derivative" not in str(rows).lower()
    for row in rows:
        assert set(
            (
                "quantity",
                "derivative",
                "backend",
                "supported",
                "not_applicable",
                "classification",
                "source_equation_ids",
                "parameter_family",
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
            "out_of_scope",
        }
    )
    assert "unsupported" not in backend_labels


def test_derivative_coverage_matrix_classifies_supported_and_out_of_scope_rows() -> None:
    state = _state()

    rows = state.derivative_coverage_matrix()
    classifications = {row["classification"] for row in rows}

    assert classifications.issubset(
        {
            "production_supported",
            "out_of_scope",
        }
    )
    for row in rows:
        if row["not_applicable"]:
            assert row["classification"] == "out_of_scope"
        else:
            assert row["supported"] is True
            assert row["classification"] == "production_supported"
