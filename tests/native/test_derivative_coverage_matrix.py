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
        "density_root",
    }.issubset(quantities)
    assert "finite_difference" not in str(rows).lower()
    for row in rows:
        assert set(
            (
                "quantity",
                "derivative",
                "backend",
                "supported",
                "not_applicable",
                "backend_unavailable_reason",
                "source_equation_ids",
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
