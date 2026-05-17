from __future__ import annotations

import epcsaft


def test_equilibrium_capabilities_expose_derivative_policy() -> None:
    policy = epcsaft.capabilities()["equilibrium"]["derivative_policy"]

    assert policy["unsupported_derivative_behavior"] == "raise"
    assert policy["auto_policy"] == "analytic_or_cppad_or_implicit_else_raise"
    assert "numerical" + "_derivative" not in {str(item).lower() for item in policy["accepted_derivative_backends"]}
    assert {
        "analytic",
        "cppad",
        "analytic_implicit",
        "cppad_implicit",
    }.issubset(set(policy["accepted_derivative_backends"]))

