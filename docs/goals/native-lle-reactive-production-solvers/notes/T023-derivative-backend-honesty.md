# T023 Worker: Accepted-Solve Derivative Backend Honesty

## Outcome

Resolved the immediate issue #116 diagnostic gap created by T022 by making the accepted electrolyte LLE Ceres solve report the derivative backend it actually uses.

The accepted Ceres solve now reports:

- `jacobian_backend = local_residual_slope`
- `derivative_backend = local_residual_slope`
- `residual_surface_jacobian_backend = cppad_implicit`
- `residual_surface_derivative_backend = cppad_implicit`

This keeps the public residual-surface T021 CppAD payload visible while avoiding a false claim that the accepted Ceres solve is directly using the T021 analytic transformed-variable Jacobian.

## Completion Boundary

This does not close issue #116 by itself. The accepted native Ceres route is installed and validated, but the issue's stronger production-Jacobian requirement still needs final audit acceptance or a follow-up analytic multi-salt repair.

## Validation

- `uv run python scripts/dev/build_epcsaft.py --profile full` passed.
- `uv run python run_pytest.py tests/native/equilibrium/test_electrolyte_lle_ceres_solver.py tests/native/equilibrium/test_electrolyte_lle_residual_jacobian.py tests/native/cppad/test_phase_state_sensitivities.py -q` passed, 5 tests.
- `uv run python run_pytest.py tests/native/equilibrium tests/native/cppad -q` passed, 52 tests.
- `git diff --check` passed.
