# T022 Worker: Ceres Accepted Electrolyte LLE Solve

## Outcome

Implemented a native Ceres trust-region accepted solve for issue #116 electrolyte LLE candidates in `src/epcsaft/native/epcsaft_equilibrium.cpp`.

The accepted production route now:

1. evaluates transformed electrolyte LLE residuals on charge-constrained formula variables;
2. routes accepted candidates through `ceres::Problem`;
3. reports `solver_backend = ceres`;
4. reports `solver_method = ceres_trust_region_residual_solve`;
5. reports Ceres trust-region, linear-solver, termination, iteration, and cost diagnostics;
6. preserves the T021 public residual/Jacobian payload with `jacobian_backend = cppad_implicit`.

## Important Follow-Up

The multi-salt distributed-ion diagnostic case exposed that the public T021 analytic transformed-variable Jacobian is not yet proven as the robust Ceres step Jacobian for all salt-basis cases. The Ceres cost currently uses a local residual slope for accepted-solve robustness while the public residual evaluator still returns the T021 analytic payload.

Do not close issue #116 until the next task either:

- repairs the analytic transformed-variable Jacobian for the multi-salt distributed-ion case and wires it directly into Ceres; or
- records a source-backed decision that the local residual slope is the intended production Jacobian surface and updates diagnostics/docs/tests honestly.

## Validation

- `uv run python scripts/dev/build_epcsaft.py --profile full` passed.
- `uv run python run_pytest.py tests/native/equilibrium/test_electrolyte_lle_ceres_solver.py tests/native/equilibrium/test_electrolyte_lle_residual_jacobian.py -q` passed, 2 tests.
- `uv run python run_pytest.py tests/native/equilibrium tests/native/cppad -q` passed, 52 tests.
- `git diff --check` passed.
