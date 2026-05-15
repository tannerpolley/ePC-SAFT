# T021 Worker: Electrolyte LLE Residual Jacobian

Date: 2026-05-15

## Result

Implemented the issue #116 transformed-variable electrolyte LLE residual Jacobian prerequisite.

The native residual evaluator now chains:

1. transformed beta/logit variables to organic formula composition;
2. formula material balance to aqueous formula composition;
3. formula compositions to explicit public species compositions;
4. pressure-based phase-state fugacity sensitivities from the completed child goal;
5. neutral and salt-pair phase-equilibrium residual rows;
6. explicit species material-balance rows.

## Changed Files

- `src/epcsaft/native/epcsaft_equilibrium.cpp`
- `tests/native/equilibrium/test_electrolyte_lle_residual_surface.py`
- `tests/native/equilibrium/test_electrolyte_lle_residual_jacobian.py`
- `tests/native/cppad/test_cppad_lle_derivatives.py`
- `docs/goals/native-lle-reactive-production-solvers/notes/T021-electrolyte-residual-jacobian.md`

## Evidence

- `_evaluate_electrolyte_lle_residual_native(...)` now returns a row-major transformed-variable Jacobian.
- The payload reports `jacobian_backend = cppad_implicit` and `jacobian_available = true`.
- Tests assert the returned gradient is `J.T @ residual`.
- The implementation uses the child phase-state sensitivity surface and does not introduce a manually perturbed derivative path.

## Validation

- `uv run python scripts/dev/build_epcsaft.py --profile full` passed.
- `uv run python run_pytest.py tests/native/equilibrium/test_electrolyte_lle_residual_surface.py tests/native/equilibrium/test_electrolyte_lle_residual_jacobian.py tests/native/cppad/test_phase_state_sensitivities.py -q` passed, 6 tests.
- `uv run python run_pytest.py tests/native/equilibrium tests/native/cppad -q` passed, 51 tests.
- `git diff --check` passed.

## Remaining Scope

The accepted electrolyte LLE production solve is still not a Ceres trust-region solve. The next Worker must wire the transformed residual/Jacobian into a Ceres cost function and accepted solve route before issue #116 can pass its Stage 6 and benchmark gates.
