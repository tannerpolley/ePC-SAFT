# T004 Worker: Association Sensitivity Extension

Date: 2026-05-15

## Result

Extended the phase-state sensitivity surface through active association solved-state sensitivities.

The extension records association outputs and residual equations with variables `[rho, x..., XA...]`, solves the association residual sensitivity system for site-fraction responses, and adds association `dmu/dv` and `dzraw/dv` terms into the existing pressure-based fugacity chain rule.

## Changed Files

- `src/epcsaft/native/epcsaft_ares.cpp`
- `tests/native/cppad/test_phase_state_sensitivities.py`
- `docs/goals/native-lle-reactive-production-solvers/subgoals/phase-state-sensitivities/notes/T004-association-sensitivity-extension.md`

## Evidence

- The active-association ionic helper now returns `supported = true` from `_native_phase_state_ln_fugacity_composition_sensitivity(...)`.
- The surface still differentiates pressure-based phase states, not fixed-density-only diagnostics.
- Tests assert the pressure closure chain rule and fugacity chain rule for both non-associating and active-association cases.
- Existing association implicit derivative contract tests still pass.

## Validation

- `uv run python scripts/dev/build_epcsaft.py --profile full` passed.
- `uv run python run_pytest.py tests/native/cppad/test_phase_state_sensitivities.py tests/native/contracts/test_association_implicit_derivative_contract.py -q` passed, 6 tests.
- `uv run python run_pytest.py tests/native/cppad tests/native/equilibrium -q` passed, 50 tests.
- `git diff --check` passed.

## Remaining Scope

The child surface now covers the pressure-based fugacity sensitivity foundation for direct-Born associating electrolyte states. Parent T017 still must wire this surface into the electrolyte LLE transformed-variable residual Jacobian and Ceres cost function.
