# T002 Worker: Phase-State Sensitivity Surface

Date: 2026-05-15

## Result

Implemented the first implicit-density phase-state sensitivity surface for pressure-based native states.

The new native surface:

- solves the selected density root at fixed temperature, pressure, phase, and composition;
- records explicit residual-Helmholtz dependence on density and composition with CppAD for non-associating states;
- derives pressure-density and fixed-density pressure-composition sensitivities;
- applies the density-root implicit chain rule at fixed pressure;
- returns row-major `d ln(phi_i) / d x_j` plus density and fixed-density diagnostic blocks.

## Changed Files

- `src/epcsaft/native/epcsaft_electrolyte.h`
- `src/epcsaft/native/epcsaft_core_internal.h`
- `src/epcsaft/native/epcsaft_ares.cpp`
- `src/epcsaft/bindings.cpp`
- `src/epcsaft/epcsaft.py`
- `tests/native/cppad/test_phase_state_sensitivities.py`

## Evidence

- `_native_phase_state_ln_fugacity_composition_sensitivity(...)` returns the native result contract used by future parent residual-Jacobian wiring.
- `ePCSAFTState.ln_fugacity_composition_derivative_result()` now reports support for pressure-based non-associating states and keeps density-based states unsupported.
- Tests assert the implicit pressure chain rule internally: `dP/dx|rho + dP/drho * drho/dx = 0`.
- Tests assert the returned total fugacity Jacobian equals the fixed-density Jacobian plus the density-chain contribution.
- The active-association case is now an explicit follow-up gate rather than a hidden completion claim.

## Validation

- `uv run python scripts/dev/build_epcsaft.py --profile full` passed.
- `uv run python run_pytest.py tests/native/cppad/test_phase_state_sensitivities.py -q` passed, 3 tests.
- `uv run python run_pytest.py tests/native/cppad tests/native/equilibrium -q` passed, 50 tests.
- `git diff --check` passed.

## Remaining Gate

Parent T017 is not unblocked yet. The issue #116 electrolyte fixtures can include associating solvents, so a Judge task must decide the next association/SSM-DS derivative slice before the child goal can report `parent_t017_unblocked: true`.
