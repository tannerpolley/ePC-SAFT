# T001 Judge: Phase-State Sensitivity Boundary

Date: 2026-05-15

## Decision

Approve an implicit-density-first Worker slice.

Reject a fixed-density-only implementation as the primary child path. The parent electrolyte LLE residual evaluates phase states by solving density at fixed temperature, pressure, phase, and composition before computing fugacity coefficients. A fixed-density Jacobian would differentiate a different state than the accepted residual uses and would not truthfully unblock parent T017.

## Approved Worker Slice

The next Worker should implement the smallest native phase-state sensitivity surface that follows the actual phase-state contract:

1. take temperature, pressure, phase, composition, and native parameter payload;
2. solve the selected density root using the existing density owner;
3. differentiate explicit pressure and fugacity expressions with respect to density and composition;
4. apply implicit density-root sensitivity at fixed pressure;
5. return a row-major `d ln(phi_i) / d x_j` payload plus density sensitivity diagnostics.

The first slice may stop at non-associating explicit-state coverage if active association requires second composition derivatives through solved site fractions. That stop is not completion; it must create a follow-up Judge or Worker gate before parent T017 resumes.

## Required Rejections

- No placeholder, stale, identity-only, or manually perturbed Jacobian.
- No derivative through density iteration history.
- No accepted parent electrolyte LLE Ceres claim until the residual evaluator consumes the returned surface.
- No broad EOS rewrite inside this first Worker.

## Worker Boundary

Allowed files:

- `src/epcsaft/native/epcsaft_core_internal.h`
- `src/epcsaft/native/epcsaft_ares.cpp`
- `src/epcsaft/native/epcsaft_Z.cpp`
- `src/epcsaft/native/epcsaft_fugcoef.cpp`
- `src/epcsaft/native/epcsaft_electrolyte.h`
- `src/epcsaft/native/epcsaft_state.cpp`
- `src/epcsaft/bindings.cpp`
- `src/epcsaft/epcsaft.py`
- `tests/native/cppad/**`
- `tests/native/equilibrium/**`
- `docs/goals/native-lle-reactive-production-solvers/subgoals/phase-state-sensitivities/notes/**`

Verify:

- `uv run python scripts/dev/build_epcsaft.py --profile full`
- `uv run python run_pytest.py tests/native/cppad tests/native/equilibrium -q`
- child GoalBuddy state checker

Stop if:

- the implementation needs second composition derivatives through active association site fractions before producing any useful real surface;
- the only available proof would be a manually perturbed derivative oracle;
- the surface cannot be wired to the same density/fugacity state used by `_evaluate_electrolyte_lle_residual_native`.
