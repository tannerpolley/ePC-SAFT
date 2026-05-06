# Chemical Equilibrium Test Coverage Handoff

## Purpose

Add focused tests that prove the new chemical-equilibrium and reactive-speciation path is coupled to ePC-SAFT thermodynamic activities/fugacity coefficients, not only solving ideal stoichiometric algebra.

This is a test-coverage task. Public APIs should not change unless a test exposes a real defect.

## Current State

The current implementation already has working coverage for:

- Native entrypoint exposure: `_core._solve_chemical_equilibrium_native`.
- Simple neutral reaction solve through `epcsaft.solve_reactive_speciation(...)`.
- Public facade routing through `mixture.equilibrium(kind="chemical_equilibrium", ...)`.
- Salt dissociation-style ionic speciation with mass, charge, and reaction residual closure.
- A MEA-like coupled ionic speciation smoke.
- Reactive stability handoff through `kind="reactive_stability"`.

Relevant files:

- `src/epcsaft/reactive_speciation.py`
- `src/epcsaft/epcsaft.py`
- `src/epcsaft/bindings.cpp`
- `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp`
- `tests/api/test_reactive_speciation.py`
- `tests/native/test_chemical_equilibrium_native.py`

## Highest-Priority Test To Add

Add a nonideal activity-coupling regression.

Goal: prove the equilibrium composition changes when the ePC-SAFT activity/fugacity model changes, even with the same reaction definition and `log_K`.

Recommended fixture:

- Neutral binary mixture, for example methanol/cyclohexane.
- Use existing helper shape from `tests/native/test_chemical_equilibrium_native.py`.
- Compute `log_K` from a target composition using ePC-SAFT fugacity-derived activities.
- Solve from a deliberately poor seed and assert the target composition is recovered.
- Then perturb a thermodynamic parameter such as `k_ij` while keeping the same `log_K`.
- Assert the solved composition shifts by a meaningful amount.

Acceptance criteria:

- First solve recovers the target composition within a tight tolerance.
- Perturbed-parameter solve still converges.
- Perturbed-parameter solution differs from the original by more than numerical noise.
- Diagnostics show the nonideal activity path, for example `activity_model == "epcsaft_neutral_fugacity_activity"`.

This catches accidental future regressions where the solver silently behaves like an ideal activity solver.

## Additional Focused Tests

### Trace-Species Robustness

Add a case where one species starts near `ReactiveSpeciationOptions.min_mole_fraction`.

Acceptance criteria:

- Solver succeeds or fails with a clear `SolutionError`.
- No negative mole fractions.
- No NaN or infinite values in `x`, residuals, activity coefficients, or diagnostics.
- If successful, residuals close within the requested tolerance.

### Inconsistent Chemistry Failure Contract

Add one intentionally invalid or inconsistent chemistry definition.

Good options:

- Unknown species in reaction stoichiometry.
- Missing material-balance total.
- Bad `initial_x` length.
- Impossible or underdetermined balance/reaction set that cannot converge.

Acceptance criteria:

- Input normalization failures raise `InputError`.
- Solver failures raise `SolutionError` with useful diagnostics.
- The error message points at the caller problem rather than exposing a native crash.

### Multi-Reaction Coupled Speciation

Add a small 4- or 5-species system with two reactions sharing species.

Acceptance criteria:

- All material balances close.
- All reaction residuals close.
- Charge residual closes when charged species are included.
- `result.to_dict()` is JSON-like with `json.dumps(..., allow_nan=False)`.

### Public Facade Parity

For the same small problem, compare:

- `epcsaft.solve_reactive_speciation(...)`
- `mixture.equilibrium(kind="chemical_equilibrium", ...)`

Acceptance criteria:

- Same composition within tolerance.
- Same mass, charge, and reaction residual behavior.
- Same native backend diagnostics.

## Suggested Implementation Order

1. Add the nonideal activity-coupling regression first.
2. Add one focused failure-contract test.
3. Add trace-species robustness if the failure-contract test does not already cover the edge.
4. Add facade parity if it can reuse the same fixture without much duplication.
5. Only add the multi-reaction test if it stays fast and clear.

Keep tests in:

- `tests/native/test_chemical_equilibrium_native.py` when the assertion is about native-backed chemistry behavior.
- `tests/api/test_reactive_speciation.py` when the assertion is about public API shape, result serialization, or error contracts.

## Validation Commands

Run the focused tests first:

```powershell
uv run python run_pytest.py tests/native/test_chemical_equilibrium_native.py tests/api/test_reactive_speciation.py -q
```

Then run the normal quick gate:

```powershell
uv run python scripts/validate_project.py quick
```

If any native files change, also run:

```powershell
uv run python scripts/build_epcsaft.py
```

Do not run full equilibrium or broad regression suites unless the implementation changes core solver behavior beyond tests.

## Non-Goals

- Do not add full MEA dataset validation.
- Do not add long optimization/fitting tests.
- Do not require plot generation.
- Do not change public API names.
- Do not make default validation depend on expensive chemical-equilibrium reproductions.

## Handoff Summary

The next agent should add small, contract-style tests proving that chemical equilibrium is thermodynamically coupled to ePC-SAFT activities/fugacity coefficients, robust near trace compositions, and clear on invalid chemistry definitions. The most valuable single test is the parameter-perturbation nonideal activity regression.
