# T002 Judge Architecture Decision

Date: 2026-05-11

## Result

Approved with constraints.

## Architecture

Implement a minimum viable native production slice for Issue #53:

- Native reactive mixed pressure/speciation regression in C++.
- Fixed-shape residuals.
- Penalty residuals for recoverable row failures.
- Ceres bounded least-squares.
- CppAD, analytic, or implicit derivatives for every residual family exposed by the supported production slice.

Do not label unsupported families as production. Unsupported derivative families must be rejected honestly at validation or capability gates.

## Ceres Option Decision

Reuse the repo-wide `EPCSAFT_ENABLE_CERES` and `EPCSAFT_USE_SYSTEM_CERES` options. Do not add a regression-specific public Ceres toggle.

T003 must broaden the option meaning from equilibrium-only to package-wide native Ceres support and report regression availability separately through runtime capabilities.

## CppAD Scope Decision

Production may ship before CppAD is threaded through every residual family in the entire package, but not before every residual family in the supported Issue #53 production slice has a non-finite-difference derivative path.

Finite differences may exist only behind explicit debug gates. No Python, SciPy, NumPy, or Ceres numeric-differentiation fallback is acceptable for production regression.

## Status Taxonomy Decision

Use a compatibility-mapped rollout with a new canonical public taxonomy:

- `converged`
- `max_iterations`
- `line_search_failed`
- `singular_jacobian`
- `all_rows_failed`
- `nonfinite_objective`
- `bounds_inconsistent`
- `invalid_input`

Compatibility mapping:

- old success -> `converged`
- exhausted budget -> `max_iterations`
- globalization failure -> `line_search_failed`
- singular linearization or `J.T @ J` -> `singular_jacobian`
- zero solved rows -> `all_rows_failed`
- nonfinite residual, objective, or Jacobian -> `nonfinite_objective`
- bad bounds -> `bounds_inconsistent`
- unsupported or invalid request -> `invalid_input`

Partial row failures belong in row diagnostics and penalty residual accounting, not as a top-level production status.

## Ordered Task Sequence

1. T003: build/dependency plumbing only. Broaden `EPCSAFT_ENABLE_CERES`; add `EPCSAFT_ENABLE_CPPAD`; expose honest capability reporting.
2. T004: native regression contracts/status/result serialization only.
3. T005: fixed-shape residual-family implementation for the minimum supported slice.
4. T006 and T007: coupled derivative-policy plus Ceres-backend tranche. No production claim before both are real.
5. T008: Python wrapper cutover only after the native tranche works.
6. T009: native benchmark fixtures and `scripts/benchmark_native_regression.py`.
7. T010: docs/capability cleanup after runtime truth is real.
8. T011 -> T012 -> T013 unchanged.

## Stop Conditions

- Stop after T003 if robust Windows/uv build plumbing for Ceres or CppAD cannot be made reproducible without a separate bootstrap decision.
- Stop after T005 if fixed-shape residual count/order cannot be guaranteed for recoverable row failures.
- Stop immediately if T006/T007 can only work by reintroducing Python/SciPy optimization, NumPy finite-difference Jacobians, Ceres numeric differentiation as production policy, or IPOPT as default/production.
- Stop for re-scope if the supported mixed pressure/speciation slice cannot get non-finite-difference derivatives without a broad whole-EOS templating rewrite.
- Do not mark any capability or docs path as production while `runtime.py` still reports `python_batched_native_solvers`, `native_hot_loop=False`, or Python finite-difference step control for the Issue #53 production route.

## Native Rebuild Coordination

- Main thread owns `_core` rebuild coordination.
- One native builder at a time.
- Workers must not run clean/repair or delete `_core` artifacts as routine validation.
- Rebuild after any CMake, `bindings.cpp`, or `native/regression/**` signature change.
- If `_core*.pyd` is locked, stop importing processes first; do not paper over with clean rebuild churn.
