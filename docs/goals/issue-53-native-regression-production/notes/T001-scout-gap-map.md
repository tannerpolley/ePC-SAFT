# T001 Scout Gap Map

Date: 2026-05-11

## Current Branch State

Branch `codex/issue-53-native-regression-production` has no implementation delta against `main` for the inspected Issue #53 files. Treat this as the current repo state versus GitHub issue #53.

## Requirement Gap Map

- Production regression must be native C++ with no Python optimizer loops or Python backend-unavailable Jacobians.
  - Status: gap.
  - Evidence: `fit_reactive_electrolyte_parameters()` in `src/epcsaft/reactive_regression.py` still owns the Gauss-Newton loop and calls `Backend_unavailable_jacobian()`.
- Ceres must be the production bounded least-squares backend.
  - Status: partial/gap.
  - Evidence: `CMakeLists.txt` has `EPCSAFT_ENABLE_CERES` and `EPCSAFT_USE_SYSTEM_CERES`, but the option is scoped to native equilibrium, default-off, and not wired to regression.
- CppAD must be the package-wide AD substrate.
  - Status: gap.
  - Evidence: no `EPCSAFT_ENABLE_CPPAD`, `EPCSAFT_USE_SYSTEM_CPPAD`, or CppAD include/build path exists. Current native AD uses Eigen autodiff in `src/epcsaft/native/epcsaft_regression.cpp`.
- Native source layout should move toward `src/epcsaft/native/regression/`.
  - Status: gap.
  - Evidence: regression is currently concentrated in `src/epcsaft/native/epcsaft_regression.cpp` with data structs in `src/epcsaft/native/epcsaft_electrolyte.h`.
- Python wrappers should validate, marshal, invoke native once, and serialize.
  - Status: partial/gap.
  - Evidence: neutral/generic wrappers already call native once; reactive regression still optimizes in Python.
- Public placeholder status `bounded_incomplete` must not appear.
  - Status: partial.
  - Evidence: current tests/docs guard against `bounded_incomplete`, but runtime status taxonomy is incomplete versus Issue #53 (`singular_jacobian`, `all_rows_failed`, `nonfinite_objective`, `bounds_inconsistent`, `invalid_input`).
- Backend unavailable must be debug-only, not production derivative policy.
  - Status: gap.
  - Evidence: reactive regression uses Python Backend unavailables in production; generic native debug/result paths still report backend-unavailable Jacobian backend.
- Mixed pressure/speciation reactive-electrolyte regression must work natively.
  - Status: partial/gap.
  - Evidence: row-level thermodynamic evaluations are native through chemical equilibrium and electrolyte bubble bindings, but the fitting optimizer remains Python-owned.
- Tests/docs/benchmarks must move to native-regression-specific artifacts.
  - Status: gap.
  - Evidence: no `tests/native/test_native_regression_types.py`, `test_native_ceres_regression.py`, `test_native_regression_autodiff.py`, `test_native_reactive_regression.py`, or `scripts/benchmark_native_regression.py` exists yet.

## Dependency And Build Risks

- Ceres on Windows/uv/Python 3.13: medium-high risk.
  - Existing FetchContent support is default-off and does not include package-local Windows dependency guardrails for Abseil/glog/gflags discovery.
- CppAD: medium risk.
  - Build integration should be manageable, but threading CppAD through residual evaluation and nested solves is the real implementation risk.
- Python 3.13/uv packaging: low-medium risk for Python, medium-high for optional native solver dependencies.
  - The likely failure point is native dependency discovery and reproducibility, not Python itself.

## Existing Reusable Native Surfaces

- Native regression structs already exist in `src/epcsaft/native/epcsaft_electrolyte.h`:
  - `PureNeutralRegressionDensityRecord`
  - `PureNeutralRegressionVLERecord`
  - `PureNeutralRegressionDebugResult`
  - `PureNeutralRegressionResult`
  - `GenericRegressionRecord`
  - `GenericRegressionDebugResult`
  - `GenericRegressionResult`
- Native row evaluators already exist for chemical equilibrium and electrolyte LLE/bubble paths:
  - `ChemicalResidualEvaluationNative`
  - `ElectrolyteLLEResidualEvaluationNative`
- Existing pybind entrypoints cover neutral/generic regression and row-level equilibrium solves, but not a native reactive-regression solve.

## Candidate Minimum Viable Slice

The smallest defensible first production slice is native Ceres solving for the existing reactive batch/context objective:

1. Add native reactive regression problem/result contracts.
2. Reuse existing native row evaluators rather than replacing thermodynamics.
3. Add one pybind entrypoint accepting prevalidated reactive rows, parameter specs, bounds, and solver options.
4. Return fixed-shape residuals, row diagnostics, and a complete status taxonomy.
5. Keep Python as API/data marshalling only.
6. Do not enable production backend-unavailable derivatives; reject unsupported derivative configurations rather than silently falling back.

## Candidate Worker Slices

- Build/dependency lane: CMake, Ceres/CppAD options, doctor/capabilities/tests.
- Native contracts lane: `src/epcsaft/native/regression/**`, status enums, residual schema, pybind serialization types.
- Native backend lane: Ceres bounded least squares and mixed pressure/speciation residual evaluation.
- Python surface lane: replace reactive Python optimizer with one-shot native binding, update capabilities and compatibility behavior.
- Tests/docs/benchmarks lane: native tests, benchmark script, docs, stale backend-unavailable production messaging cleanup.

## Ambiguities For T002 Judge

- Reuse repo-wide `EPCSAFT_ENABLE_CERES` or add regression-specific Ceres option?
- Is a native Ceres backend acceptable before CppAD is threaded through every residual family, or must production remain unavailable until all derivative paths are non-backend-unavailable?
- Should the new status taxonomy be a breaking public rename or a compatibility mapping around existing `failed_rows` behavior?

