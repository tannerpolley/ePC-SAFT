# T006/T007 Derivative And Backend Gate Receipt

Date: 2026-05-11

## Result

Done for the supported native residual-record production contract slice.

## Changes

- Added `NativeRegressionFitOptions` and `NativeRegressionFitResult`.
- Added native solve surface for already-evaluated fixed-shape residual records.
- Enforced production derivative policy:
  - `analytic` and `cppad` are accepted derivative labels for the supported slice.
  - `finite_difference` returns `invalid_input`; there is no silent production fallback.
- Added native status handling for:
  - `converged`
  - `all_rows_failed`
  - `bounds_inconsistent`
  - `invalid_input`
- Added active-bound reporting and structured result serialization.
- Exposed `_core._solve_native_regression_residual_records(...)`.
- Added `epcsaft.solve_native_regression_residual_records(...)`.
- Added runtime capability entry for `native_residual_record_regression` with `production_finite_difference_allowed=False`.
- Added tests for:
  - finite-difference production rejection;
  - analytic derivative policy acceptance;
  - bounds-inconsistent mapping;
  - fixed-shape objective result payload;
  - all-rows-failed status;
  - pressure/speciation target-family slice with `k_ij` and `d_born` parameter kinds.

## Validation

- `uv run python scripts/build_epcsaft.py --build-only --parallel 10`
  - Passed.
- `uv run python run_pytest.py tests/native/test_native_regression_autodiff.py tests/native/test_native_ceres_regression.py tests/native/test_native_reactive_regression.py tests/native/test_native_regression_types.py tests/api/test_runtime.py -q`
  - 49 passed.
- `uv run ruff check src/epcsaft/native_regression.py src/epcsaft/__init__.py src/epcsaft/runtime.py tests/native/test_native_regression_autodiff.py tests/native/test_native_ceres_regression.py tests/native/test_native_reactive_regression.py tests/native/test_native_regression_types.py tests/api/test_runtime.py`
  - Passed.
- `uv run black --check src/epcsaft/native_regression.py src/epcsaft/__init__.py src/epcsaft/runtime.py tests/native/test_native_regression_autodiff.py tests/native/test_native_ceres_regression.py tests/native/test_native_reactive_regression.py tests/native/test_native_regression_types.py tests/api/test_runtime.py`
  - Passed.

## Scope And Honesty

This tranche is not a full thermodynamic nonlinear reactive-regression optimizer. It establishes the native production contract slice for fixed-shape residual records and rejects production finite differences. The existing `ReactiveElectrolyteRegressionContext` Python-batched fitter remains marked as compatibility, not Issue #53 native production.
