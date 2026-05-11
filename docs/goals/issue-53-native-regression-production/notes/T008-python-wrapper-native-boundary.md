# T008 Python Wrapper Native Boundary Receipt

## Result

Implemented a public reactive-regression backend switch:

- `backend="native"` is now the default for `fit_reactive_electrolyte_parameters(...)`.
- The native path evaluates the structured row objective once, packs fixed-shape residual records, and calls `solve_native_regression_residual_records(...)`.
- The native path does not run the Python Gauss-Newton loop or Python finite-difference Jacobian.
- The old bounded Python optimizer remains available only as `backend="python_compat"` and reports `production_ready = false`, `python_optimizer = true`, and `finite_difference_jacobian = true`.

## Scope Note

This task removes Python optimization from the public default path, but it does not complete the full Issue #53 Ceres nonlinear thermodynamic iteration. Runtime capabilities remain honest: the mixed pressure/speciation path is a native boundary contract slice, not the finished full Ceres production backend.

## Evidence

- `src/epcsaft/reactive_regression.py`
- `src/epcsaft/runtime.py`
- `tests/api/test_reactive_regression.py`
- `tests/api/test_runtime.py`
- `docs/pages/parameter_regression.rst`

## Validation

- `uv run python run_pytest.py tests/api/test_reactive_regression.py -q`
  - 11 passed
- `uv run python scripts/build_epcsaft.py --build-only --parallel 10`
  - native rebuild passed
- `uv run python run_pytest.py tests/native/test_native_regression_autodiff.py tests/native/test_native_ceres_regression.py tests/native/test_native_reactive_regression.py tests/native/test_native_regression_types.py tests/api/test_runtime.py tests/api/test_reactive_regression.py -q`
  - 60 passed
