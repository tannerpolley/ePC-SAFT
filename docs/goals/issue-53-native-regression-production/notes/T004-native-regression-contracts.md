# T004 Native Regression Contracts Receipt

Date: 2026-05-11

## Result

Done.

## Changes

- Added native contract types under `src/epcsaft/native/regression/`.
- Added the canonical Issue #53 native regression status taxonomy:
  - `converged`
  - `max_iterations`
  - `line_search_failed`
  - `singular_jacobian`
  - `all_rows_failed`
  - `nonfinite_objective`
  - `bounds_inconsistent`
  - `invalid_input`
- Added native parameter-spec, residual-schema, row-diagnostic, and problem-contract structs.
- Exposed `_core._native_regression_contract_schema()` through pybind.
- Added `epcsaft.native_regression_contract_schema()` and `CANONICAL_NATIVE_REGRESSION_STATUSES`.
- Exported the contract helper and status tuple from the public package.
- Added native contract tests.
- Documented the contract and finite-difference production ban in `docs/pages/parameter_regression.rst`.

## Validation

- `uv run python scripts/build_epcsaft.py --build-only --parallel 10`
  - Passed after adding native contract files and pybind exposure.
- `uv run python run_pytest.py tests/native/test_native_regression_types.py tests/api/test_runtime.py -q`
  - 40 passed.
- `uv run ruff check src/epcsaft/native_regression.py src/epcsaft/__init__.py tests/native/test_native_regression_types.py tests/api/test_runtime.py`
  - Passed.
- `uv run black --check src/epcsaft/native_regression.py src/epcsaft/__init__.py tests/native/test_native_regression_types.py tests/api/test_runtime.py`
  - Passed.

## Notes

- This is contract-only. It does not implement Ceres solving or residual evaluation.
- The schema intentionally reports `production_finite_difference_allowed=False`.
