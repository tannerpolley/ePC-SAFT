# T005 Native Residual Evaluation Receipt

Date: 2026-05-11

## Result

Done.

## Changes

- Added native residual-record and residual-evaluation structs.
- Added native fixed-shape residual packing for already-evaluated target records.
- Added native penalty residual handling for recoverable row failures.
- Added row diagnostics serialization with the advertised row diagnostic fields.
- Exposed `_core._evaluate_native_regression_residual_records(...)`.
- Added `epcsaft.evaluate_native_regression_residual_records(...)`.
- Added tests proving:
  - residual order remains fixed;
  - schema carries residual indices and target families;
  - recoverable row failures produce penalty residuals;
  - row diagnostics carry failure status/message;
  - invalid scales are rejected.

## Validation

- `uv run python scripts/build_epcsaft.py --build-only --parallel 10`
  - Passed.
- `uv run python run_pytest.py tests/native/test_native_regression_types.py tests/api/test_runtime.py -q`
  - 42 passed.
- `uv run ruff check src/epcsaft/native_regression.py src/epcsaft/__init__.py tests/native/test_native_regression_types.py tests/api/test_runtime.py`
  - Passed.
- `uv run black --check src/epcsaft/native_regression.py src/epcsaft/__init__.py tests/native/test_native_regression_types.py tests/api/test_runtime.py`
  - Passed.

## Notes

- This tranche implements residual packing and diagnostics, not Ceres solving.
- Production derivatives are still blocked until the coupled T006/T007 tranche.
