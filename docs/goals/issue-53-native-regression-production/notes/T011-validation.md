# T011 Validation Receipt

## Result

Ran the required issue #53 validation set. The branch-specific native regression tests, API tests, benchmark, quick validation, docs build, ruff, black, and script syntax checks pass.

## Validation Commands

- `uv sync --no-install-project`
  - pass: resolved and checked the existing environment
- `uv run python scripts/build_epcsaft.py`
  - pass: CMake configure plus no-op Ninja build completed
- `uv run python run_pytest.py tests/native/test_native_regression_types.py tests/native/test_native_ceres_regression.py tests/native/test_native_regression_autodiff.py tests/native/test_native_reactive_regression.py -q`
  - pass: 12 passed
- `uv run python run_pytest.py tests/api/test_reactive_regression.py tests/api/test_runtime.py tests/api/test_parameter_schema.py -q`
  - pass: 54 passed
- `uv run python scripts/benchmark_native_regression.py --warmup 1 --repeat 3`
  - pass: all cases converged
- `uv run python scripts/validate_project.py quick`
  - pass: 23 passed in the quick pytest slice
- `uv run python scripts/validate_project.py docs`
  - pass: Sphinx build succeeded
- `uv run ruff check src tests scripts`
  - pass after mechanical cleanup of legacy script lint surfaced by the full gate
- `uv run black --check src tests scripts`
  - pass
- `uv run python -m py_compile ...`
  - pass for modified workflow/data scripts

## Native Benchmark Output

```text
case median_ms rows residuals params status backend
native_neutral_density_tiny 0.085 1 2 2 converged analytic_linear_native
native_binary_kij_tiny 0.076 1 3 2 converged analytic_linear_native
native_reactive_born_kij_tiny 0.087 1 4 2 converged analytic_linear_native
native_mea_pressure_speciation_35_row_surrogate 1.348 35 105 3 converged analytic_linear_native
```

## Lint Cleanup Note

The required full `ruff check src tests scripts` gate exposed older script lint in data/workflow helpers outside the native-regression files. I applied mechanical ruff fixes and narrowly retained file-level noqa only where paper filenames or legacy typing imports are intentionally preserved.
