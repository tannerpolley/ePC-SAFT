# T013 Remediation Receipt

## Result

Added analytic-sensitivity parameter movement to the native residual-record solve surface:

- `NativeRegressionResidualRecord` now accepts per-parameter `sensitivities`.
- Pybind parses `sensitivities` from Python dictionaries.
- Native solve builds and solves a small analytic normal-equation step, clamps parameters to bounds, re-evaluates residuals, and reports improved cost/residuals without production finite differences.
- Native benchmark records now include sensitivities for neutral, binary `k_ij`, reactive Born/`k_ij`, and 35-row MEA-style surrogate cases.
- Added a native test proving bounded parameter movement and objective reduction.

## Remaining Scope

This remediation improves the native solve surface but does not close the full T012 audit. Full issue #53 completion still requires native Ceres/CppAD thermodynamic row iteration for mixed pressure/speciation residuals, not Python-packed residual records.

## Validation

- `uv run python scripts/build_epcsaft.py --build-only --parallel 10`
  - pass
- `uv run python run_pytest.py tests/native/test_native_ceres_regression.py tests/workflows/test_benchmark_native_regression.py tests/api/test_reactive_regression.py -q`
  - 19 passed
- `uv run python scripts/benchmark_native_regression.py --warmup 1 --repeat 3`
  - all cases converged
- `uv run ruff check src tests scripts`
  - pass
- `uv run black --check src tests scripts`
  - pass

## Benchmark Output

```text
case median_ms rows residuals params status backend
native_neutral_density_tiny 0.055 1 2 2 converged analytic_linear_native
native_binary_kij_tiny 0.069 1 3 2 converged analytic_linear_native
native_reactive_born_kij_tiny 0.068 1 4 2 converged analytic_linear_native
native_mea_pressure_speciation_35_row_surrogate 1.144 35 105 3 converged analytic_linear_native
```
