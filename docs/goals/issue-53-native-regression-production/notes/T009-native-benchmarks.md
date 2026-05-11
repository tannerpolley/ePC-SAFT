# T009 Native Benchmark Receipt

## Result

Added a package-owned native regression benchmark harness:

- `src/epcsaft/benchmarks/native_regression.py`
- `scripts/benchmark_native_regression.py`
- `tests/workflows/test_benchmark_native_regression.py`

## Covered Cases

- `native_neutral_density_tiny`
- `native_binary_kij_tiny`
- `native_reactive_born_kij_tiny`
- `native_mea_pressure_speciation_35_row_surrogate`

The native reactive cases cover fixed-shape residuals, pressure/speciation/activity families, Born-radius parameters, binary-interaction parameters, and a 35-row public MEA-style surrogate without private downstream data.

## Benchmark Evidence

`uv run python scripts/benchmark_native_regression.py --warmup 1 --repeat 3`

```text
case median_ms rows residuals params status backend
native_neutral_density_tiny 0.036 1 2 2 converged analytic_linear_native
native_binary_kij_tiny 0.040 1 3 2 converged analytic_linear_native
native_reactive_born_kij_tiny 0.043 1 4 2 converged analytic_linear_native
native_mea_pressure_speciation_35_row_surrogate 0.670 35 105 3 converged analytic_linear_native
```

## Validation

- `uv run python run_pytest.py tests/workflows/test_benchmark_native_regression.py tests/native/test_native_reactive_regression.py -q`
  - 7 passed
- `uv run python run_pytest.py tests/workflows/test_benchmark_native_regression.py tests/native/test_native_reactive_regression.py tests/api/test_reactive_regression.py tests/api/test_runtime.py -q`
  - 55 passed
- focused `ruff check`
  - passed
- focused `black --check`
  - passed
