# T021 Capabilities, Docs, And Benchmark Receipt

## Result

Added the public reporting and benchmark layer for the new native Ceres thermodynamic regression slice:

- Runtime capabilities now expose ``native_ceres_thermodynamic_regression`` with:
  - ``native_hot_loop=True``
  - ``python_objective_used=False``
  - ``finite_difference_used=False``
  - supported slice details
  - unsupported combinations returning ``backend_unavailable``
- ``bounded_mixed_pressure_speciation_regression`` now reports a partial native Ceres thermodynamic slice rather than the older residual-record-only status.
- Added ``scripts/benchmark_native_ceres_thermo_regression.py`` and the backing benchmark module.
- README and docs now explain the Ceres-enabled benchmark, finite-difference debug gate, and current limitations.

## Important Limitation

The default local build does not enable Ceres, so the benchmark reports:

- ``status=backend_unavailable``
- ``optimizer_backend=backend_unavailable``
- ``native_hot_loop=False``

The Ceres-enabled temp build from T019 proved the same native fit path can execute Ceres with objective decrease:

- ``optimizer_backend=ceres``
- ``derivative_backend=analytic_implicit``
- ``initial_cost=0.05174377188057001``
- ``final_cost=1.1771098705162825e-19``
- ``native_hot_loop=true; python_objective_used=false; finite_difference_used=false``

## Evidence

- Focused runtime/benchmark tests:
  - `uv run python run_pytest.py tests/api/test_runtime.py tests/workflows/test_benchmark_native_regression.py tests/workflows/test_benchmark_reactive_regression.py -q`
  - Result: `52 passed in 70.93s`
- Ceres thermodynamic benchmark script on default build:
  - `uv run python scripts/benchmark_native_ceres_thermo_regression.py --warmup 1 --repeat 3`
  - Result: `reactive_speciation_logk_implicit 0.341 backend_unavailable backend_unavailable implicit False 0.0517438 0.0517438`
- Lint/format:
  - `uv run ruff check src/epcsaft/benchmarks/native_ceres_thermo_regression.py scripts/benchmark_native_ceres_thermo_regression.py src/epcsaft/runtime.py tests/workflows/test_benchmark_native_regression.py tests/api/test_runtime.py`
  - `uv run black --check src/epcsaft/benchmarks/native_ceres_thermo_regression.py scripts/benchmark_native_ceres_thermo_regression.py src/epcsaft/runtime.py tests/workflows/test_benchmark_native_regression.py tests/api/test_runtime.py`

## Files

- `src/epcsaft/benchmarks/native_ceres_thermo_regression.py`
- `scripts/benchmark_native_ceres_thermo_regression.py`
- `src/epcsaft/runtime.py`
- `tests/api/test_runtime.py`
- `tests/workflows/test_benchmark_native_regression.py`
- `README.md`
- `docs/pages/diagnostics.rst`
- `docs/pages/parameter_regression.rst`
