T029 receipt: Ceres-enabled release-path proof
==============================================

Decision
--------

The first T029 subproblem was build proof, not additional thermodynamic math:
the native thermodynamic regression slice could not be trusted as production
until this checkout could actually build and run Ceres/CppAD.

Finding
-------

The first clean Ceres/CppAD build failed at the final `_core` link on Windows
MinGW. Bundled Ceres auto-enabled LAPACK and its static library referenced
`dpotrf_`, `dpotrs_`, `dgeqrf_`, `dormqr_`, and related symbols that were not
linked into the extension.

Fix
---

Bundled FetchContent Ceres now forces `LAPACK=OFF`. This keeps the small native
least-squares slice on Ceres' Eigen-backed path and avoids a fragile external
BLAS/LAPACK dependency in the local package build.

Validation
----------

- `uv run python scripts/build_epcsaft.py --clean --enable-ceres --enable-cppad --parallel 10`
  - pass
  - configure: 37.51s
  - build: 172.85s
  - total: 211.29s
- `uv run python scripts/benchmark_native_ceres_thermo_regression.py --warmup 1 --repeat 3`
  - pass
  - `reactive_speciation_logk_implicit`
  - median: 1.399 ms
  - status: `converged`
  - backend: `ceres`
  - derivative: `analytic_implicit`
  - native_hot_loop: `True`
  - initial_cost: 0.0517438
  - final_cost: 1.17711e-19
- `uv run python scripts/doctor.py`
  - pass
  - `native_dependency_ceres: enabled=ON, found=ON, available=True`
  - `native_dependency_cppad: enabled=ON, found=ON, available=True`
- `uv run python run_pytest.py tests/workflows/test_build_epcsaft_script.py tests/native/test_native_ceres_thermodynamic_regression.py tests/api/test_runtime.py -q`
  - pass: 47 passed in 1.33s

Remaining work
--------------

This removes the local Ceres proof blocker. It does not implement the remaining
issue #53 derivative gaps:

- reactive-electrolyte bubble-pressure residual/Jacobian support;
- Born-SSM+DS `d_born` / `f_solv` sensitivities through activity/fugacity and
  density closure.
