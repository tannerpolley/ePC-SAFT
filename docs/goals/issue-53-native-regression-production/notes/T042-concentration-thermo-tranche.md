T042 receipt: concentration thermo regression tranche
====================================================

Result
------

Extended the supported native thermodynamic regression slice from ideal
mole-fraction reactive-speciation `logK` rows to concentration-standard-state
reactive-speciation `logK` rows.

What changed
------------

1. Native chemical-equilibrium residuals now support
   concentration-standard-state `jacobian_backend="autodiff"` for supported
   nonassociating states by differentiating density closure through:

   - `dP/dx`
   - `dp/drho`

2. Public reactive-speciation auto-Jacobian routing for the concentration salt
   case now uses the native autodiff path instead of reporting
   `backend_unavailable`.

3. Native thermo regression no longer hard-rejects concentration standard
   states for the existing `reaction_equilibrium_constant` / `log_equilibrium_constant`
   parameter slice.

4. The thermo benchmark harness now includes a second case:

   - `reactive_speciation_concentration_logk_implicit`

Why this matters
----------------

This closes the next honest supported-slice gap after the original pressure
derivative blocker:

- supported ideal thermo Ceres path: `cppad_implicit`
- supported concentration thermo Ceres path: `cppad_implicit`

Still not complete
------------------

This still does **not** complete the full issue / thread goal.

Remaining unsupported derivative surfaces include:

- activity-coupled chemical-equilibrium Jacobians
- reactive-electrolyte bubble-pressure differentiation in thermo regression
- non-`logK` thermo parameter sensitivities (`d_born`, `f_solv`, `k_ij`, etc.)

Validation evidence
-------------------

Focused tests:

- `uv run python run_pytest.py tests/native/test_chemical_equilibrium_native.py tests/api/test_reactive_speciation.py tests/native/test_native_ceres_thermodynamic_regression.py tests/workflows/test_benchmark_native_regression.py tests/native/test_runtime_contracts.py tests/native/test_cppad_reactive_speciation_derivatives.py -q`
  - `79 passed, 1 skipped`

Benchmark evidence:

- `uv run python scripts/benchmark_native_ceres_thermo_regression.py --case reactive_speciation_logk_implicit --warmup 1 --repeat 1`
  - `backend=ceres derivative=cppad_implicit native_hot_loop=True`
- `uv run python scripts/benchmark_native_ceres_thermo_regression.py --case reactive_speciation_concentration_logk_implicit --warmup 1 --repeat 1`
  - `backend=ceres derivative=cppad_implicit native_hot_loop=True`

Key files
---------

- `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp`
- `src/epcsaft/native/regression/thermo_regression.cpp`
- `src/epcsaft/benchmarks/native_ceres_thermo_regression.py`
- `scripts/benchmark_native_ceres_thermo_regression.py`
- `tests/native/test_chemical_equilibrium_native.py`
- `tests/native/test_native_ceres_thermodynamic_regression.py`
- `tests/api/test_reactive_speciation.py`
- `tests/workflows/test_benchmark_native_regression.py`
