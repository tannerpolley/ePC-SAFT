T041 receipt: CppAD chemical-equilibrium Jacobian tranche
=========================================================

Result
------

Completed the next native-derivative tranche after the `dP/dx` / `dp/drho`
substrate work.

What changed
------------

Added a real explicit `jacobian_backend="autodiff"` path for the native
chemical-equilibrium residual evaluator on the currently supported ideal
mole-fraction, log-species-amount variable model.

Implemented pieces:

- explicit CppAD Jacobian tape for ideal-mole-fraction chemical-equilibrium
  residuals in `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp`
- native thermo fit result backend inference so the Ceres thermo benchmark can
  distinguish `analytic_implicit` from `cppad_implicit`
- thermo benchmark request updated to exercise the explicit autodiff state
  Jacobian path

Why this matters
----------------

Before this tranche:

- the thermo Ceres benchmark was real, but it still reported
  `derivative=analytic_implicit`
- explicit `jacobian_backend="autodiff"` for chemical equilibrium raised
  immediately

After this tranche:

- explicit chemical-equilibrium autodiff Jacobians work for the supported ideal
  mole-fraction slice
- the thermo benchmark now reports `backend=ceres`,
  `derivative=cppad_implicit`, and `native_hot_loop=True`

Still not complete
------------------

This does **not** finish issue #53.

Remaining gaps still include:

- activity-coupled and concentration-coupled chemical-equilibrium Jacobians
- bubble-pressure differentiation in native thermo regression
- non-`logK` thermo parameter sensitivities such as Born/`f_solv`/`k_ij`
  through production native regression
- the older residual-record benchmark remains a separate contract benchmark and
  still honestly reports `analytic_linear_native`

Validation evidence
-------------------

Focused native/build validation passed:

- `uv run python scripts/build_epcsaft.py --build-only --parallel 10`
- `uv run python run_pytest.py tests/native/test_chemical_equilibrium_native.py tests/native/test_native_ceres_thermodynamic_regression.py tests/workflows/test_benchmark_native_regression.py tests/native/test_runtime_contracts.py tests/native/test_cppad_reactive_speciation_derivatives.py tests/native/test_cppad_eos_derivatives.py tests/native/test_cppad_bubble_derivatives.py -q`

Observed result:

- `61 passed in 3.92s`

Benchmark evidence:

- `uv run python scripts/benchmark_native_ceres_thermo_regression.py --warmup 1 --repeat 3`
  - `reactive_speciation_logk_implicit 1.518 converged ceres cppad_implicit True 0.0517438 1.17711e-19`
- `uv run python scripts/benchmark_native_regression.py --case native_reactive_born_ssmds_tiny --warmup 0 --repeat 1`
  - still `converged analytic_linear_native`, confirming that this older benchmark
    is a different residual-record contract path, not the thermo Ceres hot loop

Key files
---------

- `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp`
- `src/epcsaft/native/regression/thermo_regression.cpp`
- `src/epcsaft/benchmarks/native_ceres_thermo_regression.py`
- `tests/native/test_chemical_equilibrium_native.py`
- `tests/native/test_native_ceres_thermodynamic_regression.py`
- `tests/workflows/test_benchmark_native_regression.py`
