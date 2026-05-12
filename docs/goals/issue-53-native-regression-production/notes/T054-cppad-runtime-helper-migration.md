## T054 receipt

Objective completed:
- migrated the remaining supported public/runtime derivative helpers from older AutoDual-active paths onto native CppAD in CppAD-enabled builds, without changing the public `dadx()` auto/default contract and without claiming bubble differentiation is solved

Implemented slices:
- `src/epcsaft/native/epcsaft_Z.cpp`
  - supported pressure-composition helper now reports `cppad_composition`
  - supported pressure-density helper now reports `cppad_density`
- `src/epcsaft/native/epcsaft_regression.cpp`
  - fused-state `dpdrho` / `dlnfug_drho` helper moved to CppAD
- `src/epcsaft/native/epcsaft_parameter_setup.cpp`
  - dielectric autodiff derivative helper now uses a CppAD Jacobian
  - supported explicit contribution autodiff path now uses CppAD for HC, DISP, ION, and Born-model-1 contributions in CppAD-enabled builds
- `src/epcsaft/native/epcsaft_fugcoef.cpp`
  - supported `lnphi` helper backends now report `cppad_composition`, `cppad_density`, `cppad_parameter`
- `src/epcsaft/native/epcsaft_activity.cpp`
  - activity helper backends now report `cppad_component_activity_log_amounts` and `cppad_component_activity_parameter`
- `src/epcsaft/native/contributions/epcsaft_contrib_born.cpp`
  - added CppAD reference-solvent dielectric helper for Born-mode-1 composition derivatives

Validation:
- `uv run python scripts/build_epcsaft.py --build-only --parallel 10`
- `uv run python run_pytest.py tests/native/test_cppad_eos_derivatives.py tests/native/test_runtime_contracts.py tests/api/test_runtime.py -q`
  - `65 passed in 2.30s`
- broader runtime slice:
  - `uv run python run_pytest.py tests/native/test_runtime_contracts.py tests/api/test_runtime.py tests/native/test_chemical_equilibrium_native.py tests/workflows/test_benchmark_native_regression.py -q`
  - `90 passed in 7.23s`
- benchmark sanity:
  - `uv run python scripts/benchmark_native_regression.py`
  - converged on all listed cases, backend unchanged as `analytic_linear_native` for that residual-record contract benchmark

Not completed by T054:
- reactive-speciation native thermo regression still does not support generic `k_ij` / `l_ij` / `k_hb_ij`
- reactive-electrolyte bubble derivatives are still not differentiated
