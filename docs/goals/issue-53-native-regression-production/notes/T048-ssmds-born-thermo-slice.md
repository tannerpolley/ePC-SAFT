## T048 receipt

Implemented the first honest non-`logK` native thermo regression tranche for the supported reactive-speciation path.

### Scope completed

- Extended native `ln(phi)` composition and density derivatives to support the nonassociating ionic Born SSM+DS path, including solvent-reference Born epsilon mode on the supported `dielc_rule in {0,1}` slice.
- Extended native component-activity `d log(gamma) / d log(n)` support to that same SSM+DS slice.
- Added fixed-`rho` native `ln(phi)` parameter derivatives for:
  - `born_radius`
  - `f_solv`
- Wired those parameter derivatives into native thermo regression `R_theta` for activity-standard-state `reactive_speciation` rows.
- Corrected `fill_implicit_speciation_jacobian(...)` to evaluate the inner residual Jacobian at the converged speciation state instead of the row seed.
- Added one benchmark case proving the new Born parameter slice runs through `ceres` with `cppad_implicit`.

### Files touched

- `src/epcsaft/native/epcsaft_fugcoef.cpp`
- `src/epcsaft/native/epcsaft_activity.cpp`
- `src/epcsaft/native/epcsaft_core_internal.h`
- `src/epcsaft/native/epcsaft_electrolyte.h`
- `src/epcsaft/native/regression/thermo_regression.cpp`
- `src/epcsaft/benchmarks/native_ceres_thermo_regression.py`
- `tests/native/test_chemical_equilibrium_native.py`
- `tests/native/test_native_ceres_thermodynamic_regression.py`

### Validation

- `uv run python scripts/build_epcsaft.py --build-only --parallel 10`
  - pass
- `uv run python run_pytest.py tests/native/test_chemical_equilibrium_native.py tests/native/test_native_ceres_thermodynamic_regression.py -q`
  - `29 passed in 2.08s`
- `uv run python run_pytest.py tests/api/test_reactive_speciation.py tests/native/test_runtime_contracts.py tests/workflows/test_benchmark_native_regression.py tests/native/test_cppad_reactive_speciation_derivatives.py tests/native/test_chemical_equilibrium_native.py tests/native/test_native_ceres_thermodynamic_regression.py -q`
  - `87 passed, 1 skipped in 6.51s`
- `uv run python scripts/benchmark_native_ceres_thermo_regression.py --case reactive_speciation_activity_ssmds_born_radius_implicit --warmup 1 --repeat 1`
  - `backend=ceres derivative=cppad_implicit native_hot_loop=True`
  - cost `0.0202084 -> 0.0201522`

### Remaining truth after T048

- Native thermo regression still does **not** support generic non-`logK` dispersion/binary parameters like `k_ij` on the reactive-speciation thermo path.
- Reactive-electrolyte bubble differentiation remains a separate downstream blocker.
- This tranche closed the Born `d_born` / `f_solv` slice requested for the SSM+DS activity-standard-state path; it did not close the full issue.
