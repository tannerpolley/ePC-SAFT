## T052 result

### Outcome

Aligned the public runtime reporting layer with the derivative behavior that is already implemented.

### What changed

1. `capabilities()["equilibrium"]["reactive_speciation"]`
   - `jacobian_auto_policy` now reports:
     - `cppad_supported_else_debug_fd_or_backend_unavailable`
   - supported auto/default standard states now list:
     - `ideal_mole_fraction`
     - `concentration`
     - `mole_fraction_activity`
   - added:
     - `jacobian_auto_ideal_without_cppad = "analytic"`

2. `capabilities()["regression"]["reactive_electrolyte_batch_context"]["native_ceres_thermodynamic_regression"]`
   - `derivative_backend` now reflects the live supported slice:
     - `cppad_implicit` on CppAD-enabled builds
   - the supported slice now includes:
     - reactive-speciation `logK` rows for `ideal_mole_fraction`, `concentration`, and supported `mole_fraction_activity`
     - the supported Born-SSM+DS parameter lane for `born_radius` / `d_born` and `f_solv` / `solvation_factor`
   - `blocked_parameter_kinds` now calls out the remaining generic binary/association interaction gaps instead of stale Born parameter blockers

3. docs/tests
   - updated runtime capability assertions
   - updated diagnostics docs to describe the new default Jacobian truth and the supported thermo slice

### Validation

- `uv run python run_pytest.py tests/api/test_runtime.py -q`
  - `37 passed in 1.12s`
- `uv run python scripts/benchmark_native_ceres_thermo_regression.py --case reactive_speciation_activity_ssmds_born_radius_implicit --warmup 1 --repeat 1`
  - `backend=ceres derivative=cppad_implicit native_hot_loop=True`

### Remaining truth

- public runtime `state.dadx()` auto/default still remains analytic
- reactive-electrolyte bubble differentiation still remains the larger derivative architecture blocker
