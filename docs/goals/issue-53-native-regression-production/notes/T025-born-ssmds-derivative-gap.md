# T025 Born-SSM+DS Derivative Gap

## Decision

T025 is stopped on the task stop condition. The package should not claim native
Ceres production sensitivities for Born-SSM+DS parameters until those
sensitivities are analytic, CppAD-backed, or implicit over a scalar-templated
native residual path. A native parameter-perturbation shortcut would still be a
backend-unavailable derivative and is not acceptable for Issue #53 production.

## Scope Correction

Reactive electrolyte regression coverage should focus on the Born-SSM+DS
parameters that belong to the Born model:

- `d_born` / `born_radius`
- `f_solv` / `solvation_factor`

Generic binary interaction parameters belong to the binary regression lane:

- `k_ij`
- `l_ij`
- `k_hb_ij`

Those should be regressed through `fit_binary_pair(...)` against direct binary
composition data, not through a reactive-speciation Born sensitivity test.

## Evidence

- `src/epcsaft/native/regression/thermo_regression.cpp` now accepts
  `solvation_factor` / `f_solv` in the native thermodynamic parameter adapter,
  alongside `born_radius` / `d_born`.
- `src/epcsaft/native/regression/regression_types.cpp` now exposes
  `solvation_factor` in the native regression contract.
- `src/epcsaft/runtime.py` reports the missing Born-SSM+DS derivative call graph
  and points generic binary `k_ij`/`l_ij`/`k_hb_ij` work to `fit_binary_pair(...)`.
- `tests/native/test_native_ceres_thermodynamic_regression.py` asserts
  Born-SSM+DS `d_born` and `f_solv` report `backend_unavailable` for the current
  Ceres thermodynamic slice instead of falling back to Backend unavailables.
- `src/epcsaft/benchmarks/native_regression.py` keeps `native_binary_kij_tiny`
  as the generic binary fixture and changes the reactive fixture to
  `native_reactive_born_ssmds_tiny`.

## Missing Production Derivative Path

The scalar-templated or analytic derivative path still needed is:

`NativeThermoCeresCostFunction::Evaluate` ->
`evaluate_native_thermo_regression_rows` ->
`chemical_equilibrium_native` / `evaluate_chemical_equilibrium_residual_native` ->
`activity_coefficients` ->
`ePCSAFTStateNative::activity_coefficient_native` ->
`residual_chemical_potential_result_cpp` ->
`composition_contribution_result_cpp` ->
`ares_contributions_cpp` ->
`born_intermediate_state_cpp` / `dadx_born_cpp` for Born-SSM+DS `d_born` and
`f_solv`, including pressure-closure sensitivities through `solve_density_scoped`
when activity or concentration standard states require density solves.

## Validation

- `uv run python scripts/build_epcsaft.py --build-only --parallel 10`: pass
- `uv run python run_pytest.py tests/native/test_native_ceres_thermodynamic_regression.py tests/native/test_cppad_eos_derivatives.py tests/api/test_runtime.py -q`: 42 passed
- `uv run python run_pytest.py tests/native/test_native_ceres_thermodynamic_regression.py tests/native/test_native_reactive_regression.py tests/native/test_native_regression_types.py tests/api/test_runtime.py tests/workflows/test_benchmark_native_regression.py -q`: 55 passed
- `uv run ruff check src/epcsaft/runtime.py src/epcsaft/benchmarks/native_regression.py tests/native/test_native_ceres_thermodynamic_regression.py tests/native/test_native_reactive_regression.py tests/native/test_native_regression_types.py tests/api/test_runtime.py tests/workflows/test_benchmark_native_regression.py`: pass
- `uv run black --check src/epcsaft/runtime.py src/epcsaft/benchmarks/native_regression.py tests/native/test_native_ceres_thermodynamic_regression.py tests/native/test_native_reactive_regression.py tests/native/test_native_regression_types.py tests/api/test_runtime.py tests/workflows/test_benchmark_native_regression.py`: pass

