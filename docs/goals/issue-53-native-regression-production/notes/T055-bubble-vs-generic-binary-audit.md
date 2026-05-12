## T055 Bubble vs Generic Binary Audit

Date: 2026-05-11

### Decision

The next honest derivative tranche after `T054` is **reactive-electrolyte bubble differentiation**, not generic reactive-speciation `k_ij` / `l_ij` / `k_hb_ij` sensitivities.

### Why generic reactive binary sensitivities are not the next tranche

1. The current native thermodynamic regression path already applies `k_ij` values through `apply_native_thermo_parameters(...)` in
   `src/epcsaft/native/regression/thermo_regression.cpp`, but the derivative adapter for
   `reactive_speciation` rows is intentionally narrower: it calls
   `component_activity_parameter_derivative_result_cpp(...)`, which currently derives only the supported Born/SSM+DS activity slice.

2. Runtime capabilities explicitly still route generic binary interaction work elsewhere:
   - `src/epcsaft/runtime.py`
     - `generic_binary_regression_path = use fit_binary_pair(...) for k_ij, l_ij, and k_hb_ij against direct binary VLE composition data`
     - `blocked_parameter_kinds.k_ij`
     - `blocked_parameter_kinds.l_ij`
     - `blocked_parameter_kinds.k_hb_ij`

3. The user explicitly narrowed the intended role of these parameters: `d_born` and `f_solv` are the reactive Born/SSM+DS regression focus, while `k_ij` and other binary interaction parameters should remain generic binary-regression parameters against direct binary composition data.

4. Wiring `k_ij` through the reactive-speciation thermo path is possible in principle, but it would be additive scope, not the highest-value blocker remaining for issue #53 after the supported CppAD runtime/helper migration and Born tranche.

### Why bubble is now the next blocker

1. Runtime capability metadata still records bubble derivatives as missing:
   - `src/epcsaft/runtime.py`
     - `missing_bubble_derivative_residuals`
     - `unsupported_status = backend_unavailable`

2. The native thermo regression runtime still enters the bubble solve path in
   `src/epcsaft/native/regression/thermo_regression.cpp` but does not provide the derivative call graph needed to move bubble rows into production Ceres/CppAD support.

3. Existing tests confirm only the implicit sensitivity helper contract, not production bubble regression derivatives:
   - `tests/native/test_cppad_bubble_derivatives.py`
   - `tests/api/test_runtime.py`

### Next tranche recommendation

Create the next worker tranche around **bubble derivative architecture**, starting with a narrow supported slice and an explicit unknown/residual layout:

1. define the supported bubble residual vector and unknown ordering
2. isolate the fixed-liquid / vapor-composition / solved-pressure derivative call graph
3. implement the first CppAD/implicit bubble regression slice only after that layout is locked

### Suggested validation slice for the next tranche

- `uv run python run_pytest.py tests/native/test_cppad_bubble_derivatives.py tests/native/test_native_ceres_thermodynamic_regression.py tests/api/test_runtime.py -q`
- any new package-owned bubble benchmark or derivative-check script added by the tranche
