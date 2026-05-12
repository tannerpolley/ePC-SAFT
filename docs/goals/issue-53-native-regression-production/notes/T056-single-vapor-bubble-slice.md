## T056 Single-Vapor Bubble Slice

Date: 2026-05-11

### Outcome

Implemented the first real native reactive-electrolyte bubble derivative slice in the thermodynamic Ceres hot loop.

### Supported slice added

- `row_mode = reactive_electrolyte_bubble`
- exactly one neutral vapor species
- `pressure` targets only
- supported Born/SSM+DS parameter kinds already on the native thermo path:
  - `born_radius`
  - `born_diameter`
  - `f_solv`
  - `solvation_factor`

### Implementation

In `src/epcsaft/native/regression/thermo_regression.cpp`:

1. Added bubble-row structure gating with `supported_single_vapor_bubble_row(...)`.
2. Extended `thermo_derivative_supported(...)` so the narrow bubble slice is allowed while broader bubble combinations stay blocked.
3. Added `fill_single_vapor_bubble_jacobian(...)`:
   - solves the inner bubble problem with `electrolyte_bubble_pressure_native(...)`
   - reconstructs liquid and vapor native states at the solved pressure
   - uses native `dlnphi/drho` and `dP/drho` helpers to form the single implicit log-pressure state derivative
   - uses native `dlnphi/dtheta` helpers for the supported Born/SSM+DS parameter kinds
   - writes the pressure-target Jacobian directly into the Ceres hot loop with no backend-unavailable production fallback
4. Routed `NativeThermoCeresCostFunction::Evaluate(...)` to this Jacobian builder for bubble rows.
5. Updated row diagnostics so the supported bubble slice reports `derivative_backend = autodiff` instead of `not_differentiated`.

### Capability/documentation alignment

Updated `src/epcsaft/runtime.py` to reflect:

- a real supported bubble slice now exists
- the broader bubble problem is still incomplete
- generic reactive `k_ij` / `l_ij` / `k_hb_ij` remains outside this thermo slice

### Validation

Passed:

- `uv run python scripts/build_epcsaft.py --build-only --parallel 10`
- `uv run python run_pytest.py tests/native/test_cppad_bubble_derivatives.py tests/native/test_native_ceres_thermodynamic_regression.py tests/api/test_runtime.py -q`
  - `50 passed`
- `uv run python run_pytest.py tests/workflows/test_benchmark_native_regression.py -q`
  - `7 passed`

Also probed a real local single-vapor pressure-fit case:

- `reactive_electrolyte_bubble`
- vapor species `["H2O"]`
- parameter `H2O.f_solv`
- result: `status=converged`, `optimizer_backend=ceres`, `derivative_backend=cppad_implicit`

### Remaining gap after T056

This does **not** finish bubble differentiation in general.

Still unsupported:

- multi-vapor bubble unknowns and normalization
- vapor-composition targets through solved `y`
- broader bubble derivative call graphs beyond the single-vapor pressure-only slice
- any bubble slice that would require fake Backend unavailables or unsupported parameter kinds

