## T057 Post-Bubble-Slice Audit

Date: 2026-05-11

### Decision

After `T056`, the next real derivative blocker is **not** the old pressure-composition substrate and **not** generic reactive binary sensitivities. The next blocker is the broader bubble state:

- multi-vapor composition unknowns
- normalization coupling across vapor species
- vapor-composition targets through solved `y`

### Verified current state

Supported now:

- `reactive_speciation` thermo rows on the existing ideal / concentration / supported activity slices
- `reactive_electrolyte_bubble` thermo rows for:
  - exactly one neutral vapor species
  - pressure targets only
  - supported Born/SSM+DS parameter kinds

Still unsupported:

- multi-vapor `reactive_electrolyte_bubble` thermo rows
- `vapor_composition` bubble targets in the thermo regression path
- broader bubble unknown/residual systems beyond the single scalar log-pressure slice

### Evidence

- `src/epcsaft/native/regression/thermo_regression.cpp`
  - `supported_single_vapor_bubble_row(...)`
  - bubble gating still explicitly restricts to one vapor species and pressure targets
- `src/epcsaft/runtime.py`
  - `production_blockers`
  - `missing_bubble_derivative_residuals`
- `tests/api/test_runtime.py`
  - confirms the remaining missing bubble residual list starts with multi-vapor composition unknowns
- `tests/native/test_native_ceres_thermodynamic_regression.py`
  - proves the single-vapor pressure slice runs as `ceres + cppad_implicit`

### Recommendation

Next worker tranche should target:

1. a multi-vapor bubble state Jacobian
2. solved-`y` sensitivity propagation
3. at least one real `vapor_composition` target in the native thermo regression path

Do not reopen:

- generic reactive `k_ij` thermo routing
- already-finished runtime/helper CppAD substrate work
- the single-vapor pressure-only bubble slice
