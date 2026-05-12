# T026 Reactive-Electrolyte Bubble Derivative Gap

## Decision

T026 is stopped on the task stop condition. `reactive_electrolyte_bubble` rows
are native-evaluated, but they are value-only regression rows today. They do not
have a production Ceres derivative path, and adding a backend-unavailable pressure
or vapor-composition sensitivity would violate the Issue #53 backend-unavailable
debug-only policy.

## Evidence

- `src/epcsaft/native/regression/thermo_regression.cpp` rejects non-`reactive_speciation`
  rows in the Ceres derivative gate.
- The same file evaluates `reactive_electrolyte_bubble` rows through
  `electrolyte_bubble_pressure_native(...)`, but row diagnostics explicitly set
  `derivative_backend = "not_differentiated"`.
- The Ceres cost function fills Jacobians only through
  `fill_implicit_speciation_jacobian(...)`, which uses
  `evaluate_chemical_equilibrium_residual_native(...)` and
  `solve_native_implicit_sensitivity(...)`.
- `src/epcsaft/native/epcsaft_equilibrium.cpp` has a native bubble solver and
  best-point diagnostics, but no exported bubble residual/Jacobian payload that
  Ceres can consume.
- `tests/native/test_native_ceres_reactive_pressure_speciation.py` now asserts
  that bubble pressure row evaluation works and that fitting such a row reports
  `backend_unavailable`.

## Missing Residual Equations

The missing production implicit system needs residuals and Jacobians for:

- log-pressure unknown for bubble pressure;
- vapor-composition unknowns plus vapor normalization;
- liquid fixed-composition fugacity/activity state;
- vapor submixture fugacity state;
- fugacity equality residuals
  `log(y_i) + log(phi_i^vap) - log(x_i) - log(phi_i^liq)`;
- pressure target residual through the solved pressure.

That residual interface also needs parameter sensitivities for the fitted
parameter families and continuation variables before it can be wired into
`solve_native_implicit_sensitivity(...)` and Ceres.

## Validation

- `uv run python run_pytest.py tests/native/test_native_ceres_reactive_pressure_speciation.py tests/native/test_cppad_bubble_derivatives.py tests/api/test_runtime.py -q`: 41 passed

