# T018 Implicit Sensitivity Receipt

## Result

Added a native implicit-sensitivity primitive for nested solves:

- `src/epcsaft/native/regression/implicit_sensitivity.h`
- `src/epcsaft/native/regression/implicit_sensitivity.cpp`
- `_core._solve_native_implicit_sensitivity(...)`

The native solve implements:

```text
R(u, theta) = 0
u_theta = -R_u^{-1} R_theta
```

It reports explicit failure states (`invalid_input`, `backend_unavailable`, `nonfinite_objective`, `singular_jacobian`) and always reports `finite_difference_used=false`.

## Important Limit

This is the native linear algebra primitive T019 can use. The reactive speciation and electrolyte bubble row evaluator still needs to expose its converged inner residual Jacobians before full production sensitivities are wired through real row solves.

## Validation

- `uv run python scripts/build_epcsaft.py --build-only --parallel 10`
  - pass
- `uv run python run_pytest.py tests/native/test_cppad_reactive_speciation_derivatives.py tests/native/test_cppad_bubble_derivatives.py -q`
  - pass; `4 passed`
- `uv run ruff check tests/native/test_cppad_reactive_speciation_derivatives.py tests/native/test_cppad_bubble_derivatives.py`
  - pass
- `uv run black --check tests/native/test_cppad_reactive_speciation_derivatives.py tests/native/test_cppad_bubble_derivatives.py`
  - pass
