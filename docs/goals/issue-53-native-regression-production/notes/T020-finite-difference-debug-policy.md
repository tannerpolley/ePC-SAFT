# T020 Finite-Difference Debug Policy Receipt

## Result

Finite-difference derivative use is now debug-gated for the touched production equilibrium/speciation paths:

- Native chemical-equilibrium explicit ``jacobian_backend="finite_difference"`` now requires ``EPCSAFT_ALLOW_FINITE_DIFFERENCE_DEBUG=1``.
- Native chemical-equilibrium ``jacobian_backend="auto"`` no longer silently falls back to finite differences for concentration- or activity-coupled standard states. It reports ``backend_unavailable`` unless the debug gate is enabled.
- Native electrolyte LLE residual finite-difference Jacobian evaluation is also debug-gated.
- Runtime capabilities now report ``jacobian_auto_policy="analytic_ideal_else_backend_unavailable"``, ``finite_difference_requires_explicit_request=True``, and the debug gate name.
- Tests that intentionally exercise finite-difference diagnostics opt in with the debug environment variable.

## Evidence

- Native build:
  - `uv run python scripts/build_epcsaft.py --build-only --parallel 10`
- Focused finite-difference policy/API/native tests:
  - `uv run python run_pytest.py tests/api/test_finite_difference_debug_policy.py tests/api/test_reactive_speciation.py tests/api/test_runtime.py tests/native/test_chemical_equilibrium_native.py tests/native/test_equilibrium_native_contracts.py tests/equilibrium/test_lle.py -q`
  - Result: `118 passed, 1 skipped in 17.63s`
- Lint/format:
  - `uv run ruff check src/epcsaft/runtime.py tests/api/test_finite_difference_debug_policy.py tests/api/test_reactive_speciation.py tests/api/test_runtime.py tests/native/test_chemical_equilibrium_native.py tests/native/test_equilibrium_native_contracts.py tests/equilibrium/test_lle.py`
  - `uv run black --check src/epcsaft/runtime.py tests/api/test_finite_difference_debug_policy.py tests/api/test_reactive_speciation.py tests/api/test_runtime.py tests/native/test_chemical_equilibrium_native.py tests/native/test_equilibrium_native_contracts.py tests/equilibrium/test_lle.py`

## Files

- `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp`
- `src/epcsaft/native/epcsaft_equilibrium.cpp`
- `src/epcsaft/runtime.py`
- `tests/api/test_finite_difference_debug_policy.py`
- `tests/api/test_reactive_speciation.py`
- `tests/api/test_runtime.py`
- `tests/native/test_chemical_equilibrium_native.py`
- `tests/native/test_equilibrium_native_contracts.py`
- `tests/equilibrium/test_lle.py`
- `docs/pages/electrolyte_vle_reactive_workflow.rst`
- `docs/pages/equilibrium_cookbook.rst`
- `docs/pages/parameter_regression.rst`
