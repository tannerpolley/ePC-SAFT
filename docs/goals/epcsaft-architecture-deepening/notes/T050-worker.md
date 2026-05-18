# T050 Worker Receipt

## Result

Done.

## Summary

`capabilities()` now derives its Ipopt public routes, equilibrium route payloads, problem-object class list, and derivative coverage rows from registered capability evidence in `src/epcsaft/runtime.py`.

The exposed capability payload also includes a `capability_evidence` summary tying public routes, problem-object classes, regression keys, and derivative row count back to the registered evidence. The existing public capability claims and route order remain intact.

## Changed Files

- `src/epcsaft/runtime.py`
- `tests/api/runtime/test_runtime_capabilities_dependency_gates.py`

## Verification

- `uv run python run_pytest.py tests/api/runtime/test_runtime_capabilities_dependency_gates.py tests/api/runtime/test_runtime_exports_and_metadata.py -q` -> `15 passed`
- `uv run python run_pytest.py tests/api/runtime tests/api/package tests/native/contracts/test_derivative_coverage_matrix.py tests/native/contracts/test_property_derivative_backend_contract.py -q` -> `55 passed`
- `uv run python scripts/dev/check_text_gates.py` -> passed
- `git diff --check` -> passed

## Boundaries

- No unsupported production capability was promoted.
- Reactive batch regression still reports diagnostic residual context only, not a production Ceres optimizer.
- The Ipopt public route list remains evidence-derived and unchanged in order.
