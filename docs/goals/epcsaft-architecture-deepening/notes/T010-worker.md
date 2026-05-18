# T010 Worker Receipt

Result: done

## Summary

Deepened typed **Equilibrium Problem** objects so they own route classification and route-diagnostic threading for their public solve path. Problem objects now attach `equilibrium_route` and `route_reason` on success and preserve that route context when a native solver raises `SolutionError`.

## Changed Files

- `src/epcsaft/equilibrium.py`
- `tests/equilibrium/core/test_vle.py`

## Behavior

- Added internal `_solve_problem_route(...)` support in `equilibrium.py`.
- Updated `TPFlash`, `StabilityAnalysis`, `BubblePoint`, `DewPoint`, `LLEProblem`, `ElectrolyteLLEProblem`, and `ElectrolyteBubblePoint` to solve through the typed route gate.
- Preserved public string API compatibility.
- Added direct `solve_equilibrium(epcsaft.TPFlash(...))` coverage proving typed problem route diagnostics.

## Verification

Passed:

- `uv run python run_pytest.py tests/equilibrium/core/test_vle.py tests/equilibrium/core/test_api.py tests/equilibrium/core/test_lle.py tests/api/package/test_downstream_integration_smokes.py -q` -> `42 passed`
- `uv run python run_pytest.py tests/api/reactive/test_reactive_phase_equilibrium_problem_routes_native.py tests/api/reactive/test_staged_reactive_route_not_production.py -q` -> `4 passed`
- `uv run python scripts/dev/check_text_gates.py` -> passed
- `git diff --check` -> passed

Timeout-limited inherited from T020:

- `tests/api/equilibrium/test_electrolyte_lle_problem_native_ipopt.py` still exceeds 5 minutes in this worktree.

## Remaining Risk

The slice does not replace the full `mixture.equilibrium(kind=...)` string dispatcher. It deepens the typed problem-object path and keeps the existing string API stable. Further dispatcher consolidation can be considered only if a later Judge asks for it as part of a larger public-API cleanup.
