# T012 Worker Receipt

Result: done

## Summary

Finished the **Equilibrium Problem** follow-up by making explicit public string requests build typed `EquilibriumProblem` objects and then run through `solve_equilibrium(...)`. This removes the duplicated explicit-route ownership from `ePCSAFTMixture.equilibrium(...)` for non-reactive routes while preserving public string API compatibility.

## Changed Files

- `src/epcsaft/equilibrium.py`
- `src/epcsaft/epcsaft.py`
- `tests/equilibrium/core/test_api.py`

## Behavior

- Added `equilibrium_problem_from_request(...)` in `equilibrium.py`.
- Added typed conversion for bubble, dew, TP flash, neutral LLE, electrolyte LLE, electrolyte bubble, neutral stability, and electrolyte stability requests.
- Updated `BubblePoint` to represent fixed-temperature or fixed-pressure bubble-point problems.
- Changed `ePCSAFTMixture.equilibrium(...)` explicit non-reactive branches to create typed problem objects and call `solve_equilibrium(...)`.
- Added coverage proving explicit string requests dispatch through typed problem objects.

## Verification

Passed:

- `uv run python run_pytest.py tests/equilibrium/core -q` -> `73 passed`
- `uv run python run_pytest.py tests/api/package/test_downstream_integration_smokes.py tests/api/reactive/test_reactive_phase_equilibrium_problem_routes_native.py tests/api/reactive/test_staged_reactive_route_not_production.py -q` -> `8 passed`
- `uv run python run_pytest.py tests/api/equilibrium --collect-only -q` -> `1 test collected`
- `uv run python scripts/dev/check_text_gates.py` -> passed
- `git diff --check` -> passed

Timeout-limited:

- `uv run python run_pytest.py tests/equilibrium/core tests/api/equilibrium -q` timed out because `tests/api/equilibrium/test_electrolyte_lle_problem_native_ipopt.py` exceeds the current command window in this worktree.

## Remaining Risk

The `auto` and reactive branches still need route-specific branching because they infer routes or require reaction-specific validation. They remain public adapters around typed or specialized problem flows and should not block the accepted Equilibrium Problem slice.
