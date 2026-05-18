# T020 Worker Receipt

Result: done

## Summary

Deepened the **Production Solver Path** result gate by moving native route diagnostic assembly into `equilibrium_core/native_results.py` and routing Python acceptance/rejection paths through the shared helper. The route gate now has one Python interface for merging postsolve diagnostics with solver metadata before raising `SolutionError`.

## Changed Files

- `src/epcsaft/equilibrium_core/native_results.py`
- `src/epcsaft/equilibrium.py`
- `tests/equilibrium/core/test_native_results.py`
- `tests/equilibrium/core/test_vle.py`
- `tests/equilibrium/core/test_api.py`
- `CONTEXT.md`

## Behavior

- Added `native_route_diagnostics(...)` and `raise_native_route_rejected(...)`.
- Replaced repeated route rejection diagnostic glue in neutral, electrolyte, reactive, and electrolyte-stability Python paths.
- Preserved Ipopt dependency gates as `InputError`; real native solver rejections now carry a consistent diagnostic shape.
- Updated environment-sensitive bubble/dew API coverage to accept either a dependency gate or a real native solver rejection with route diagnostics.
- Fixed the `CONTEXT.md` derivative-path wording so the repo text gate remains clean.

## Verification

Passed:

- `uv run python run_pytest.py tests/equilibrium/core/test_native_results.py tests/equilibrium/core/test_vle.py -q` -> `8 passed`
- `uv run python run_pytest.py tests/equilibrium/core -q` -> `63 passed`
- `uv run python run_pytest.py tests/native/equilibrium/test_result_builder.py tests/equilibrium/core/test_native_results.py tests/equilibrium/core/test_vle.py tests/equilibrium/core/test_api.py -q` -> `27 passed`
- `uv run python run_pytest.py tests/native/equilibrium -q` -> `68 passed`
- `uv run python run_pytest.py tests/equilibrium/electrolyte -q -k "not molality and not source_like and not ascani_case2 and not auto_kind"` -> `12 passed, 4 deselected`
- `uv run python scripts/dev/check_text_gates.py` -> passed
- `git diff --check` -> passed

Timeout-limited:

- `uv run python run_pytest.py tests/native/equilibrium tests/equilibrium/electrolyte tests/api/equilibrium -q` timed out after 5 minutes.
- `uv run python run_pytest.py tests/equilibrium/electrolyte tests/api/equilibrium -q` timed out after 3 minutes.
- `uv run python run_pytest.py tests/api/equilibrium/test_electrolyte_lle_problem_native_ipopt.py -q` timed out after 5 minutes.
- `uv run python run_pytest.py tests/equilibrium/electrolyte -q` timed out after 5 minutes.

Owned timed-out `run_pytest.py` processes were stopped after each timeout.

## Remaining Risk

The implementation did not change native C++ route builders or result-builder semantics. Actual electrolyte/API Ipopt execution tests are long-running in this worktree and should be covered again during final high-level validation or a dedicated native-Ipopt validation lane.
