# T030 Worker Receipt

## Result

Done.

## Summary

Target-family compilation now has one shared summary compiler in `src/epcsaft/regression.py`. Generic regression keeps its existing `target_family_summaries` payload shape, `TargetDataset` can report row-family counts and summaries directly, and reactive regression uses the same compiler for both context schema diagnostics and evaluated residual-family evidence.

The reactive path remains a structured residual evaluator. The new diagnostics report residual counts and residual block norms by target family, but do not claim a native Ceres reactive production optimizer.

## Changed Files

- `src/epcsaft/regression.py`
- `src/epcsaft/reactive_regression.py`
- `tests/api/regression/test_regression_problem_schema.py`
- `tests/api/reactive/test_reactive_regression_execution.py`

## Verification

- `uv run python run_pytest.py tests/api/regression/test_regression_problem_schema.py tests/api/reactive/test_reactive_regression_execution.py -q` -> `10 passed`
- `uv run python run_pytest.py tests/api/regression tests/api/reactive/test_reactive_regression_setup.py tests/api/reactive/test_reactive_regression_execution.py tests/api/reactive/test_reactive_regression_diagnostics.py tests/native/ceres -q` -> `45 passed`
- `uv run python scripts/dev/check_text_gates.py` -> passed
- `git diff --check` -> passed

## Boundaries

- No downstream-specific regression surface was added.
- No reactive optimizer capability claim was broadened.
- Ceres availability was concrete in this worktree for the T030 native regression test lane.
