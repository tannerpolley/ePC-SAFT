# T090 Judge Receipt

## Decision

Approved.

## Rationale

The five architecture candidates have executable implementation changes and focused tests, not documentation-only scaffolding.

- T010/T012: typed Equilibrium Problem dispatch and route diagnostics are implemented in `equilibrium.py` and `epcsaft.py`.
- T020: native route rejection diagnostics are centralized in `equilibrium_core/native_results.py`.
- T030: Target Dataset and reactive regression share target-family summary compilation.
- T040: `ParameterSet` owns runtime compilation through `to_runtime_dict()` and dataset loading adapts through `ParameterSet.from_dataset()`.
- T050: capabilities derive Ipopt routes, problem classes, regression keys, and derivative row counts from registered evidence.

Reactive regression and staged paths remain bounded as residual/diagnostic contexts, not production optimizers.

## Missing Evidence

- Final validation, cleanup hook, clean git status, and local commit are still required.
- The known long-running native Ipopt API/electrolyte lane remains timeout-limited and must not be used as proof unless rerun successfully or explicitly recorded as bounded risk.

## Final Validation Plan

- `uv run python scripts/dev/validate_project.py quick`
- `uv run python scripts/dev/check_text_gates.py`
- `git diff --check`
- `pwsh.exe -NoProfile -ExecutionPolicy Bypass -File "$env:USERPROFILE\.codex\hooks\codex-cleanup.ps1" -RepoRoot .`
- `git status --short --branch`
- Create a local commit if validation passes.
