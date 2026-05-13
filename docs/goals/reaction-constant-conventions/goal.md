# General reaction and equilibrium-constant convention layer

## Objective

Prepare and, once dependency A is merged, implement the general reaction-constant convention layer for `epcsaft` without adding application-specific APIs or finite-difference derivatives.

## Original Request

Run Goal Prep for task D, "General reaction and equilibrium-constant convention layer," on branch `codex/reaction-constant-conventions`; use the local live GoalBuddy board; do not start implementation until dependency A passes.

## Intake Summary

- Input shape: `existing_plan`
- Audience: `epcsaft` maintainers and downstream package consumers needing generic reaction/equilibrium convention support
- Authority: `requested`
- Proof type: `test`
- Completion proof: Dependency A is satisfied, implementation is completed on this branch, focused tests and quick validation pass, and a focused PR documents summary, tests, limitations, and next dependencies.
- Likely misfire: Treating this as an application-specific MEA/lithium workflow, defaulting reaction-constant fitting, or bypassing the dependency gate.
- Blind spots considered: Dependency A is currently not merged into `origin/main`; public API names must stay generic; derivative handling must not use finite difference; PR #56 must not be used as a base or modified.
- Existing plan facts: Task D depends on A, uses branch `codex/reaction-constant-conventions`, has `auto_start_after_gate: false`, and requires `dependency_gate.yaml` plus `PREPARED_WAITING` when dependencies are missing.

## Goal Kind

`existing_plan`

## Current Tranche

Prepared ready tranche: maintain a valid GoalBuddy board and dependency-gate record for task D. Dependency A is merged into `origin/main`, the assigned branch is correct, and the branch fast-forwarded cleanly to `origin/main`.

The operator explicitly overrode the loaded task prompt on 2026-05-13: `auto_start_after_gate: true`, `watcher_mode: bounded`, do not stop at `PREPARED_READY`, and start implementation immediately after the gate passes.

## Non-Negotiable Constraints

- Current branch must be `codex/reaction-constant-conventions` before Goal Prep writes files.
- Do not create local `.worktrees/`.
- Do not start implementation until dependency A passes.
- Do not use PR #56 as a base, modify it, or close unrelated issues.
- Do not add application-specific public APIs for MEA, lithium extraction, absorption columns, extraction efficiency, distribution coefficient, selectivity, or similar downstream workflows.
- Use generic package concepts such as `equilibrium(...)`, `regress_parameters(...)`, `ReactionSet(...)`, `EquilibriumProblem(...)`, `RegressionProblem(...)`, `TargetDataset(...)`, `PhaseSpec(...)`, and `ParameterSet(...)`.
- No finite difference.
- Do not tape iterative solver loops as production derivatives.
- Do not force reaction-constant fitting as the default.
- Validate, when implementation is eventually allowed, with:
  - `uv run python run_pytest.py tests/api/test_reactive_staged_workflow_contract.py tests/api/test_reaction_constant_conventions.py -q`
  - `uv run python scripts/validate_project.py quick`

## Stop Rule

Stop while dependency status is `PREPARED_WAITING`.

After dependency A passes, proceed into the bounded Scout/Judge/Worker implementation path because the operator explicitly set `auto_start_after_gate: true`.

Implementation is authorized by the operator's 2026-05-13 override. Keep watcher mode bounded and preserve the package-scope, derivative, and PR #56 constraints.

## Canonical Board

Machine truth lives at:

`docs/goals/reaction-constant-conventions/state.yaml`

The dependency gate lives at:

`docs/goals/reaction-constant-conventions/dependency_gate.yaml`

If this charter and `state.yaml` disagree, `state.yaml` wins for task status, active task, receipts, verification freshness, and completion truth.

## Run Command

```text
/goal Follow docs/goals/reaction-constant-conventions/goal.md.
```

## PM Loop

On every `/goal` continuation:

1. Read this charter.
2. Read `state.yaml`.
3. Read `dependency_gate.yaml`.
4. Re-run the dependency gate from the task prompt before any implementation or planning beyond board maintenance.
5. Work only on the active board task.
6. If dependency A is still missing or the branch cannot rebase cleanly on `origin/main`, update the gate receipt and stop with `PREPARED_WAITING` or `BLOCKED_DEPENDENCY_OR_REBASE`.
7. If the gate passes, activate bounded Scout/Judge/Worker tasks and continue without asking.
8. Keep implementation within the selected Worker `allowed_files`.
