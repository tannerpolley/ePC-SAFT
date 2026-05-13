# Generic Target-Row and Dataset Schema

## Objective

Execute Task E: create application-neutral regression and validation target-row schemas for future generic regression and equilibrium workflows.

## Original Request

Run Goal Prep for Task E, Generic target-row and dataset schema, using the local live GoalBuddy board; do not start implementation until dependency A passes.

## Intake Summary

- Input shape: `existing_plan`
- Audience: `epcsaft` package maintainers and future roadmap agents
- Authority: `requested`
- Proof type: `test`
- Completion proof: Generic `TargetDataset`/target-row schema behavior is implemented and validated by `uv run python run_pytest.py tests/api/test_regression_problem_schema.py tests/api/test_regression_api.py -q` and `uv run python scripts/validate_project.py docs`.
- Likely misfire: Treating this as an application-specific schema or optimizer task instead of a generic package schema task gated by dependency A.
- Blind spots considered: dependency A was initially unresolved during prep; Task E now uses bounded watcher auto-run with `auto_start_after_gate: true`; no finite-difference derivative path may be introduced; no PR #56 base or application-specific public APIs are allowed.
- Existing plan facts: Task key E; branch `codex/generic-target-row-schema`; dependency A is Task A / PR #98; bounded watcher mode; `auto_start_after_gate: true`; do not stop at `PREPARED_READY`.

## Goal Kind

`existing_plan`

## Current Tranche

The current tranche is bounded watcher auto-run plus Task E implementation. Dependency A is satisfied by PR #98 merged into `origin/main`. The branch has been fast-forwarded to `origin/main`, and implementation proceeds immediately without asking.

## Non-Negotiable Constraints

- Stay on branch `codex/generic-target-row-schema` before writing or executing Task E work.
- Do not create local `.worktrees/`.
- Do not use PR #56 as a base.
- Do not add application-specific public APIs.
- Keep `epcsaft` general-purpose.
- No finite difference.
- Do not implement optimizer internals in Task E.
- Use generic concepts such as `TargetDataset`, `RegressionProblem`, `EquilibriumProblem`, `ReactionSet`, `PhaseSpec`, and `ParameterSet`.
- Verify dependency A before implementation: PR #98 merged, dependency merge commit present in `origin/main`, correct branch, and clean rebase or approved `git merge --ff-only origin/main` fallback when rebase is blocked by tool approval policy and the branch has no unique commits.

## Stop Rule

Stop with `BLOCKED_REBASE_CONFLICT` if the branch cannot be rebased or fast-forwarded cleanly to `origin/main`.

Do not stop after planning, discovery, or Judge selection if the dependency gate later passes and a safe Worker task can be activated by an explicit follow-up run.

## Canonical Board

Machine truth lives at:

`docs/goals/generic-target-row-schema/state.yaml`

If this charter and `state.yaml` disagree, `state.yaml` wins for task status, active task, receipts, verification freshness, and completion truth.

## Run Command

```text
/goal Follow docs/goals/generic-target-row-schema/goal.md.
```

## PM Loop

On every `/goal` continuation:

1. Read this charter.
2. Read `state.yaml`.
3. Read `dependency_gate.yaml`.
4. Run the bundled GoalBuddy update checker when available and mention a newer version without blocking.
5. Work only on the active board task.
6. Re-check dependency A / PR #98 before implementation.
7. If the gate fails or the branch cannot fast-forward cleanly, write a compact receipt and stop with `BLOCKED_REBASE_CONFLICT`.
8. If the gate passes, continue implementation automatically because `auto_start_after_gate` is true.
