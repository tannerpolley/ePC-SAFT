# Derivative backend completion audit and coverage matrix hard gate

## Objective

Audit and harden derivative coverage and runtime capability reporting so backend derivative gaps are explicit before new equilibrium or regression implementation starts.

## Original Request

Run Goal Prep for task key A, "Derivative backend completion audit and coverage matrix hard gate," on branch `codex/backend-coverage-hard-gate`, using the local live GoalBuddy board. Do not start implementation until the dependency gate in the task prompt passes.

## Intake Summary

- Input shape: `existing_plan`
- Audience: ePC-SAFT package maintainers and downstream agents that depend on truthful derivative backend capabilities.
- Authority: `requested`
- Proof type: `test`
- Completion proof: derivative coverage and runtime capability semantics are updated, targeted derivative/capability tests and quick validation pass, a focused PR is opened, and remaining blockers are documented as follow-up work.
- Likely misfire: treating this as a broad solver implementation task instead of a hard-gate audit and matrix/capability contract task.
- Blind spots considered: dependency gate must be re-run before implementation; no finite difference is allowed; public APIs must remain generic; PR #56 must not be used as a base or modified; any missing dependency or rebase conflict stops implementation.
- Existing plan facts: task A comes from `docs/roadmaps/agent_prompts/index.yaml`; prompt file is `docs/roadmaps/agent_prompts/A_derivative_backend_completion_audit_and_coverage_matrix_hard_gate.md`; dependencies are none; `auto_start_after_gate` is true; assigned branch is `codex/backend-coverage-hard-gate`.

## Goal Kind

`existing_plan`

## Current Tranche

Prepare the GoalBuddy board, re-run the dependency/rebase gate before implementation, then complete successive safe verified slices until task A is done: expand derivative coverage and runtime capability row-family matrices, classify derivative paths as `production_supported`, `blocker`, or `out_of_scope`, align capability reporting with coverage, add tests for coverage semantics, and produce a blocker list if gaps remain.

## Non-Negotiable Constraints

- Work on branch `codex/backend-coverage-hard-gate`.
- Do not create local `.worktrees/`.
- Do not use PR #56 as a base.
- Do not modify PR #56.
- Do not add application-specific public APIs.
- No finite difference.
- Allowed derivative backend labels are `analytic`, `cppad`, `analytic_implicit`, `cppad_implicit`, `legacy_eigen_forward` only for validated legacy/local paths, and `backend_unavailable` only for explicitly out-of-scope workflows.
- Do not tape iterative solver loops as production derivatives.
- Do not implement large solver features for this task.
- Stop with `BLOCKED_DEPENDENCY_OR_REBASE` if the dependency gate or rebase fails.

## Stop Rule

Stop only when a final audit proves the full original outcome is complete.

Do not stop after planning, discovery, or Judge selection if a safe Worker task can be activated.

Do not stop after a single verified Worker slice when the broader owner outcome still has safe local follow-up slices. After each slice audit, advance the board to the next highest-leverage safe Worker task and continue.

Do not stop because a slice needs owner input, credentials, production access, destructive operations, or policy decisions. Mark that exact slice blocked with a receipt, create the smallest safe follow-up or workaround task, and continue all local, non-destructive work that can still move the goal toward the full outcome.

## Canonical Board

Machine truth lives at:

`docs/goals/derivative-backend-coverage-hard-gate/state.yaml`

If this charter and `state.yaml` disagree, `state.yaml` wins for task status, active task, receipts, verification freshness, and completion truth.

## Run Command

```text
/goal Follow docs/goals/derivative-backend-coverage-hard-gate/goal.md.
```

## PM Loop

On every `/goal` continuation:

1. Read this charter.
2. Read `state.yaml`.
3. Run the bundled GoalBuddy update checker when available and mention a newer version without blocking.
4. Re-check the intake, prompt, dependency gate, branch, constraints, and likely misfire.
5. Work only on the active board task.
6. Assign Scout, Judge, Worker, or PM according to the task.
7. Write a compact task receipt.
8. Update the board.
9. If Judge selected a safe Worker task with `allowed_files`, `verify`, and `stop_if`, activate it and continue unless blocked.
10. Treat a slice audit as a checkpoint, not completion, unless it explicitly proves the full original outcome is complete.
11. Finish only with a Judge/PM audit receipt that maps receipts and verification back to the original user outcome and records `full_outcome_complete: true`.
