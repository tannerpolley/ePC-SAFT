# Explicit CppAD Parameter Derivatives

## Objective

Prepare Task B, wait only until dependency A is merged, then auto-start implementation under a bounded watcher. Implement or verify explicit CppAD/analytic parameter derivatives for EOS and property APIs without finite differences or application-specific public APIs.

## Original Request

Run Goal Prep for Task B: Explicit CppAD parameter derivatives for EOS/property APIs, using the local live GoalBuddy board, on branch `codex/cppad-explicit-parameter-derivatives`, and do not start implementation until dependency A passes.

## Intake Summary

- Input shape: `existing_plan`
- Audience: ePC-SAFT package maintainers and downstream regression/equilibrium users
- Authority: `requested`
- Proof type: `test`
- Completion proof: Task B has a focused PR with generic explicit parameter derivative support, CppAD-enabled build passes, targeted native CppAD derivative tests pass, `validate_project.py quick` passes, and the PR body records issue link, summary, tests, limitations, and next dependencies.
- Likely misfire: Starting implementation before Task A is merged, adding application-specific public APIs, or using finite-difference derivatives while claiming CppAD coverage.
- Blind spots considered: Dependency A must be merged into `origin/main`; this branch must rebase cleanly before code work; solved-state implicit sensitivities are explicitly out of scope; PR #56 must not be used as a base.
- Existing plan facts: Task key `B`; branch `codex/cppad-explicit-parameter-derivatives`; dependency `A`; `watcher_mode: bounded`; `auto_start_after_gate: true`; `poll_interval_seconds: 120`; `max_wait_minutes: 480`; package must stay general-purpose; CppAD is the default for explicit algebraic derivatives; no finite differences; no solver-loop tapes in production derivatives.

## Goal Kind

`existing_plan`

## Current Tranche

Implementation tranche: Task A / PR #98 is merged into `origin/main`, the dependency gate has passed, and bounded watcher auto-start is enabled. Continue from the active board task without stopping at `PREPARED_READY`.

## Coordinator Scope Repair

The worktree coordinator repaired the local GoalBuddy automation paths on 2026-05-13. Auxiliary files now belong under `docs/goals/explicit-cppad-parameter-derivatives/notes/`.

Task B must finish with `PR_OPENED` or `BLOCKED_SCOPE_GAP`; it must not stop with only a vague `k_hb_ij` limitation and no PR. If full active-association property/regression derivatives through `k_hb_ij` require implicit association site-fraction sensitivity, classify that as `blocker_requires_implicit_association_sensitivity`, include `k_hb_ij` in coverage/capability, add a test preventing silent omission or overclaim, and identify Task C as the owner.

## Non-Negotiable Constraints

- Do not implement while dependency A is unmerged.
- Do not create local `.worktrees/`.
- Do not use PR #56 as a base.
- Do not modify PR #56.
- Do not add application-specific public APIs.
- Do not use finite differences.
- Do not implement implicit solved-state sensitivities in Task B.
- Do not tape iterative solver loops as production derivatives.
- Keep public package APIs generic.
- Validate in proportion to the active task and run the repo cleanup hook before completion reporting.

## Stop Rule

If dependency A is not merged, write `PREPARED_WAITING` and stop.

When the dependency gate passes, do not stop at `PREPARED_READY`; auto-start implementation. Stop only for a documented blocker such as `BLOCKED_REBASE_CONFLICT`, `BLOCKED_MISSING_INPUT`, repeated unclear verification failure, or when a final audit proves the full original outcome is complete.

## Canonical Board

Machine truth lives at:

`docs/goals/explicit-cppad-parameter-derivatives/state.yaml`

If this charter and `state.yaml` disagree, `state.yaml` wins for task status, active task, receipts, verification freshness, and completion truth.

## Dependency Gate

Gate details live at:

`docs/goals/explicit-cppad-parameter-derivatives/dependency_gate.yaml`

## Run Command

```text
/goal Follow docs/goals/explicit-cppad-parameter-derivatives/goal.md.
```

## PM Loop

On every `/goal` continuation:

1. Read this charter.
2. Read `state.yaml`.
3. Read `dependency_gate.yaml`.
4. Run the bundled GoalBuddy update checker when available and mention a newer version without blocking.
5. Re-check that dependency A is merged, `origin/main` contains the dependency merge, the current branch is correct, and the branch rebases cleanly on `origin/main`.
6. Work only on the active board task.
7. Write a compact task receipt.
8. Update the board.
9. Do not activate implementation work until the dependency gate passes.
10. With `watcher_mode: bounded` and `auto_start_after_gate: true`, continue implementation immediately after the gate passes instead of recording `PREPARED_READY`.
