# Generic reactive LLE and chemical phase equilibrium

## Objective

Prepare the goal board and then implement a generic reactive LLE and chemical phase equilibrium path for `epcsaft` only after the dependency gate passes.

## Original Request

Run Goal Prep for task key J: Generic reactive LLE and chemical phase equilibrium. Use the local live GoalBuddy board in Codex, do not start implementation until the dependency gate passes, and record the current issue scope for GitHub issue #93 during `/goal` execution.

## Intake Summary

- Input shape: `existing_plan`
- Audience: `epcsaft` maintainers and downstream modeling agents
- Authority: `requested`
- Proof type: `artifact`
- Completion proof: GitHub issue #93 is closed by a focused merged PR, the branch rebases cleanly on `origin/main`, required validation passes, and the final GoalBuddy audit records `full_outcome_complete: true`.
- Likely misfire: Starting implementation before the dependency gate passes, narrowing the work to an application-specific workflow, or introducing finite-difference derivatives.
- Blind spots considered:
  - The current issue scope must be read from GitHub during `/goal` execution and recorded in `docs/goals/generic-reactive-lle/notes/issue_scope.md`.
  - Dependencies I, D, and F must all be satisfied before implementation.
  - The branch must stay generic and avoid application-specific public APIs.
  - No finite difference is allowed.
- Existing plan facts:
  - Task key J maps to `docs/roadmaps/agent_prompts/J_generic_reactive_lle_and_chemical_phase_equilibrium.md`.
  - Assigned branch is `codex/generic-reactive-lle`.
  - GitHub issue is #93.
  - Dependencies are I, D, and F.
  - The task uses bounded watcher mode with `auto_start_after_gate: true`, `poll_interval_seconds: 120`, and `max_wait_minutes: 480`.
  - Do not create local `.worktrees/`.
  - Do not add application-specific public APIs.
  - No finite difference.

## Goal Kind

`existing_plan`

## Current Tranche

Prepare the local GoalBuddy board, capture issue scope during `/goal`, verify the dependency gate, and then run the gated implementation continuously until the full owner outcome is complete.

## Non-Negotiable Constraints

- Current branch must be `codex/generic-reactive-lle` before Goal Prep writes files.
- Do not create local `.worktrees/`.
- Do not start implementation until the dependency gate passes.
- Do not add application-specific public APIs.
- No finite difference.
- Record the current issue scope in `docs/goals/generic-reactive-lle/notes/issue_scope.md` during `/goal` execution.
- Use the local live GoalBuddy board in Codex.

## Stop Rule

Stop only when a final audit proves the full original outcome is complete.

Do not stop after planning, discovery, or Judge selection if the user asked for working software or automation and a safe Worker task can be activated.

Do not stop after a single verified Worker slice when the broader owner outcome still has safe local follow-up slices.

Do not stop because a slice needs owner input, credentials, production access, destructive operations, or policy decisions. Mark that exact slice blocked with a receipt, create the smallest safe follow-up or workaround task, and continue all local, non-destructive work that can still move the goal toward the full outcome.

## Canonical Board

Machine truth lives at:

`docs/goals/generic-reactive-lle/state.yaml`

If this charter and `state.yaml` disagree, `state.yaml` wins for task status, active task, receipts, verification freshness, and completion truth.

## Run Command

```text
/goal Follow docs/goals/generic-reactive-lle/goal.md.
```

## PM Loop

On every `/goal` continuation:

1. Read this charter.
2. Read `state.yaml`.
3. Read `notes/dependency_gate.yaml`.
4. Re-run the dependency gate from the task prompt before any implementation or planning beyond board maintenance.
5. Work only on the active board task.
6. If dependency I, D, or F is still missing or the branch cannot rebase cleanly on `origin/main`, update the gate receipt and stop with `PREPARED_WAITING` or `BLOCKED_DEPENDENCY_OR_REBASE`.
7. If the gate passes, activate bounded Scout/Judge/Worker tasks and continue without asking.
8. Keep implementation within the selected Worker `allowed_files`.
