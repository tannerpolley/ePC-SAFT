# Generic electrolyte LLE with distributed ions

## Objective

Prepare and execute the issue #92 tranche for generic electrolyte LLE with distributed ions on the generic LLE foundation, while keeping `epcsaft` general-purpose and avoiding application-specific public APIs.

## Original Request

Run Goal Prep for task I: Generic electrolyte LLE with distributed ions, use the local live GoalBuddy board in Codex, and do not start implementation until the dependency gate in the task prompt passes.

## Intake Summary

- Input shape: `existing_plan`
- Audience: repo maintainers and future Codex agents
- Authority: `requested`
- Proof type: `artifact`
- Completion proof: a focused PR closes issue #92 after the required validation passes, with the board and receipts preserving the issue scope, dependency gate, and limitations truthfully.
- Likely misfire: starting implementation before the upstream gate opens, or narrowing the work to board prep only without preserving the actual issue scope and constraints.
- Blind spots considered: dependency H (#91) is still open, the task must remain generic, no finite difference is allowed, and Ascani 2022 Case Study 2 may need documentation instead of a forced pass.
- Existing plan facts: assigned branch `codex/generic-electrolyte-lle`; dependencies `H, C`; local live GoalBuddy board requested; bounded watcher mode required with `auto_start_after_gate: true`, `poll_interval_seconds: 120`, and `max_wait_minutes: 480`; record current issue scope in `docs/goals/<slug>/notes/issue_scope.md`; do not create local `.worktrees/`; do not add application-specific public APIs; no finite difference; open a focused draft PR with `Closes #92`; self-review against `origin/main`; merge only when checks and GoalBuddy audit are complete; clean up branch after merge.
- Existing plan facts: assigned branch `codex/generic-electrolyte-lle`; dependencies `H, C`; local live GoalBuddy board requested; bounded watcher mode required with `auto_start_after_gate: true`, `poll_interval_seconds: 120`, and `max_wait_minutes: 480`; record current issue scope in `docs/goals/<slug>/notes/issue_scope.md`; keep the dependency gate at `docs/goals/<slug>/notes/dependency_gate.yaml`; keep the watcher at `docs/goals/<slug>/notes/watch_dependency.ps1`; do not create local `.worktrees/`; do not add application-specific public APIs; no finite difference; open a focused draft PR with `Closes #92`; self-review against `origin/main`; merge only when checks and GoalBuddy audit are complete; clean up branch after merge.

## Goal Kind

`existing_plan`

## Current Tranche

Prepare the goal root, keep the dependency gate truthful, and wait on dependency H before implementation. Once the gate opens, execute the issue #92 tranche continuously until the full original outcome is complete, with the generic API and derivative constraints preserved.

## Non-Negotiable Constraints

- Keep `epcsaft` general-purpose.
- Do not add application-specific public APIs.
- No finite difference.
- Do not start implementation until dependencies `H` and `C` are merged into `origin/main` and the branch rebases cleanly.
- Use the local live GoalBuddy board and bounded dependency watcher.

## Stop Rule

Stop only when a final audit proves the full original outcome is complete.

Do not stop after planning, discovery, or gate setup if the user asked for working software and a safe Worker task can be activated.

Do not stop after a single verified Worker slice when the broader owner outcome still has safe local follow-up slices.

Do not stop because a slice needs owner input, credentials, production access, destructive operations, or policy decisions. Mark that exact slice blocked with a receipt, create the smallest safe follow-up or workaround task, and continue all local, non-destructive work that can still move the goal toward the full outcome.

## Canonical Board

Machine truth lives at:

`docs/goals/generic-electrolyte-lle/state.yaml`

Dependency automation lives at:

`docs/goals/generic-electrolyte-lle/notes/dependency_gate.yaml`
`docs/goals/generic-electrolyte-lle/notes/watch_dependency.ps1`

If this charter and `state.yaml` disagree, `state.yaml` wins for task status, active task, receipts, verification freshness, and completion truth.

## Run Command

```text
/goal Follow docs/goals/generic-electrolyte-lle/goal.md.
```

## PM Loop

On every `/goal` continuation:

1. Read this charter.
2. Read `state.yaml`.
3. Run the bundled GoalBuddy update checker when available and mention a newer version without blocking.
4. Re-check the intake: original request, input shape, authority, proof, blind spots, existing plan facts, and likely misfire.
5. Work only on the active board task.
6. Assign Scout, Judge, Worker, or PM according to the task.
7. Write a compact task receipt.
8. Update the board.
9. If Judge selected a safe Worker task with `allowed_files`, `verify`, and `stop_if`, activate it and continue unless blocked.
10. If a problem, suggestion, or follow-up should become a repo artifact, create an approved issue/PR or ask the operator whether to create one.
11. Treat a slice audit as a checkpoint, not completion, unless it explicitly proves the full original outcome is complete.
12. Finish only with a Judge/PM audit receipt that maps receipts and verification back to the original user outcome and records `full_outcome_complete: true`.

Issue and PR handoffs are supporting artifacts. `state.yaml` remains authoritative, and every external artifact decision must be recorded in a task receipt.
