# Generic VLE/Fugacity-Equilibrium Solver for Volatile Neutral Species

## Objective

Prepare and execute Task G for a generic VLE/fugacity-equilibrium solver for volatile neutral species on branch `codex/generic-vle-fugacity-equilibrium`.

## Original Request

Run Goal Prep for task key G: Generic VLE/fugacity-equilibrium solver for volatile neutral species. Use the local live GoalBuddy board in Codex. Read GitHub issue #90, record the current issue scope in `docs/goals/generic-vle-fugacity-equilibrium/notes/issue_scope.md`, and do not start implementation until the dependency gate passes.

## Intake Summary

- Input shape: `existing_plan`
- Audience: `epcsaft` maintainers and downstream modeling agents
- Authority: `requested`
- Proof type: `test`
- Completion proof: a focused generic VLE/fugacity implementation is validated, diagnostics report the route used, no finite-difference fallback exists, and a focused PR closes issue #90.
- Likely misfire: turning the work into an application-specific API, assuming vapor-phase ions without explicit model support, or falling back to finite differences.
- Blind spots considered:
  - Dependency C is already merged into `origin/main` via PR #104, so the gate is open.
  - The issue scope includes direct volatile partial pressure from liquid fugacity plus bubble/dew/TP flash routes and route diagnostics.
  - The package must stay generic and application-neutral.
  - No finite difference is allowed.
  - Solved-state sensitivities should use the approved implicit or analytic derivative policy, not taped solver loops.
- Existing plan facts:
  - Task key: G.
  - Prompt file: `docs/roadmaps/agent_prompts/G_generic_vle_fugacity_equilibrium_solver_for_volatile_neutral_species.md`.
  - Assigned branch: `codex/generic-vle-fugacity-equilibrium`.
  - Dependency: C.
  - Dependency PR #104 is merged into `origin/main`.
  - Current branch is on `origin/main` and rebased cleanly.
  - Use bounded watcher mode and auto-start after gate.
  - Do not create local `.worktrees/`.
  - Do not add application-specific public APIs.
  - No finite difference.

## Goal Kind

`existing_plan`

## Current Tranche

Dependency C is satisfied, so the board is ready for bounded Scout/Judge/Worker execution on the assigned branch. Continue into the first generic VLE/fugacity-equilibrium implementation slice without asking again.

## Non-Negotiable Constraints

- Stay on branch `codex/generic-vle-fugacity-equilibrium`.
- Do not create local `.worktrees/`.
- Do not add application-specific public APIs.
- Keep `epcsaft` generic and reusable.
- No finite difference.
- Do not assume ions distribute to vapor unless the model explicitly includes that behavior.
- Use the approved derivative policy for solved states and explicit algebraic derivatives.
- Validate with the task-listed pytest and project quick-validation commands once implementation lands.

## Stop Rule

Stop with `BLOCKED_DEPENDENCY_OR_REBASE` if the branch cannot rebase or fast-forward cleanly to `origin/main`.

## Canonical Board

Machine truth lives at:

`docs/goals/generic-vle-fugacity-equilibrium/state.yaml`

If this charter and `state.yaml` disagree, `state.yaml` wins for task status, active task, receipts, verification freshness, and completion truth.

## Run Command

```text
/goal Follow docs/goals/generic-vle-fugacity-equilibrium/goal.md.
```

## PM Loop

On every `/goal` continuation:

1. Read this charter.
2. Read `state.yaml`.
3. Read `notes/dependency_gate.yaml`.
4. Re-run the dependency gate from the task prompt before implementation or broader planning.
5. Work only on the active board task.
6. If the gate fails or the branch cannot rebase cleanly, update the gate receipt and stop with `PREPARED_WAITING` or `BLOCKED_DEPENDENCY_OR_REBASE`.
7. If the gate passes, continue implementation automatically because `auto_start_after_gate: true`.
8. Keep the implementation within the selected Worker `allowed_files`.
