# Generic implicit sensitivity framework for solved states

## Original Request

Run Goal Prep for task key C from `docs/roadmaps/agent_prompts/index.yaml`: Generic implicit sensitivity framework for solved states.

## Current Status

IMPLEMENTED_VALIDATED

The bounded watcher observed Task B as merged, rebased `codex/generic-implicit-sensitivity-framework` onto `origin/main`, and stopped the heartbeat automation. Task C implementation has started and completed for the first generic package-owned slice.

Watcher settings:

- `watcher_mode: bounded`
- `auto_start_after_gate: true`
- `poll_interval_seconds: 120`
- `max_wait_minutes: 480`
- heartbeat automation: `task-c-dependency-watcher` (deleted after gate pass)

## Goal

Prepare the branch and GoalBuddy board for Task C so bounded watcher runs can resume safely after Task B is merged, update/fast-forward this branch to `origin/main`, and start implementation without asking.

Implement reusable implicit-sensitivity machinery for solved internal states while keeping `epcsaft` general-purpose. The completed slice adds a generic `ImplicitSolveResult`, an implicit Jacobian solve helper, explicit `backend_unavailable` payloads, runtime blocker metadata, and reactive-speciation solved-state diagnostics.

## Non-Negotiable Constraints

- Current branch must be `codex/generic-implicit-sensitivity-framework`.
- Do not start implementation until Task B is merged and the dependency gate passes.
- Do not create local `.worktrees/`.
- Do not use PR #56 as a base.
- Do not modify PR #56.
- Do not add application-specific public APIs.
- No finite difference.
- Do not tape iterative solver loops as production derivatives.
- Do not merge broad equilibrium rewrites into this issue.
- Use only generic package concepts and APIs.

## Dependency Gate

Before implementation, rerun the gate:

```powershell
git fetch origin --prune
git branch --show-current
git status --short
git rebase origin/main
```

Also verify:

- Task A / PR #98 remains merged;
- Task B is merged;
- `origin/main` contains the dependency merge commits;
- current branch is `codex/generic-implicit-sensitivity-framework`;
- branch rebases cleanly on `origin/main`.

If any dependency is missing or the rebase conflicts, stop with status `BLOCKED_DEPENDENCY_OR_REBASE`.

Gate result for this run: `GATE_PASSED`. The watcher rebased the assigned branch onto `origin/main` after Task B merged.

Do not stop at `PREPARED_READY`. With `auto_start_after_gate: true`, once Task B merges and the branch is updated, continue directly into implementation without asking.

## Bounded Watcher

When Task B is not merged, run:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File docs/goals/generic-implicit-sensitivity-framework/notes/watch_dependency.ps1 -PollIntervalSeconds 120 -MaxWaitMinutes 480
```

The watcher should poll every 120 seconds and time out after 480 minutes with `PREPARED_WAITING`. When Task B merges, rerun the dependency gate, update/fast-forward the branch to `origin/main`, and continue into implementation without asking.

The Codex thread heartbeat `task-c-dependency-watcher` was deleted after the gate passed so it will not keep polling stale instructions.

## Unblocked Implementation Scope

- association site fractions;
- density root;
- speciation solve;
- VLE root;
- LLE phase split;
- reactive LLE solve;
- `ImplicitSolveResult` carrying state, residual, jacobians, sensitivity, backend, status, and diagnostics.

Allowed solved-state derivative backends are `analytic_implicit` and `cppad_implicit`. Existing validated local paths may keep `legacy_eigen_forward`; explicitly out-of-scope workflows may report `backend_unavailable`.

## Validation When Unblocked

```powershell
uv run python run_pytest.py tests/native/test_association_implicit_derivative_contract.py tests/native/test_reactive_speciation_implicit_sensitivity.py -q
```

```powershell
uv run python scripts/validate_project.py quick
```

## Canonical Board

The board truth is `docs/goals/generic-implicit-sensitivity-framework/state.yaml`.

Start later with:

```text
/goal Follow docs/goals/generic-implicit-sensitivity-framework/goal.md.
```

