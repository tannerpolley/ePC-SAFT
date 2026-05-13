# Generic speciation solver using ePC-SAFT activities

## Original Request

Run Goal Prep for task key F from `docs/roadmaps/agent_prompts/index.yaml`: Generic speciation solver using ePC-SAFT activities.

## Current Status

PREPARED_READY

The assigned branch is checked out, the roadmap and issue were read, the dependency issues are closed, and the dependency gate is recorded as passed. Implementation is not started yet in Goal Prep.

Watcher settings:

- `watcher_mode: bounded`
- `auto_start_after_gate: true`
- `poll_interval_seconds: 120`
- `max_wait_minutes: 480`

## Goal

Prepare the branch and GoalBuddy board for Task F so implementation can start on the assigned branch without asking, while keeping `epcsaft` general-purpose and avoiding application-specific public APIs.

Implement generic liquid reactive/speciation solving with ideal/apparent and activity-based modes using ePC-SAFT activities where requested, while preserving the package-wide derivative policy and diagnostic honesty.

## Non-Negotiable Constraints

- Current branch must be `codex/generic-activity-speciation`.
- Do not start implementation until the dependency gate passes.
- Do not create local `.worktrees/`.
- Do not add application-specific public APIs.
- No finite difference.
- Do not fit reaction constants by default.
- Keep scope generic and package-owned.

## Dependency Gate

Before implementation, the gate was checked using the issue and dependency closure state. The required prerequisite issues are closed, so the gate is satisfied.

Required checks:

```powershell
git fetch origin --prune
gh issue view 89 --repo tannerpolley/ePC-SAFT --json number,title,state,body,url
git branch --show-current
git status --short
git rebase origin/main
```

Dependency status:

- Task C / issue #86: closed
- Task D / issue #87: closed

If a future recheck finds a missing dependency or a rebase conflict, stop with status `BLOCKED_DEPENDENCY_OR_REBASE`.

## Scope

- ideal/apparent mode
- activity-based mode
- fixed K mode
- fitted K mode support where schema exists
- diagnostics for residual blocks and derivative status
- ePC-SAFT activities/fugacities where requested

## Do Not Do

- do not implement MEA-specific speciation APIs
- do not fit reaction constants by default

## Validation

```powershell
uv run python run_pytest.py tests/api/test_reactive_staged_workflow_contract.py tests/native/test_reactive_speciation_implicit_sensitivity.py -q
```

```powershell
uv run python scripts/validate_project.py quick
```

