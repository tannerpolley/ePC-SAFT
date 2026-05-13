# General Reactive/Electrolyte Roadmap Setup For ePC-SAFT Agents

## Objective

Prepare and execute a docs/issues/project/branch setup tranche for the next ePC-SAFT general reactive/electrolyte roadmap, using the provided setup bundle as the source of truth and without implementing package source code.

## Original Request

General reactive/electrolyte roadmap setup for ePC-SAFT agents. Use local live GoalBuddy board. Do not start `/goal` during prep. Do not edit package source files during Goal Prep. Create the goal control files, then print the exact `/goal Follow docs/goals/<slug>/goal.md.` command once and stop.

## Intake Summary

- Input shape: `existing_plan`
- Audience: ePC-SAFT Codex agents and the repository maintainer
- Authority: `requested`
- Proof type: `artifact`
- Completion proof: the setup PR contains the requested setup files from `C:\Users\Tanner\Downloads\epcsaft_agent_setup_bundle_branch_bootstrap.zip`, validates with `uv run python scripts/validate_project.py docs`, modifies no package source files, and either creates the requested GitHub issues/project/branches where available or records manual follow-up requirements.
- Likely misfire: treating the roadmap setup as permission to implement reactive/electrolyte package source changes, creating implementation branches before the setup PR is merged, or improvising GitHub Project behavior when tooling is unavailable.
- Blind spots considered: GitHub Project tool availability, `gh` authentication, setup bundle path/content validity, whether merge permissions are available, and keeping local work clean without repo-local `.worktrees/`.
- Existing plan facts: use branch `codex/general-reactive-electrolyte-roadmap-setup`; run branch bootstrap before writing setup files; use local live GoalBuddy board; create only `docs/goals/<slug>/goal.md`, `docs/goals/<slug>/state.yaml`, and `docs/goals/<slug>/notes/` during Goal Prep; during `/goal`, create/update root setup prompts, branch caveat docs, roadmap docs, agent prompts, issue drafts, watcher templates, issues A-M if available, the GitHub Project if available, setup PR, validation, and only post-merge placeholder branches.

## Goal Kind

`existing_plan`

## Current Tranche

Validate and operationalize the provided setup plan, install the setup-bundle roadmap artifacts into the repository, create supported GitHub tracking artifacts, open or merge the setup PR according to available permissions, and stop before implementation-source work. If the setup PR cannot be merged, do not create implementation branches.

## Non-Negotiable Constraints

- Current branch must be `codex/general-reactive-electrolyte-roadmap-setup`.
- Use `C:\Users\Tanner\Downloads\epcsaft_agent_setup_bundle_branch_bootstrap.zip` as the source of truth for setup content.
- Do not implement package source code.
- Do not edit package source files during Goal Prep.
- Do not create local `.worktrees/`.
- Do not ask the user to choose among options or ask repeated confirmation questions.
- If a missing value can be defaulted safely, use the safest conservative default and record it in `notes/assumptions.md`.
- If a missing value prevents safe work, stop with status `BLOCKED_MISSING_INPUT`.
- Do not create placeholder implementation branches until the setup PR is merged into `origin/main`.
- Do not close parent roadmap issues unless explicitly instructed.

## Stop Rule

Stop only when a final audit proves the full original outcome is complete.

Do not stop after planning, discovery, or Judge selection if a safe Worker or PM task can advance the setup goal.

Do not stop after creating docs alone if the plan still has safe GitHub issue, PR, validation, or branch-gating tasks available.

Do not stop because a task needs unavailable GitHub Project tooling or merge permission. Mark that exact slice blocked with a receipt, complete all safe local setup work, and report the manual follow-up.

## Canonical Board

Machine truth lives at:

`docs/goals/general-reactive-electrolyte-roadmap-setup/state.yaml`

If this charter and `state.yaml` disagree, `state.yaml` wins for task status, active task, receipts, verification freshness, and completion truth.

## Run Command

```text
/goal Follow docs/goals/general-reactive-electrolyte-roadmap-setup/goal.md.
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
