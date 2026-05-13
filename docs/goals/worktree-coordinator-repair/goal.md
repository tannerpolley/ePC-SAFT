# Worktree Coordinator Repair

## Original Request

Inspect and repair all active ePC-SAFT roadmap Codex worktrees using GoalBuddy, but do not start `/goal` execution during prep.

## Interpreted Outcome

The next `/goal` run has a clear, constrained board for auditing every active roadmap worktree under `C:\Users\Tanner\.codex\worktrees`, repairing local GoalBuddy automation files, patching roadmap prompt policy only when needed, validating changed docs, and reporting exact status without touching package source code or creating repo-local `.worktrees/`.

## Goal Mode

- Kind: recovery
- Input shape: existing_plan
- Audience: Tanner and the next GoalBuddy PM thread
- Proof type: artifact and review

## Non-Negotiable Constraints

- Do not start `/goal` execution during Goal Prep.
- Do not modify package source code.
- Do not create local `.worktrees/`.
- Do not use PR #56 as a base.
- Do not add application-specific public APIs.
- Do not invent transcript content when Codex transcript/log files are unavailable.
- If only inspecting and writing local coordination notes, stay on `main`.
- If editing repo roadmap or prompt docs with intent to commit, use `codex/worktree-coordinator-repair`.
- Use the local live GoalBuddy board in Codex.
- Do not ask repeated confirmation questions; use the defaults from the prompt.

## Required Inputs For Execution

- Attached handoff: `C:\Users\Tanner\Downloads\epcsaft_worktree_repair_coordinator_bundle.zip`
- Worktree root: `C:\Users\Tanner\.codex\worktrees`
- Repo files, if present:
  - `docs/roadmaps/general_reactive_electrolyte_equilibrium_readiness.md`
  - `docs/roadmaps/agent_dependency_plan.md`
  - `docs/roadmaps/agent_prompts/index.yaml`
  - `docs/roadmaps/watcher_templates/dependency_gate_template.yaml`
  - `docs/roadmaps/watcher_templates/bounded_dependency_watcher.ps1`

## Existing Plan Facts To Preserve

- Inventory every worktree under `C:\Users\Tanner\.codex\worktrees`.
- Identify task A-M for each active worktree.
- Read GoalBuddy files, `state.yaml`, `goal.md`, notes, dependency gates, watch scripts, git status, branch, and PR state.
- Read Codex transcript/log files only when available.
- Repair local GoalBuddy automation files where needed:
  - move `dependency_gate.yaml` into `notes/`
  - create or update `notes/watch_dependency.ps1`
  - set `watcher_mode: bounded`
  - set `auto_start_after_gate: true`
  - set `poll_interval_seconds: 120`
  - set `max_wait_minutes: 480`
  - remove `PREPARED_READY` as terminal behavior after implementation has started
- Repair Task B scope handling:
  - no vague `k_hb_ij` recorded limitation
  - require `PR_OPENED` or `BLOCKED_SCOPE_GAP`
  - include `k_hb_ij` in coverage/capability as blocker or requiring implicit association sensitivity if not implemented
  - identify Task C as owner of implicit association sensitivity
- Review Task D PR #99 for truthfulness; do not auto-merge unless safe and authorized by available tool or policy.
- Patch repo prompt policy if needed:
  - `index.yaml` A-M `auto_start_after_gate: true`
  - watcher mode bounded
  - dependency gates under `notes/`
  - scope contract policy
  - no silent narrowing
- If repo files changed, validate docs with `uv run python scripts/validate_project.py docs`.
- Open a docs-only PR if repo prompt policy changed and normal repo/GitHub policy allows it.

## Likely Misfire To Avoid

GoalBuddy could succeed at a narrow local-file cleanup while failing to map all active worktrees, truthfully classify B/C/D/E, inspect PR #99, or produce the required `docs/roadmaps/worktree_repair_status.md` report.

## Completion Proof

Completion requires `docs/roadmaps/worktree_repair_status.md` plus task receipts showing the worktree inventory, A-M mapping, B/C/D/E status, PR #99 status, repaired local goal files, patched repo files if any, validation results when repo files changed, and PR URL if opened.

## Starter Command

```text
/goal Follow docs/goals/worktree-coordinator-repair/goal.md
```
