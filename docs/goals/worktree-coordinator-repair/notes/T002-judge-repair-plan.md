# T002 Judge Repair Plan

## Decision

Proceed with local coordination repairs, then patch repo prompt policy on a docs-only branch.

## Local GoalBuddy Repairs Allowed

Allowed local repair targets:

- `C:\Users\Tanner\.codex\worktrees\4a52\ePC-SAFT\docs\goals\derivative-backend-coverage-hard-gate\**`
- `C:\Users\Tanner\.codex\worktrees\964e\ePC-SAFT\docs\goals\explicit-cppad-parameter-derivatives\**`
- `C:\Users\Tanner\.codex\worktrees\b752\ePC-SAFT\docs\goals\generic-implicit-sensitivity-framework\**`
- `C:\Users\Tanner\.codex\worktrees\6443\ePC-SAFT\docs\goals\reaction-constant-conventions\**`
- `C:\Users\Tanner\.codex\worktrees\e188\ePC-SAFT\docs\goals\generic-target-row-schema\**`

Required local behavior:

- auxiliary dependency files under `notes/`
- watcher mode bounded
- `auto_start_after_gate: true`
- `poll_interval_seconds: 120`
- `max_wait_minutes: 480`
- no `PREPARED_READY` terminal behavior after implementation starts

Do not edit package source in those worktrees.

## Task B Scope Decision

Task B cannot remain at "recorded limitation, no PR." Its local board must say the next PM action is either:

- `PR_OPENED`, if the existing implementation can be opened with truthful coverage/capability language, or
- `BLOCKED_SCOPE_GAP`, if it cannot truthfully represent `k_hb_ij`.

The Task B repair must explicitly preserve:

- `k_hb_ij` appears in coverage/capability.
- If full active-association derivatives are unavailable, classify it as `blocker_requires_implicit_association_sensitivity`.
- Task C owns implicit association-site-fraction sensitivity.

## Repo Policy Decision

The main checkout still has stale prompt-policy content. A docs-only branch is required before repo prompt-policy edits:

```text
codex/worktree-coordinator-repair
```

Allowed repo files:

- `docs/roadmaps/agent_prompts/index.yaml`
- `docs/roadmaps/agent_prompts/*.md`
- `docs/roadmaps/watcher_templates/dependency_gate_template.yaml`
- `docs/roadmaps/watcher_templates/bounded_dependency_watcher.ps1`
- `docs/roadmaps/agent_dependency_plan.md`
- `docs/roadmaps/worktree_repair_status.md`
- policy docs under `docs/roadmaps/`

Do not touch package source code.

## PR #99 Decision

PR #99 is open and mergeable with green required checks. Its body and tests truthfully state that molality/apparent convention routes are defined but unsupported in the native residual and raise `backend_unavailable`. No application-specific API or finite-difference route was found in the inspected diff surface. Recommend merge after Tanner/normal repo policy approval; do not auto-merge from this coordinator run.
