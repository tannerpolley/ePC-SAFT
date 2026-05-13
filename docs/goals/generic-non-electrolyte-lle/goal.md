# Generic non-electrolyte LLE benchmark and solver hardening

## Original Request

Run Goal Prep for task key H from `docs/roadmaps/agent_prompts/index.yaml`: Generic non-electrolyte LLE benchmark and solver hardening.

## Current Status

GATE_PASSED

Dependency C is merged into `origin/main`, and the assigned branch rebased cleanly. The board is prepared to continue with bounded execution defaults.

Watcher settings:

- `watcher_mode: bounded`
- `auto_start_after_gate: true`
- `poll_interval_seconds: 120`
- `max_wait_minutes: 480`

## Goal

Prepare the board for Task H and keep the work generic: prove ordinary two-liquid-phase splitting before any electrolyte layering, while keeping `epcsaft` general-purpose.

## Non-Negotiable Constraints

- Current branch must be `codex/generic-non-electrolyte-lle`.
- Do not create local `.worktrees/`.
- Do not add application-specific public APIs.
- No finite difference.
- Do not include electrolyte accounting in this issue.
- Do not force a poor benchmark to pass.
- Use only generic concepts and outputs.

## Dependency Gate

Verified before implementation work:

- Issue #86 is closed.
- PR #104 is merged.
- Commit `974e6a232025d9d305297590228c0d2131fbc4fe` is present in the repo history.
- `git rebase origin/main` completed cleanly on `codex/generic-non-electrolyte-lle`.

## Scope

- phase split
- fugacity equality
- stability checks or anti-trivial-solution strategy
- clear phase diagnostics
- simple literature or repo benchmark

## Out of Scope

- electrolyte accounting
- application-specific public APIs
- finite difference
- forcing an unsuitable benchmark to pass

## Start Signal

The board is ready for the next `/goal` run:

```text
/goal Follow docs/goals/generic-non-electrolyte-lle/goal.md.
```
