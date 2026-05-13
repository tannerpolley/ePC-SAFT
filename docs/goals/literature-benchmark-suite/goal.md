# Literature Benchmark Suite

## Objective

Inventory and then build generic literature benchmarks for package-level confidence while keeping `epcsaft` general-purpose and free of application-specific public APIs.

## Original Request

Task L: Literature benchmark suite.

Use the local live GoalBuddy board in Codex. Read GitHub issue #95, record the current issue scope in `docs/goals/literature-benchmark-suite/notes/issue_scope.md`, and prepare the goal on the assigned branch `codex/literature-benchmark-suite`.

## Intake Summary

- Input shape: existing_plan
- Audience: package maintainers and downstream users
- Authority: requested
- Proof type: artifact
- Completion proof: the goal board is prepared on the assigned branch, the issue scope is recorded, and the later `/goal` run can proceed with a live local board and bounded dependency gate.
- Likely misfire: start implementation too early, narrow the suite to a single downstream use case, or introduce application-specific public APIs.
- Blind spots considered: the benchmark tranche must remain generic, dependency gating must stay bounded, and any missing upstream solver/regression support should be treated as a blocker for implementation slices rather than papered over.
- Existing plan facts: branch bootstrap to `codex/literature-benchmark-suite`; GitHub issue #95; dependency list is `none`; local board must be used; watcher settings are bounded with `auto_start_after_gate: true`, `poll_interval_seconds: 120`, and `max_wait_minutes: 480`; validation is `uv run python scripts/validate_project.py docs` plus the listed regression tests; no finite difference; no application-specific public APIs; no local `.worktrees/`; focused PR must include `Closes #95`.

## Goal Kind

`existing_plan`

## Current Tranche

Prepare the board for a continuous execution run that first inventories the literature benchmark surface, then advances through the first safe generic benchmark slice once the dependency gate is confirmed open.

## Non-Negotiable Constraints

- Keep `epcsaft` general-purpose.
- Do not add application-specific public APIs.
- No finite difference.
- Use explicit derivative backends only.
- Keep the task board and receipts honest about blockers, gaps, and scope.
- Do not create local `.worktrees/`.

## Scope Anchors

- MEA simple workflow benchmark
- MDEA ePC-SAFT benchmark
- Figiel 2025 SSM+DS Born benchmark
- Held 2014 revised ePC-SAFT benchmark
- non-electrolyte LLE benchmark
- Ascani 2022 electrolyte LLE benchmark
- Ascani 2023 reactive LLE benchmark
- Khudaida salting-out LLE benchmark
- Hubach/Yu lithium-related equilibrium benchmark

## Stop Rule

Stop only when a final audit proves the full original outcome is complete.

Do not stop after planning, discovery, or Judge selection if a safe Worker task can be activated later in `/goal`.

## Canonical Board

Machine truth lives at:

`docs/goals/literature-benchmark-suite/state.yaml`

If this charter and `state.yaml` disagree, `state.yaml` wins for task status, active task, receipts, verification freshness, and completion truth.

## Run Command

```text
/goal Follow docs/goals/literature-benchmark-suite/goal.md.
```

