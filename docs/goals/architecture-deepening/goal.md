# Architecture Deepening

## Objective

Implement the six architecture-deepening candidates identified by the improve-codebase-architecture pass, carrying each one to a maintainable code/test/doc finish without broadening public capability claims beyond evidence.

## Original Request

`/goal` followed by `Implement all of them`, referring to the six architecture candidates from the prior codebase architecture review.

## Intake Summary

- Input shape: `existing_plan`
- Audience: ePC-SAFT package maintainers and downstream package consumers.
- Authority: `requested`
- Proof type: `test`
- Completion proof: all six candidates have corresponding code/test/doc updates or explicit blocked receipts, focused tests pass for each touched boundary, the GoalBuddy final audit maps receipts back to the original six candidates, cleanup passes, and tracked repo changes are committed locally.
- Likely misfire: only producing a board, wrappers, or documentation while leaving the duplicated behavior paths and evidence boundaries essentially unchanged.
- Blind spots considered: native `_core` rebuild coordination, public API compatibility, downstream capability claims, broad validation cost, and keeping reactive/electrolyte claims evidence-backed.
- Existing plan facts: preserve and implement the six candidates: equilibrium request/result normalization; native route variableization and residual core; Ipopt solve metadata and diagnostics; regression evidence and target-family summaries; ParameterSet payload boundary; capability evidence and validation lane registry.

## Goal Kind

`existing_plan`

## Current Tranche

Continuously execute safe, verified implementation slices until all six candidates are complete or a specific candidate has a durable blocker receipt. Start with the lowest-risk Python architecture boundaries, then move through validation/capability evidence, then native Ipopt diagnostics and native route variableization. Main thread owns `_core` rebuilds and cross-slice integration.

## Non-Negotiable Constraints

- Preserve current public APIs unless a local test-backed simplification proves the migration path.
- Do not add Codex-specific tracked public files.
- Do not broaden `capabilities()` claims without executable evidence.
- Coordinate native rebuilds through the main thread; do not clean or delete `_core` while agents may import it.
- Use focused validation first, then `uv run python scripts/dev/validate_project.py quick` before closeout if the tree is healthy enough.
- Keep docs concise and current; remove obsolete workflow text if a workflow is replaced.

## Stop Rule

Stop only when a final audit proves the full original outcome is complete.

Do not stop after planning, discovery, or Judge selection if a safe Worker task can be activated.

Do not stop after a single verified Worker package when the broader owner outcome still has safe local follow-up work. Advance the board to the next highest-leverage safe Worker package and continue unless a phase, risk, rejected-verification, ambiguity, or final-completion review is due.

Do not create one Worker/Judge pair per repeated file, table, route, or helper. Put repeated same-shape work into one Worker package and review the package as a whole.

## Slice Sizing

Safe means bounded, explicit, verified, and reversible. It does not mean tiny.

A good task is the largest safe useful slice.

A Worker should finish the whole assigned slice. A Judge should judge the whole slice. The PM should reorient the board when tasks are safe but not moving the outcome.

## Canonical Board

Machine truth lives at:

`docs/goals/architecture-deepening/state.yaml`

If this charter and `state.yaml` disagree, `state.yaml` wins for task status, active task, receipts, verification freshness, and completion truth.

## Run Command

```text
/goal Follow docs/goals/architecture-deepening/goal.md.
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
9. If safe local work remains, choose the next largest reversible Worker package and continue unless blocked.
10. Review at phase, risk, rejected-verification, ambiguity, or final-completion boundaries; do not review every small Worker by habit.
11. Finish only with a Judge/PM audit receipt that maps receipts and verification back to the original user outcome and records `full_outcome_complete: true`.
