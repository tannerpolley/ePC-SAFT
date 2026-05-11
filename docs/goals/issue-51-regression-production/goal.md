# Issue #51 Regression Production Readiness

## Objective

Fully implement GitHub issue #51: finish canonical parameter-schema polish and make bounded mixed pressure/speciation reactive-electrolyte regression production-grade, with real tests, benchmark evidence, docs, capabilities updates, PR reporting, and no public `bounded_incomplete` status.

## Original Request

Look at the new issue posted and prepare for its full implementation.

## Intake Summary

- Input shape: `existing_plan`
- Audience: ePC-SAFT package users, downstream manuscript/modeling agents, and downstream MEA-style pressure/speciation regression workflows.
- Authority: `requested`
- Proof type: `test | metric | artifact`
- Completion proof: issue #51 is implemented on a branch, verified locally with required tests/docs/lint/build/benchmarks, PR is created with required reporting sections and benchmark evidence, merged only when checks pass, and no public package payload emits `bounded_incomplete`.
- Likely misfire: treating this as documentation or relabeling only, while the mixed pressure/speciation bounded regression path remains placeholder-like, untested, or still objectively incomplete.
- Blind spots considered: molar-mass unit safety can silently corrupt electrolyte workflows; benchmark defaults can become too slow; the optimizer must not pretend to be a global mathematical optimizer; downstream `bounded_incomplete` may originate outside this package but still requires a clear package status contract.
- Existing plan facts: preserve issue #51 phases 1-6, non-goals, validation commands, benchmark requirements, PR reporting sections, and definition of done.

## Goal Kind

`existing_plan`

## Current Tranche

Complete issue #51 end to end. The safe first slice is read-only evidence mapping of current `main` against every issue requirement, followed by Judge validation of implementation order and disjoint write scopes. After that, implement and verify each slice until the full issue definition of done is met.

## Non-Negotiable Constraints

- Do not modify or remove `analyses/2014_held/...`.
- Do not undo PR #50 architecture changes.
- Do not remove legacy dictionary parameter support or existing top-level imports.
- Do not silently relabel incomplete regression as production.
- Do not make IPOPT automatic.
- Do not require private downstream MEA-Thermodynamics files for package tests.
- Do not migrate ePC-SAFT electrolyte/Born handling into FeOs.
- Preserve native-first thermodynamic behavior; Python may orchestrate/batch but should not fake thermodynamic solver completion.
- Use real package solver calls for mixed pressure/speciation regression tests.
- Keep benchmark and docs claims precise: “bounded multi-row mixed pressure/speciation least-squares regression,” not global optimization in the mathematical sense.

## Required Issue #51 Scope

1. Fix canonical molar-mass units and safeguards.
2. Add `create_parameter_template(..., schema="canonical")` while preserving `schema="legacy"`.
3. Add `epcsaft.eos.Mixture` and `epcsaft.eos.State`.
4. Make bounded mixed pressure/speciation reactive-electrolyte regression production-grade, including precise statuses and no `bounded_incomplete`.
5. Add a 35-row pressure/speciation package benchmark or downstream validation evidence.
6. Update docs, capabilities, tests, benchmark reporting, and PR reporting standards.

## Stop Rule

Stop only when a final audit proves the full original issue #51 outcome is complete.

Do not stop after planning, discovery, or Judge selection if safe implementation tasks remain.

Do not stop after one verified slice if any issue #51 definition-of-done item remains safe to implement locally.

If a slice requires downstream-private files or external input, mark only that slice blocked, create the smallest package-owned surrogate or status-contract task that still advances the issue, and keep going.

## Canonical Board

Machine truth lives at:

`docs/goals/issue-51-regression-production/state.yaml`

If this charter and `state.yaml` disagree, `state.yaml` wins for task status, active task, receipts, verification freshness, and completion truth.

## Run Command

```text
/goal Follow docs/goals/issue-51-regression-production/goal.md.
```

## PM Loop

On every `/goal` continuation:

1. Read this charter.
2. Read `state.yaml`.
3. Run the bundled GoalBuddy update checker when available and mention a newer version without blocking.
4. Re-check issue #51 and preserve any new comments or edits.
5. Work only on the active board task.
6. Assign Scout, Judge, Worker, or PM according to the task.
7. Write a compact task receipt.
8. Update the board.
9. If Judge selected a safe Worker task with `allowed_files`, `verify`, and `stop_if`, activate it and continue unless blocked.
10. Treat a slice audit as a checkpoint, not completion, unless it explicitly proves every issue #51 definition-of-done item.
11. Finish only with a Judge/PM audit receipt that maps receipts and verification back to issue #51 and records `full_outcome_complete: true`.
