# T003 Judge Decision: Issue #116 Source-Edit Readiness

Date: 2026-05-15

Decision: approved_for_source_edits_on_issue_116

## Rationale

The dependency gate passed, issue-required intake notes exist, and the current route map is concrete enough to assign bounded source work. The existing implementation is not complete because accepted electrolyte LLE still depends on hand-coded optimizer paths and unavailable residual sensitivities.

## Approved Worker Sequence

1. T004 may edit only native equilibrium files, native equilibrium subdirectory files, native tests, and the issue #116 variable-transform note. It must establish solver ownership, explicit distributed-ion basis reporting, transformed feasible variables, and charge/material diagnostics.
2. T005 may add residual blocks, production Jacobian support, and the Ceres trust-region solve. It must fail if the accepted route can still report old hand-coded optimizer behavior as production.
3. T006 may wire the generic Python route and required distributed-ion plus salting-out benchmark tests. TPD and g-hat may seed, check stability, or support acceptance only.
4. T007 must run the issue #116 validation ladder and route-guard searches, then repair only narrow failures inside #116 scope.

## Required Stop Conditions

- Stop if a required production residual cannot be represented in one solved state.
- Stop if the accepted solve cannot expose a real Jacobian provenance.
- Stop if benchmark data would need fabrication.
- Stop if source changes need files outside the active Worker allowed_files.
- Stop if route guards leave an accepted production path on the old hand-coded optimizer labels.

## Verification Minimum For #116 Audit

Before T008 may unlock issue #117, receipts must show a Ceres accepted solve, available Jacobian diagnostics, material and charge accounting, distributed-ion and salting-out benchmark assertions, generic Python API proof, and the targeted validation commands from the board.
