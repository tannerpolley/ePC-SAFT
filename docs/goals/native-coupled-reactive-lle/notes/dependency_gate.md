# Issue #117 Dependency Gate

Date: 2026-05-15

Issue: https://github.com/tannerpolley/ePC-SAFT/issues/117

## Shared Gate Result

- Branch: `#116-117`
- HEAD before source edits: `4f1efb046b4614e9a941ded5371349d9246924bd`
- `origin/main` after `git fetch origin --prune`: `f3ce3cc332bc06e8504ece122e478c8396bc27b8`
- Fast-forward gate: `git merge --ff-only origin/main` reported `Already up to date.`
- Worktree before source edits: clean by `git status --short`.
- Required dependency PRs are merged and ancestors of this branch:
  - PR #123, merge commit `5391b4011f46f6817081a81faec0f3c5f407580d`
  - PR #124, merge commit `b8cc04b095c0cd9d53d7de08c9315d4040e691ff`
  - PR #125, merge commit `acf6fc74c90734f2cbcab4f79dfa73a801297de6`

## Issue #117 Gate Decision

Issue #117 source edits are not open yet. The coupled reactive LLE production solver depends on issue #116 because reactive electrolyte LLE must be able to use the production electrolyte LLE solver as a subcomponent, seed, comparator, or acceptance guard. The current branch may implement both issues, but #116 must be completed and audited first.

## Current Blocking Evidence

- The current reactive phase workflow is staged: `solve_reactive_staged_equilibrium` solves chemical equilibrium first, then calls a phase route.
- `ReactivePhaseEquilibriumProblem.solve` delegates to `mixture.reactive_staged_equilibrium(...)`.
- Current diagnostics intentionally report `workflow = chemical_equilibrium_then_phase_equilibrium`, `reactive_workflow_class = staged`, `coupling_level = staged_not_full_simultaneous_nlp`, and `full_simultaneous_reactive_nlp = False`.
- Native chemical equilibrium from issue #115 exists as a homogeneous chemical-equilibrium route. It is not a coupled reactive phase-equilibrium solver.
- Issue #116 is not complete at this gate because electrolyte LLE residual sensitivities still report `not_available`, and accepted production electrolyte LLE is not yet a Ceres trust-region residual solve with a production Jacobian.

## Allowed Work Before #116 Audit

- Read-only route mapping for issue #117.
- Intake notes and fixture selection.
- Board-state updates and receipts.
- No reactive source, test, API, or benchmark implementation unless the board explicitly chooses a same-branch path that completes and audits #116 first.

## Unlock Condition

Issue #117 may start only after the issue #116 audit records that accepted electrolyte LLE production results come from a native Ceres trust-region residual solve with real Jacobian diagnostics, public generic Python route proof, distributed-ion benchmark proof, salting-out benchmark proof, and passing targeted validation.
