# Native LLE and Reactive Production Solvers

## Objective

Execute GitHub issues #116 and #117 as one dependency-gated GoalBuddy tranche: first replace the electrolyte LLE accepted production route with a Ceres trust-region residual solver using production Jacobians, then replace staged reactive phase equilibrium with one coupled native reactive phase-equilibrium residual solve.

## Original Request

Read issue #116 and #117 on GitHub and set up a GoalBuddy board for all phases to follow issue #116 and #117 explicitly in all of its steps.

## Intake Summary

- Input shape: `existing_plan`
- Audience: package maintainers and downstream ePC-SAFT consumers
- Authority: `requested`
- Proof type: `artifact`
- Completion proof: the board contains explicit, ordered tasks for every issue #116 and #117 stage, with #117 gated behind #116 production-solver completion, issue-specific validation commands, and a final audit against both definitions of done.
- Likely misfire: the board could drift into inventory, diagnostics, staged routes, synthetic fixtures, or capability-label closure instead of requiring production native Ceres solves with real residuals, Jacobians, benchmarks, and public generic APIs.
- Blind spots considered: current branch may not fast-forward to `origin/main`; dependency PR ancestry must be proven; #117 cannot begin until #116 is production-complete or implemented first in the same branch; GoalBuddy files must avoid banned exact backend-policy tokens from the roadmap.
- Existing plan facts: issue #116 is the native distributed-ion electrolyte LLE production solver plan; issue #117 is the native coupled reactive phase-equilibrium plan; both require Ceres trust-region residual solving on transformed variables, production Jacobians, generic APIs, benchmark proof, route guards, and full validation.

## Goal Kind

`existing_plan`

## Current Tranche

Prepare and then execute the full issue #116 and #117 implementation sequence. The board must start with dependency and intake proof, make no source edits before the issue-required intake notes exist, implement #116 through its production solver and benchmark validation, audit #116, then proceed into #117 only after the dependency gate is satisfied.

## Non-Negotiable Constraints

- Follow issue #116 and #117 bodies as controlling plans.
- Read `docs/roadmaps/FULL_ROADMAP.md` before execution and apply its completion standard.
- Keep public APIs generic and avoid application-specific MEA, lithium extraction, absorber, selectivity, efficiency, or distribution-coefficient APIs.
- Accepted production results must come from mature native solver backends with analytic, CppAD, analytic-implicit, CppAD-implicit, or Ceres-with-CppAD Jacobian provenance.
- Do not close with inventory, manifests, schema-only support, diagnostic-only routes, staged-only routes, synthetic-only fixtures, mocked payloads, documented limitations, capability labels, hand-coded simplex/Powell/Nelder-Mead production solvers, or derivative-approximation Jacobians.
- Preserve issue #116 before #117: #117 may use #116 only as an initialization/subcomponent source after #116 production behavior is implemented and verified.
- Keep package behavior upstream and generic; downstream repos own downstream metrics.
- Use the repo validation workflow from local `AGENTS.md` and the issue-specific commands.
- Avoid writing the exact banned backend-policy tokens named in `docs/roadmaps/FULL_ROADMAP.md` into committed goal, source, test, or documentation text.

## Stop Rule

Stop only when a final Judge or PM audit proves both issue definitions of done are complete. Do not stop after planning, discovery, partial diagnostics, one passing narrow test, or an honest limitation note.

## Canonical Board

Machine truth lives at:

`docs/goals/native-lle-reactive-production-solvers/state.yaml`

If this charter and `state.yaml` disagree, `state.yaml` wins for task status, active task, receipts, verification freshness, and completion truth.

## Run Command

```text
/goal Follow docs/goals/native-lle-reactive-production-solvers/goal.md.
```

## PM Loop

On every `/goal` continuation, read this charter, read `state.yaml`, refresh the GitHub issue bodies if needed, work only on the active task, write a compact receipt, update the board, and continue to the next largest safe work package until the full issue #116 and #117 outcome is complete.
