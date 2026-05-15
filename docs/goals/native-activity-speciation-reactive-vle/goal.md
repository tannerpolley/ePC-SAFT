# Native Activity Speciation Reactive VLE

## Objective

Prepare and execute the full issue #115 vertical implementation: native coupled activity speciation, implicit solved-state sensitivities, volatile-neutral pressure/VLE, generic Python API exposure, and a pressure/speciation benchmark that proves the production path end to end.

## Original Request

Read GitHub issue #115 and set up a GoalBuddy board for all phases to follow issue #115 explicitly in all steps. Use JetBrains MCP where smart and helpful.

## Intake Summary

- Input shape: `existing_plan`
- Audience: ePC-SAFT package maintainers and downstream MEA-Thermodynamics consumers blocked on generic activity-coupled speciation.
- Authority: `requested`
- Proof type: `test`
- Completion proof: issue #115 remains open until the native backend, derivative/sensitivity path, Python API surface, pytest suite, docs, and benchmark proof all pass together through generic package APIs.
- Likely misfire: completing inventory, schema, diagnostics, staged helpers, synthetic payloads, or an "already supported" assertion without the exact passing tests required by issue #115.
- Blind spots considered: banned-token hygiene in GoalBuddy files, generic package-scope API boundaries, downstream evidence from the issue comment, native-regression no-edit boundary, and the need to keep Python out of production solve loops.
- Existing plan facts: GitHub issue #115 defines Stage 0 through Stage 5, owned paths, banned implementation escape routes, required derivative mechanisms, validation commands, completion definition, and a downstream MEA-Thermodynamics blocker comment.

## Goal Kind

`existing_plan`

## Current Tranche

This tranche is the whole issue #115 vertical implementation. The board must preserve the issue order and completion contract:

1. Stage 0 intake gate and evidence note.
2. Stage 1 native coupled activity speciation.
3. Stage 2 implicit sensitivity for solved composition.
4. Stage 3 volatile-neutral VLE / partial-pressure path.
5. Stage 4 generic Python API exposure.
6. Stage 5 benchmark proof with repo-contained fixture data.
7. Required validation and final audit against issue #115.

## Non-Negotiable Constraints

- Follow issue #115 at every phase boundary and in every task receipt.
- Do not edit native regression code for this issue.
- Do not add downstream-application-specific public APIs.
- Do not satisfy required production behavior with Python-owned production optimizer loops, Python-owned production nonlinear solver loops, numeric-differencing Jacobians, derivative-free production minimizers, staged-only production equilibrium, mocked equilibrium payloads, or application-specific public APIs.
- Production derivative support must be analytic, CppAD, analytic implicit sensitivity, or CppAD implicit sensitivity.
- For solved internal states, use implicit sensitivity rather than differentiating through iterative loops.
- Activity/fugacity terms must be evaluated inside the nonlinear residual, not from stale external activity vectors.
- Charged species remain liquid-only unless the problem explicitly defines a charged vapor model.
- Public API concepts must stay generic, such as `ReactiveEquilibriumProblem`, `equilibrium(...)`, and `solve_reactive_speciation(...)`.
- The old missing-backend status token and old numeric-differencing token family must not appear as contiguous literals in committed repository text, including GoalBuddy files.
- Use post-PR-113 build and validation commands from issue #115.
- Treat `state.yaml` as board truth if it differs from this charter.

## Stop Rule

Stop only when a final Judge or PM audit proves the full issue #115 completion definition is satisfied. Do not stop after discovery, a working scaffold, a staged route, one passing focused test, or a draft PR that honestly reports incompleteness.

## Slice Sizing

Each Worker task should be the largest safe useful phase slice that preserves the issue boundary, allowed files, and validation path. Judge tasks occur at Stage 0 readiness, high-risk native/derivative boundaries, benchmark proof, and final completion.

## Canonical Board

Machine truth lives at:

`docs/goals/native-activity-speciation-reactive-vle/state.yaml`

## Run Command

```text
/goal Follow docs/goals/native-activity-speciation-reactive-vle/goal.md.
```

## PM Loop

On every `/goal` continuation:

1. Read this charter.
2. Read `state.yaml`.
3. Read GitHub issue #115 and its downstream blocker comment before source edits.
4. Work only on the active board task.
5. Use JetBrains MCP for semantic navigation when tracing public API, native bindings, call hierarchy, or whether code is truly unused.
6. Write receipts that map the result back to the exact issue #115 stage and completion contract.
7. Advance to the next required stage until the final audit proves the full issue is complete.
