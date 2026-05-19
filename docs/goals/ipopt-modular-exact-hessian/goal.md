# Modular Ipopt Exact Hessian

## Objective

Implement a modular exact-Hessian architecture for the native Ipopt equilibrium stack. The work keeps the existing `NlpProblem` Hessian callback as the only Ipopt-facing interface, moves route Hessian algebra into shared second-order carriers plus one Lagrangian Hessian assembler, and makes exact Hessians the default accepted path across native Ipopt equilibrium routes.

## Original Request

"PLEASE IMPLEMENT THIS PLAN: Modular Ipopt Exact-Hessian Refactor Plan"

## Intake Summary

- Input shape: `existing_plan`
- Audience: `ePC-SAFT` package maintainers and downstream users of generic equilibrium and reactive equilibrium APIs.
- Authority: `requested`
- Proof type: `test`
- Completion proof: all covered native Ipopt route families expose and use exact Hessians by default, representative route fixtures are accepted, Ipopt diagnostics prove `eval_h` execution, selected proof runs capture Ipopt logs, focused validation passes, cleanup passes, and tracked changes are locally committed.
- Likely misfire: adding another route-specific Hessian implementation while leaving the shared derivative architecture incomplete, or setting provider flags without proving Ipopt consumed exact Hessian callbacks on accepted route fixtures.

## Coverage Target

Coverage includes every native Ipopt equilibrium backend/internal path used by public or route-helper surfaces:

- reactive speciation standard-state routes;
- neutral TP flash and neutral LLE flash;
- neutral bubble and dew routes;
- neutral and electrolyte stability routes;
- electrolyte LLE and electrolyte bubble-pressure routes;
- reactive LLE, reactive electrolyte LLE, and reactive two-phase internals.

## Non-Negotiable Constraints

- Keep the package generic; do not add case-study-specific solver branches.
- Keep public Python APIs stable. Additive diagnostics are allowed; public derivative objects are not part of this goal.
- Exact Hessians mean the actual Ipopt Lagrangian Hessian, including objective, nonlinear constraint, residual, multiplier, and transform terms needed by the route.
- Route classes may declare variables, constraints, residuals, and transforms, but lower-triangle extraction, multiplier weighting, sparsity handling, and Lagrangian Hessian assembly must live in shared infrastructure.
- Existing exact-Hessian routes must migrate to the shared assembler instead of remaining permanent one-off implementations.
- Do not use surrogate derivative checks as proof or implementation.
- Preserve explicit `limited-memory` Hessian mode as a truthful opt-out.
- Explicit `hessian_mode="exact"` may fail for invalid inputs or a build without Ipopt, not because a route lacks a provider.
- Do not run the known-freezing broad `tests/equilibrium -q` sweep.
- Keep `_core` rebuild coordination on the main thread only.
- Do not touch unrelated dirty tracked files unless they directly block this goal.

## Acceptance Gates

1. Architecture gate: route files no longer own bespoke lower-triangle Hessian extraction or multiplier-weighted Lagrangian assembly.
2. Provider gate: covered routes report `exact_hessian_available == true`, `hessian_approximation == "exact"`, `hessian_backend != "limited-memory"`, and `eval_h_calls > 0` for default `auto` and explicit `exact`.
3. Accepted-fixture gate: each route family has at least one accepted solve fixture using default exact Hessians.
4. Hydrocarbon proof gate: the first neutral VLE proof uses the methane/ethane/propane workbook benchmark from `workbooks/PC-SAFT Calculations - Hydrocarbon Basis.xlsm`; the workbook saves a bubble-pressure VLE point at `T=233.15 K`, `P=1,276,369.4735856401 Pa`, liquid `x=[0.1, 0.3, 0.6]`, vapor `y≈[0.7246628928, 0.2029319137, 0.07240519344]`, and phase-specific liquid/vapor density roots.
5. Electrolyte proof gate: electrolyte LLE remains an accepted exact-Hessian proof route with diagnostics.
6. Log gate: representative VLE and electrolyte LLE proof commands capture Ipopt logs while tests assert structured diagnostics and iteration history.
7. Opt-out gate: explicit `limited-memory` still runs where supported and reports `eval_h_calls == 0`.

## Stop Rule

Stop only after a final audit states which route families satisfy the accepted default-exact gate, which proof commands ran, whether any required derivative tier remains incomplete, and whether the goal can truthfully be marked done.

If any lower derivative tier cannot support exact Hessians for a covered route, keep the goal active and record the exact owner gap rather than substituting limited-memory mode or a compatibility flag.

## Canonical Board

Machine truth lives at:

`docs/goals/ipopt-modular-exact-hessian/state.yaml`

If this charter and `state.yaml` disagree, `state.yaml` wins for task status, active task, receipts, verification freshness, and completion truth.

## Run Command

```text
/goal Follow docs/goals/ipopt-modular-exact-hessian/goal.md.
```
