# Ipopt Exact Hessian

## Objective

Implement exact Hessians everywhere currently feasible in the native Ipopt equilibrium stack, with explicit priority on making electrolyte LLE run successfully in exact-Hessian mode and exposing exact-Hessian usage as the default route behavior when a verified provider exists. If electrolyte-LLE exact Hessians require deeper thermodynamic derivatives, implementing that full derivative tier is part of this original goal rather than a separate scope increase.

## Original Request

"/goal

Fully implement exact hessians everywhere possible with the explicit goal of getting it to work with electrolyte LLE and we can see how it improved. It has to by default expecDoctor"

## Intake Summary

- Input shape: `specific`
- Audience: `ePC-SAFT` package maintainers and downstream users of generic equilibrium and reactive equilibrium APIs.
- Authority: `requested`
- Proof type: `test`
- Completion proof: electrolyte LLE succeeds through the public/native Ipopt route in exact-Hessian mode with diagnostics proving exact Hessian use, the default Hessian selection prefers exact on covered routes without silently downgrading explicit requests, additional feasible current-public routes gain exact Hessian coverage where implemented, focused validation passes, cleanup passes, and tracked changes are committed locally.
- Likely misfire: adapter-level exact-Hessian plumbing or diagnostics improve while the real electrolyte-LLE route still rejects exact Hessians, or exact mode works only when manually forced while default behavior still hides limited-memory selection truth.
- Existing plan facts:
  - `docs/plans/ipopt_improvement_plan.md` ranks exact Hessian support as a score-`9.5` item.
  - The adapter already supports exact-Hessian callbacks, diagnostics, and loud rejection when an NLP does not provide a Hessian.
  - Real public route-family tests still expect exact-Hessian rejection for electrolyte LLE, stability, neutral/reactive LLE, reactive electrolyte phase-equilibrium, and bubble/dew routes until providers are implemented.
  - Electrolyte LLE is the must-win route for this goal; broader coverage should expand only where the current route formulation makes exact Hessians genuinely implementable.
  - The original goal already implied adding the full higher-order derivative tier required by electrolyte-LLE exact Hessians if that is what the route needs to work.

## Goal Kind

`specific`

## Current Tranche

Complete the exact-Hessian tranche in safe vertical slices:

1. finish the GoalBuddy board and preserve the verified gap map;
2. implement the full higher-order thermodynamic derivative tier required by the electrolyte-LLE Lagrangian Hessian;
3. wire exact-Hessian provider coverage into electrolyte LLE and make default Hessian selection choose exact when coverage exists;
4. extend exact-Hessian providers to other currently feasible public Ipopt routes unlocked by the new derivative tier without weakening loud-failure semantics on uncovered routes;
5. validate native/public evidence, update capability/reporting truth, and close with a local commit.

## Non-Negotiable Constraints

- Electrolyte LLE exact-Hessian support is the primary success condition.
- Keep the package generic; no case-study-specific solver branches.
- Preserve loud failure for `hessian_mode="exact"` on routes that still lack verified providers.
- Preserve explicit `limited-memory` selection and accurate diagnostics on every solve.
- Do not use approximate derivative substitutes.
- If a new derivative tier is added, route providers must consume it end-to-end instead of leaving the new surface dark.
- Do not overclaim route capability in `epcsaft.capabilities()` or docs.
- Do not rerun the known freezing broad `tests/equilibrium -q` sweep.
- Keep `_core` rebuild coordination on the main thread only.
- Do not touch unrelated dirty tracked files already present in the worktree unless the user asks.

## Stop Rule

Stop only after a final audit shows whether the full derivative tier landed, whether electrolyte LLE exact Hessians now work, which additional routes gained exact-Hessian coverage, what remains impossible or unimplemented, and the exact verification/cleanup evidence behind those claims.

If electrolyte LLE exact Hessian still cannot be completed in this tranche, record the exact remaining derivative or route-ownership gap and keep the goal active rather than polishing the closeout.

## Canonical Board

Machine truth lives at:

`docs/goals/ipopt-exact-hessian/state.yaml`

If this charter and `state.yaml` disagree, `state.yaml` wins for task status, active task, receipts, verification freshness, and completion truth.

## Run Command

```text
/goal Follow docs/goals/ipopt-exact-hessian/goal.md.
```
