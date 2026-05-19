# Ipopt Score 8 Plus

## Objective

Implement the remaining Ipopt improvement-plan items with score `>= 8.0` so the full score-`>= 8.0` set is covered in production code and verified through the package's native/public equilibrium routes.

This goal starts after the completed tranche-1 work on branch `ipopt` at commit `7095aeea`, which already covered the score `>= 9.0` items:

- per-iteration diagnostics
- real scaling diagnostics
- warm-start continuation state
- Hessian-mode plumbing and loud exact-Hessian rejection when no provider exists

## Original Request

"/goal Now implement all the ipopt changes from the document that have a score of 8 or higher"

## Intake Summary

- Input shape: `existing_plan`
- Audience: `ePC-SAFT` package maintainers and downstream consumers using generic equilibrium/speciation APIs.
- Authority: `requested`
- Proof type: `test`
- Completion proof: the remaining score-`8.x` items from `docs/plans/ipopt_improvement_plan.md` are implemented, exercised through focused route/public tests, broad validation passes, cleanup passes, and the tracked changes are committed locally.
- Likely misfire: only exposing more option names without route-owned behavior, or claiming score-`>= 8.0` completion while sparse Jacobians, linear-solver reporting, tolerance families, or initial-point strategy remain stubbed or route-limited without tests.
- Existing plan facts:
  - Tranche 1 is already complete on `ipopt`.
  - Remaining score `>= 8.0` items are initial-point strategy, linear solver selection/reporting, sparse Jacobian structures, and separate tolerance families.
  - Keep the package generic; no Khudaida-only solver behavior.
  - Preserve exact derivative/Jacobian requirements and honest capability boundaries.

## Goal Kind

`existing_plan`

## Current Tranche

Implement and verify the remaining score-`8.x` Ipopt features in safe vertical slices:

1. map the current adapter/route/test gaps;
2. add route-owned initial-point strategy behavior;
3. expose and verify linear-solver plus tolerance-family controls/diagnostics;
4. replace dense-by-default Jacobian structures where sparsity is known and verifiable;
5. run a final audit against the score-`>= 8.0` plan surface.

## Non-Negotiable Constraints

- Keep the package generic; no downstream-specific public API or behavior.
- Do not broaden `epcsaft.capabilities()` or docs beyond executable evidence.
- Do not use approximate derivative substitutes.
- Do not silently degrade an explicitly requested mode or control.
- Keep `_core` rebuild coordination on the main thread only.
- Preserve the completed `docs/goals/ipopt-tranche1/` board as historical truth; this goal supersedes it for the remaining score-`8.x` work only.

## Stop Rule

Stop only when a final audit maps current receipts, verification, cleanup, and git status back to the score-`>= 8.0` items in `docs/plans/ipopt_improvement_plan.md`.

If any score-`8.x` item cannot be completed in this goal, record the exact missing route/backend gap and leave the goal active rather than polishing the closeout.

## Canonical Board

Machine truth lives at:

`docs/goals/ipopt-score8plus/state.yaml`

If this charter and `state.yaml` disagree, `state.yaml` wins for task status, active task, receipts, verification freshness, and completion truth.

## Run Command

```text
/goal Follow docs/goals/ipopt-score8plus/goal.md.
```
