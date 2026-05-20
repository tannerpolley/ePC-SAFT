# Ipopt Tranche 1

## Objective

Implement the score `>= 9.0` Ipopt improvement tranche for native equilibrium NLP routes: bounded per-iteration diagnostics, real scaling diagnostics, warm-start continuation state, and exact Hessian support with explicit diagnostics and loud rejection when an exact Hessian is requested but unavailable.

## Original Request

"PLEASE IMPLEMENT THIS PLAN: Ipopt Tranche 1 Shipshape Prep"

## Intake Summary

- Input shape: `existing_plan`
- Audience: ePC-SAFT package maintainers and downstream consumers using generic native equilibrium APIs.
- Authority: `requested`
- Proof type: `test`
- Completion proof: focused native/Python tests prove diagnostics, scaling, warm-start validation, and Hessian-mode behavior; Ipopt-enabled build/doctor proof is run or a concrete blocker is reported; cleanup passes; tracked tranche changes are committed locally.
- Likely misfire: adding diagnostic fields or option names that look complete while routes still silently use limited-memory Hessians, unvalidated scaling, or incompatible continuation state.
- Blind spots considered: exact Ipopt Hessians require Lagrangian Hessians, not just objective Hessians; fast-profile builds do not prove Ipopt; route capability claims must remain generic and evidence-backed; generated Graphify output should not bloat the repo.
- Existing plan facts: implement per-iteration diagnostics first, then scaling, warm starts, and exact Hessian support; preserve limited-memory as an explicit selectable mode; default `hessian_mode="auto"` may choose limited-memory when exact route coverage is unavailable, but `hessian_mode="exact"` must not silently fall back.

## Goal Kind

`existing_plan`

## Current Tranche

Execute successive safe verified slices until the tranche behavior exists through the native adapter and public equilibrium APIs. Use test-first vertical slices where practical: adapter diagnostics/scaling first, warm-start state second, public option/result propagation third, route Hessian-mode gates last.

## Non-Negotiable Constraints

- Keep the package generic; no Khudaida-only or downstream-specific public behavior.
- Do not broaden `epcsaft.capabilities()` or documentation claims beyond executable evidence.
- Preserve exact gradient and exact Jacobian requirements.
- Do not use approximate derivative substitutes.
- Do not silently fall back from `hessian_mode="exact"` to limited-memory.
- Keep Graphify generated outputs local/untracked unless explicitly requested later.
- `codex/algorithm-registry` and its worktree/branch are out of scope for this tranche.
- Coordinate native `_core` rebuilds in the main thread only.

## Stop Rule

Stop only when a final audit maps current receipts, focused tests, build/doctor state, cleanup, and git status back to the original tranche request.

If exact Hessian coverage cannot be made all-route complete in this tranche, record the exact unsupported route/backend gap and ensure `hessian_mode="exact"` fails loudly with tests rather than silently selecting limited-memory.

## Canonical Board

Machine truth lives at:

`docs/goals/ipopt-tranche1/state.yaml`

If this charter and `state.yaml` disagree, `state.yaml` wins for task status, active task, receipts, verification freshness, and completion truth.

## Run Command

```text
/goal Follow docs/goals/ipopt-tranche1/goal.md.
```
