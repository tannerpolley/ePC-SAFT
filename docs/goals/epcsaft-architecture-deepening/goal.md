# ePC-SAFT Architecture Deepening

## Objective

Implement the five architecture deepening changes identified in the architecture review so the `epcsaft` package has deeper, more testable modules around **Equilibrium Problem**, **Production Solver Path**, **Target Dataset**, **Parameter Family**, and **Capability Contract** behavior.

## Original Request

"Lets prep a goal to implement all of these changes."

## Intake Summary

- Input shape: `existing_plan`
- Audience: package maintainers and downstream consumers that depend on generic ePC-SAFT workflows.
- Authority: `requested`
- Proof type: `test`
- Completion proof: all five architecture slices are implemented, tests prove the public generic contracts, `uv run python scripts/dev/check_text_gates.py` and an appropriate validation ladder pass, cleanup passes, and tracked changes are committed locally.
- Likely misfire: GoalBuddy could produce documentation-only architecture notes or shallow wrappers while leaving the old dispatch/result/parameter/capability complexity spread across callers.
- Blind spots considered: native `_core` rebuild coordination, Ipopt/Ceres optional build modes, public compatibility for downstream repos, banned text gates, no overclaiming in capabilities, and preserving the staged-vs-production reactive distinction.
- Existing plan facts: the architecture review named five candidates and recommended starting with **Equilibrium Problem** because it is the widest leverage point.

## Goal Kind

`existing_plan`

## Current Tranche

Complete successive safe, verified architecture slices until the full five-candidate outcome is done. Start by validating the plan against the live branch, `CONTEXT.md`, `docs/roadmaps/FULL_ROADMAP.md`, current tests, and current module ownership. Then implement the largest safe useful slice at a time, beginning with **Equilibrium Problem** unless Judge finds a stronger dependency order.

## Non-Negotiable Constraints

- Keep the package generic; do not add downstream-specific public interfaces or metrics.
- Preserve the package completion standard in `docs/roadmaps/FULL_ROADMAP.md`.
- Keep staged reactive workflows explicitly separate from coupled **Reactive LLE Problem** and **Reactive Phase Equilibrium** production claims.
- Do not broaden **Capability Contract** claims without executable tests and evidence.
- Do not introduce banned exact backend or derivative tokens.
- Coordinate native `_core` rebuilds through the PM thread; avoid clean rebuilds while tests or agents may import `_core`.
- Use `uv run python run_pytest.py <focused-targets> -q` for focused tests and `uv run python scripts/dev/validate_project.py quick` for normal high-level validation when the slice warrants it.
- Finish tracked repo changes with local commits unless the operator explicitly says not to.

## Stop Rule

Stop only when a final audit proves the full original outcome is complete.

Do not stop after planning, discovery, or Judge selection if a safe Worker task can be activated.

Do not stop after a single verified Worker package while the broader five-candidate outcome still has safe local follow-up work.

Do not create one Worker/Judge pair per repeated file, route, or helper. Put repeated same-shape work into one Worker package and review the package as a whole.

## Slice Sizing

Safe means bounded, explicit, verified, and reversible. It does not mean tiny.

A good task is the largest safe useful slice. For this goal, a good Worker task should make one architecture concept materially deeper and prove it through public package contracts or focused native/backend tests.

## Canonical Board

Machine truth lives at:

`docs/goals/epcsaft-architecture-deepening/state.yaml`

If this charter and `state.yaml` disagree, `state.yaml` wins for task status, active task, receipts, verification freshness, and completion truth.

## Run Command

```text
/goal Follow docs/goals/epcsaft-architecture-deepening/goal.md.
```

## PM Loop

On every `/goal` continuation:

1. Read this charter.
2. Read `state.yaml`.
3. Run the bundled GoalBuddy update checker when available and mention a newer version without blocking.
4. Re-check the intake, likely misfire, constraints, live branch, and current dirty state.
5. Work only on the active board task.
6. Assign Scout, Judge, Worker, or PM according to the task.
7. Write a compact task receipt.
8. Update the board.
9. If safe local work remains, choose the next largest reversible Worker package and continue unless blocked.
10. Review at phase, risk, rejected-verification, ambiguity, or final-completion boundaries.
11. Finish only with a Judge/PM audit receipt that maps receipts and verification back to the original five architecture candidates and records `full_outcome_complete: true`.
