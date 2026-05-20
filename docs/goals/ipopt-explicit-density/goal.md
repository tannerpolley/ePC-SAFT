# Ipopt Explicit Density Goal

## Objective

Make density explicit throughout native Ipopt equilibrium paths where density currently should be an NLP variable instead of an internal pressure-root solve. Focus first on remaining pressure-root, stability, and reactive paths. Normal state/property calls may keep pressure-root density solving because those are not Ipopt NLP evaluations.

## Requirements

- Keep the existing phase-volume formulation where density is already explicit through `rho = sum(n_i) / V`.
- Find every native Ipopt equilibrium path that still calls pressure-root density solves during objective, constraint, Jacobian, or Hessian evaluations.
- Refactor eligible paths to include density or phase volume as an NLP variable with pressure consistency as an explicit nonlinear constraint.
- Preserve exact Hessian support and the shared Lagrangian Hessian assembler.
- Do not lift association site fractions `X_A` into the NLP in this tranche.
- Association implicit handling may be improved only after the explicit-density work is complete enough that remaining failures are clearly association-owned.
- Preserve normal public property/state behavior: `state(T, P, x)` can still solve density outside Ipopt; `state(T, rho, x)` remains direct.

## Completion Proof

The goal is complete when all native Ipopt equilibrium routes where density can and should be explicit have no pressure-root density solve inside Ipopt evaluation callbacks, exact/default Hessian tests pass for the affected route families, electrolyte LLE remains accepted, focused validation passes, cleanup passes, and the work is committed locally.

## Non-Goals

- Do not convert normal property evaluations to explicit-density APIs.
- Do not add route-specific fallback flags.
- Do not hide unsupported derivative paths behind limited-memory fallbacks.
- Do not broaden package capability claims without route evidence.
