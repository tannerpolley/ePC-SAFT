# T005 Blocker: Electrolyte LLE Residual/Jacobian Surface

Date: 2026-05-15

## Blocked Task

T005: issue #116 Stages 4-6, production residual blocks, production Jacobian, and Ceres trust-region residual solve.

## Blocker

The accepted electrolyte LLE production route cannot honestly move to a Ceres trust-region residual solve until the native code has a real residual Jacobian with respect to the transformed electrolyte LLE variables. The current implementation can evaluate the residual at a double-valued candidate, but the residual sensitivity route still reports unavailable.

## Current Evidence

- Current required native build has CppAD and Ceres compiled.
- Ceres support exists for regression problems, but not for electrolyte LLE phase-equilibrium residuals.
- Existing CppAD surfaces cover contribution and regression derivative paths, but there is no complete native residual/Jacobian surface for electrolyte LLE transformed variables, density closure, and fugacity residuals.
- `newton_step` and `_evaluate_electrolyte_lle_residual_native` still report unavailable electrolyte LLE residual sensitivities.
- Implementing T005 by filling a placeholder Jacobian or by relying on a manual numeric perturbation Jacobian would violate issue #116 and the roadmap.

## Decision

Stop T005 source edits until a Judge/PM boundary decision expands or redesigns the implementation package around the missing residual/Jacobian surface. The next active task is a boundary-resolution task, not a source patch.

## Required Next Decision

Decide whether issue #116 T005 needs:

- a new native residual/Jacobian helper under `src/epcsaft/native/equilibrium/**`,
- expanded allowed files in fugacity, density, or composition-derivative native surfaces,
- a child subgoal dedicated to electrolyte LLE CppAD/implicit sensitivities, or
- a revised issue #116 plan if Ceres production LLE cannot be achieved inside this branch.
