# T019 Judge Decision: Create Phase-State Sensitivity Child Goal

Date: 2026-05-15

Decision: create_derivative_foundation_child_goal

## Rationale

Issue #116 cannot complete the production electrolyte LLE Ceres solve until the native phase-state sensitivity chain exists. That work crosses density closure, fugacity coefficients, residual chemical potentials, composition derivative terms, CppAD support, and equilibrium residual assembly. It is too broad to fold into the original T017 Worker without losing reviewability or risking a false production Jacobian.

## Decision

Create a child goal under this board for the derivative foundation:

`docs/goals/native-lle-reactive-production-solvers/subgoals/phase-state-sensitivities/state.yaml`

The child goal must produce a native phase-state sensitivity surface before T017 can resume. The parent issue #116 work remains blocked at the Ceres/Jacobian slice until that child goal completes.

## Required Child Goal Proof

- Native API returns density and fugacity coefficient sensitivities for liquid phase states.
- Sensitivities include implicit density dependence or explicitly documented fixed-density scope with a later implicit task.
- The electrolyte LLE residual evaluator can consume those sensitivities to produce a transformed-variable Jacobian.
- Tests fail if the Jacobian is placeholder, stale, identity-only, or disconnected from the actual residual state.
- No accepted production electrolyte LLE result claims Ceres completion until the child proof is integrated.
