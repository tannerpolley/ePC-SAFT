# Native Phase-State Sensitivities For Electrolyte LLE

Parent goal: `docs/goals/native-lle-reactive-production-solvers/goal.md`

This child goal unblocks issue #116 by implementing the native phase-state sensitivity foundation needed for the electrolyte LLE transformed-variable residual Jacobian and Ceres solve.

Completion requires a verified native surface that can support the electrolyte LLE residual Jacobian without placeholder, stale, identity-only, or derivative-approximation derivatives.
