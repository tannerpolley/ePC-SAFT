# T015 Judge Decision: Split T005 Residual/Jacobian Work

Date: 2026-05-15

Decision: split_t005_into_residual_surface_then_jacobian_ceres

## Rationale

The original T005 bundled three separate risks: residual-surface ownership, Jacobian availability, and Ceres solve ownership. Current code still throws from `_evaluate_electrolyte_lle_residual_native`, so the next useful source slice is to make that evaluator return a real native residual payload on transformed electrolyte LLE variables.

This does not complete issue #116 Stages 4-6. It creates the surface that the Jacobian and Ceres tasks must use.

## Approved Next Worker

T016 may implement the native electrolyte LLE residual evaluation surface only. It must:

- reuse the existing transformed electrolyte variables from T004,
- return finite residuals, phase compositions, phase fraction, densities, material balance, charge balance, phase distance, and Gibbs diagnostics,
- report residual block sizes and basis diagnostics,
- keep `jacobian_available = false` until a real Jacobian exists,
- avoid Ceres production-route claims.

## Deferred Worker

T017 must handle the production Jacobian and Ceres solve. It must not proceed from placeholder, identity, stale, or manual numeric perturbation Jacobians. It may need expanded native files beyond `epcsaft_equilibrium.cpp` if the density/fugacity sensitivity chain needs new ownership.

## Stop Conditions

- Stop if the residual evaluator would need to call the old accepted predictive solve loop.
- Stop if a test or diagnostic would imply Ceres production before Ceres owns the accepted solve.
- Stop if Jacobian diagnostics cannot be truthfully kept separate from residual-evaluator availability.
