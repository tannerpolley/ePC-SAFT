# T016 Residual Surface Note

Date: 2026-05-15

## Scope

T016 implements the native electrolyte LLE residual evaluation surface that T005 lacked. It does not claim the production Jacobian or Ceres accepted-solve route.

## Implemented Surface

`_evaluate_electrolyte_lle_residual_native` now returns a payload for transformed electrolyte LLE variables instead of throwing immediately. The payload includes:

- transformed variables and bounds,
- phase-equilibrium residual block,
- material-balance residual block,
- aqueous and organic phase compositions,
- phase fugacity coefficients,
- phase densities,
- organic phase fraction,
- material-balance, charge-balance, phase-distance, Gibbs, and objective diagnostics,
- residual block sizes and basis diagnostics.

## Deliberately Not Claimed

- `jacobian_available` remains false.
- `jacobian_backend` remains `not_available`.
- The accepted production electrolyte LLE route is still not a Ceres trust-region solve.

## Next Required Work

T017 must implement the real transformed-variable Jacobian and Ceres ownership. Issue #116 is still incomplete until T017, T006, T007, and T008 pass.
