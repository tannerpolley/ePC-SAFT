T037 receipt: concentration-standard-state derivative blocker
===========================================================

Result
------

Stop on implementation for this tranche as currently scoped.

Blocked objective
-----------------

`T037` aimed to implement native concentration-standard-state
chemical-equilibrium derivatives and then expose that support to native
thermodynamic regression.

Verified blocker
----------------

Concentration-standard-state reaction residuals require density sensitivity:

- for `STANDARD_STATE_CONCENTRATION`, native chemical equilibrium uses
  `species_activity = x_i * molar_density`
  in `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp`.

That means the residual derivatives need, at minimum:

1. `d log(rho) / d log(n)` for the inner chemical-equilibrium Jacobian; and
2. `d log(rho) / d theta` for regression parameter sensitivities.

At constant `T` and `P`, those density sensitivities require non-finite-
difference pressure-composition and pressure-parameter derivatives through the
density closure.

What is missing today
---------------------

There is no existing native pressure-composition derivative surface:

- `CompressibilityFactorResult` only exposes scalar `z` contribution totals;
- `CompositionContributionResult` exposes `dadx`, `ares`, and `z_raw` term
  values, but not `dP/dx` or `dZ/dx`;
- `solve_density_scoped(...)` and the current density validity path use a
  finite-difference `dpdrho` check inside `epcsaft_density.cpp`;
- there is no package-owned native `dpdx`, `dPdx`, or parameterized pressure
  derivative callback that `chemical_equilibrium.cpp` can reuse for production
  derivatives.

Because of that, the concentration-standard-state derivative tranche cannot be
finished honestly without first adding a new native pressure-composition
derivative substrate.

Why this is a real blocker
--------------------------

Using finite differences here would violate the issue #53 production derivative
policy. Pretending concentration derivatives are available without a real
pressure-composition derivative path would be false.

Minimum next prerequisite
-------------------------

Add a native finite-difference-free pressure derivative substrate for the
general EOS state, sufficient to recover:

- `dp/drho` from a production derivative path;
- `dP/dx` at fixed `T`, `rho`, and composition;
- later, `dP/dtheta` for fitted thermodynamic parameters.

Only after that can concentration-standard-state chemical-equilibrium
derivatives be implemented honestly.
