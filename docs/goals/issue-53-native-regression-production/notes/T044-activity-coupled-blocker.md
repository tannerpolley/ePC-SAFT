T044 receipt: activity-coupled derivative blocker
=================================================

Result
------

Stopped the first activity-coupled worker tranche at its `stop_if` boundary.

The current native code does **not** yet expose the backend-unavailable-free
state-level derivative substrate needed to differentiate
`mole_fraction_activity` reaction residuals honestly.

What was checked
----------------

1. Activity-coupled reaction residuals in
   `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp`
   use:

   - `log(x_i * gamma_i)`

2. Neutral activity coefficients are built from fugacity coefficients:

   - `gamma_i = exp(ln_phi_i - ln_phi_ref_i)`

3. Fugacity coefficients are built from residual chemical potentials in
   `src/epcsaft/native/epcsaft_fugcoef.cpp`.

4. The current value-level fugacity path depends on
   `ResidualChemicalPotentialResult` /
   `CompositionContributionResult`, which only expose:

   - residual Helmholtz contribution values
   - first composition derivatives `dadx`
   - `sum_x_dadx`
   - compressibility-related value terms

Why this blocks T044
--------------------

To differentiate activity-coupled reaction residuals with respect to the
chemical-equilibrium variables without Backend unavailables, the package needs a
backend-unavailable-free derivative of:

- `ln_gamma_i`
- equivalently `ln_phi_i`

with respect to the active composition variables under density closure.

The current native derivative substrate does **not** expose that.

For the neutral-fugacity activity path, differentiating `ln_phi_i` requires a
state-level derivative surface beyond the currently exposed value + first
composition-derivative bundle.

For the ionic component-activity path, the requirement is at least as broad,
because the reference-state activity construction is more complex than the
neutral fugacity-reference case.

Concrete missing substrate
--------------------------

What is missing is not just a small glue change in
`epcsaft_chemical_equilibrium.cpp`. The next architecture slice needs one of:

- a native backend-unavailable-free `d ln(phi_i) / d log(n_j)` or equivalent
  state-level fugacity-derivative result for supported states, or
- a broader second-derivative/Hessian-capable EOS composition-derivative
  substrate that can be used to build that surface honestly

Why bubble and non-logK still stay downstream
---------------------------------------------

- Bubble rows remain explicitly `not_differentiated`.
- Non-`logK` thermodynamic parameters (`d_born`, `f_solv`, `k_ij`) still need
  the activity-coupled derivative path before they become meaningful in the
  production thermo-regression hot loop.

Evidence
--------

- `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp`
- `src/epcsaft/native/epcsaft_fugcoef.cpp`
- `src/epcsaft/native/epcsaft_mu.cpp`
- `src/epcsaft/native/epcsaft_electrolyte.h`
- `src/epcsaft/native/regression/thermo_regression.cpp`

