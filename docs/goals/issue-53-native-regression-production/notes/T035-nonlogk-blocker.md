T035 receipt: non-logK reactive-speciation blocker
===============================================

Result
------

Stop on implementation for this task as scoped.

Blocked claim
-------------

The active T035 objective assumed that the first real non-`logK`
reactive-speciation parameter-sensitivity path could be `d_born`, `f_solv`, or
generic binary parameters such as `k_ij` inside the currently supported native
thermodynamic regression slice.

That assumption is false for the current derivative-capable slice.

Verified reason
---------------

In native chemical equilibrium, the reaction residual for
`standard_state == STANDARD_STATE_IDEAL_MOLE_FRACTION` is:

- `-log_k + sum(nu_i * log(x_i))`

and does not use density or activity coefficients.

File evidence:

- `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp`
  - `evaluate_chemical(...)`
  - for `STANDARD_STATE_IDEAL_MOLE_FRACTION`, `species_activity = x_i`

The current thermo-regression derivative gate only accepts that ideal
mole-fraction standard-state slice:

- `src/epcsaft/native/regression/thermo_regression.cpp`
  - `thermo_derivative_supported(...)`
  - rejects any reaction standard state other than `1`

Therefore:

- `d_born`, `f_solv`, and `k_ij` do not enter the supported reactive-speciation
  residual equations for this slice;
- speciation-target sensitivities for those parameters are identically zero in
  the current supported slice;
- any attempt to “enable” those parameters here would be fake progress unless
  the native chemical-equilibrium derivative backend first supports
  concentration- or activity-coupled standard states.

What this means for issue #53
-----------------------------

The next honest prerequisite is not more glue in thermo regression. It is a new
native chemical-equilibrium derivative tranche for nonideal standard states:

1. concentration-coupled standard states, where density enters the residual;
2. mole-fraction-activity standard states, where activity coefficients enter
   the residual;
3. only after that can native thermodynamic regression honestly expose
   `d_born`, `f_solv`, or related parameter sensitivities through
   reactive-speciation rows.

Validation evidence
-------------------

A focused native test now codifies the invariance of the current ideal
mole-fraction slice with respect to `d_born` and `f_solv`.
