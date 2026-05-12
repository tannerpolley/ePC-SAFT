T045 receipt: next architecture after activity-coupled blocker
==============================================================

Decision
--------

Continue issue #53, but narrow the next implementation tranche to the most
relevant supported activity path first:

- build a backend-unavailable-free **state-level activity / fugacity derivative
  substrate for the ionic `epcsaft_component_activity` path**
- then use that to unlock activity-coupled `reactive_speciation` `logK` rows
  in native thermo regression

Do **not** jump to bubble-pressure differentiation or non-`logK` parameter
sensitivities before that substrate exists.

Why this is the next honest step
--------------------------------

1. The current worker blocker is real:

   - activity-coupled residuals already evaluate, but the native code does not
     yet expose a backend-unavailable-free state-level derivative of
     `ln(gamma)` / `ln(phi)` under density closure

2. The ionic component-activity path is the best next slice within issue #53:

   - it is directly relevant to MEA-like electrolyte regression
   - it is a smaller activity-reference problem than the neutral
     component-rich fugacity-reference path
   - it is the prerequisite for meaningful `d_born`, `f_solv`, and `k_ij`
     thermodynamic sensitivities in the native thermo-regression hot loop

3. Bubble rows stay downstream:

   - they are still explicitly `not_differentiated`
   - they need a differentiated outer bubble residual, not just the missing
     activity-coupled speciation Jacobian

Architecture target
-------------------

The next worker should not try to solve all activity-coupled cases at once.
It should target:

- supported nonassociating ionic states only
- `epcsaft_component_activity` only
- `reactive_speciation` `logK` rows only

If that slice still requires a broad second-derivative EOS rewrite, the next
blocker should be recorded explicitly at the fugacity/activity derivative
substrate level.

What remains out of scope for that tranche
------------------------------------------

- neutral `epcsaft_neutral_fugacity_activity` derivative support
- reactive-electrolyte bubble differentiation
- native non-`logK` parameter sensitivities

Evidence
--------

- `docs/goals/issue-53-native-regression-production/notes/T044-activity-coupled-blocker.md`
- `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp`
- `src/epcsaft/native/epcsaft_activity.cpp`
- `src/epcsaft/native/regression/thermo_regression.cpp`
- `tests/native/test_chemical_equilibrium_native.py`

