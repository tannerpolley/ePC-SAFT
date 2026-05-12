T043 receipt: activity-coupled derivative prerequisite
=====================================================

Result
------

After the supported ideal and concentration `cppad_implicit` thermo slices,
the next highest-value remaining production derivative gap is
**activity-coupled standard states**.

Decision
--------

Choose the next implementation target as:

- finite-difference-free native activity-coupled chemical-equilibrium Jacobians
  for the supported nonassociating state slice

Do **not** choose bubble-pressure differentiation or non-`logK` thermodynamic
parameter sensitivities first.

Why this is the next honest target
----------------------------------

1. `thermo_regression.cpp` still rejects activity-coupled standard states
   directly:

   - `implicit reactive-speciation derivatives do not yet support
     activity-coupled reaction standard states.`

2. `thermo_regression.cpp` still rejects non-`logK` thermodynamic parameters
   directly:

   - `native Ceres thermodynamic derivatives currently support reaction
     log-equilibrium constants only.`

3. Bubble rows are still explicitly undifferentiated:

   - row diagnostics say `derivative_backend = "not_differentiated"`
   - the native thermo fit still rejects bubble rows because the current
     derivative loop supports `reactive_speciation` rows only

4. Non-`logK` parameters such as `d_born`, `f_solv`, and `k_ij` only become
   meaningful in the production thermodynamic regression path once the
   reaction residual can differentiate through activity-coupled standard
   states. The supported ideal and concentration `logK` slices do not provide
   that path.

Implication
-----------

The real next prerequisite is a state-level activity derivative substrate for
supported nonassociating states, so that activity-coupled reaction residuals
can stop failing closed as `backend_unavailable`.

That worker tranche should aim to produce:

- native activity-coupled chemical-equilibrium Jacobians for supported
  nonassociating states
- native thermo regression support for activity-coupled `logK` rows

What still remains after that
-----------------------------

Even after an activity-coupled `logK` slice lands, the following would still
remain separate work:

- reactive-electrolyte bubble-pressure differentiation
- non-`logK` thermodynamic parameter sensitivities (`d_born`, `f_solv`,
  `k_ij`, etc.)

Evidence
--------

- `src/epcsaft/native/regression/thermo_regression.cpp`
- `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp`
- `tests/native/test_native_ceres_reactive_pressure_speciation.py`
- `tests/native/test_native_ceres_thermodynamic_regression.py`
