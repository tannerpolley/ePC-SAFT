T030 receipt: derivative tranche stop condition
===============================================

Decision
--------

Stop on implementation for this tranche. The remaining issue #53 derivative
work is real but not safely small enough for the current branch without a broad
native scalar-templating or implicit-derivative design pass.

Bubble-pressure derivative gap
------------------------------

Current state:

- `evaluate_native_thermo_regression_rows(...)` evaluates
  `reactive_electrolyte_bubble` rows in native C++.
- The row diagnostic explicitly reports `derivative_backend =
  "not_differentiated"`.
- `thermo_derivative_supported(...)` rejects bubble rows before Ceres fitting
  and returns `backend_unavailable`.

Required production work:

- define residual variables for log-pressure and vapor composition;
- expose residual assembly for fugacity equalities and vapor normalization;
- differentiate liquid fixed-composition fugacity/activity state with respect
  to fitted parameters;
- differentiate vapor submixture fugacity state with respect to pressure,
  composition, and fitted parameters;
- solve the implicit sensitivity system and map it into target residual
  derivatives.

Born-SSM+DS derivative gap
--------------------------

Current state:

- `apply_native_thermo_parameters(...)` can apply `d_born` and `f_solv` into
  native mixture args.
- `thermo_derivative_supported(...)` currently allows only reaction
  log-equilibrium constants for the native thermodynamic Ceres fit.
- `contribution_dadx_autodiff_cpp(...)` still throws for SSM/DS Born
  composition derivatives (`born_model == 2`).

Required production work:

- implement scalar-templated or CppAD-backed SSM/DS Born contribution
  derivatives;
- carry those parameter sensitivities through residual chemical potential,
  fugacity/activity coefficients, and any pressure-density closure used by the
  row;
- extend `fill_implicit_speciation_jacobian(...)` to include non-logK
  `r_theta` terms for `d_born` and `f_solv`;
- add tests proving no backend-unavailable production derivatives are used.

Reason not to patch further here
--------------------------------

Both paths require more than a local gate change. Enabling parameter kinds in
`thermo_derivative_supported(...)` without a correct `r_theta` chain would
produce false Ceres confidence. Using Backend unavailables would violate the issue
#53 production rule.

Validation already completed before this stop
---------------------------------------------

- Ceres/CppAD local build now succeeds after disabling LAPACK in bundled Ceres.
- Native thermo logK benchmark now proves Ceres execution with
  `native_hot_loop=True` and objective reduction to numerical zero.
- Tests still codify Born and bubble gaps as `backend_unavailable`.

