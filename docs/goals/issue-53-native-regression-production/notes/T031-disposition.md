T031 disposition
================

Decision
--------

Use the current branch as a partial PR candidate only if the user explicitly
wants a partial issue #53 tranche published. Do not close issue #53.

Why
---

The branch now contains independently useful production progress:

- native thermodynamic regression is reachable from the public reactive
  regression wrapper for the supported speciation/logK slice;
- bundled Ceres builds on Windows/MinGW without LAPACK link failures;
- a Ceres-enabled local benchmark proves `optimizer_backend=ceres`,
  `derivative=analytic_implicit`, `native_hot_loop=True`, and objective
  reduction to numerical zero;
- docs/tests explicitly prevent overclaiming the remaining Born and bubble
  derivative gaps.

But issue #53 is still not complete because mixed pressure/speciation reactive
electrolyte production fitting needs additional native derivative architecture.

Required partial PR wording
---------------------------

Recommended title:

`Issue #53 partial slice: native thermo regression routing and Ceres build proof`

The PR body must say:

- this does not close issue #53;
- supported production slice is reactive speciation rows with
  ideal-mole-fraction reaction standard states, linear speciation targets, and
  reaction logK parameters;
- Born-SSM+DS `d_born`/`f_solv` and reactive-electrolyte bubble-pressure
  derivatives remain `backend_unavailable`;
- Ceres/CppAD proof was run locally with the validation commands in T029.

Next implementation plan
------------------------

A new broad native-derivative tranche should start with design, not a direct
patch:

1. Choose bubble residual/Jacobian or Born-SSM+DS scalar-templating as the next
   single owning path.
2. Define the residual variables and sensitivity contract.
3. Add native tests that fail unless derivatives are analytic/implicit/CppAD,
   never production finite differences.
4. Only then widen `thermo_derivative_supported(...)`.
