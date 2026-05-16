# T018 Scout: Phase-State Sensitivity Map

Date: 2026-05-15

## Required Jacobian Chain

Issue #116 production Ceres needs the Jacobian of the electrolyte LLE residual vector with respect to transformed variables:

1. transformed variables to beta/formula composition,
2. formula composition to explicit public species composition,
3. explicit composition to density closure,
4. density and composition to fugacity coefficients,
5. fugacity coefficients to neutral and salt-pair phase-equilibrium residuals,
6. formula/material transform to material residual rows if retained.

## Current Owners

- Transformed electrolyte variables and residual evaluation live in `src/epcsaft/native/epcsaft_equilibrium.cpp`.
- Phase-state evaluation is `phase_state(...)` in `epcsaft_equilibrium.cpp`; it calls `mixture->solve_density_scoped(...)`, then `mixture->state(...)`, then `ln_fugacity_coefficient()`.
- Density root scanning and solve diagnostics live in `src/epcsaft/native/epcsaft_density.cpp` via `density_solve_report_cpp(...)`.
- Fugacity coefficient assembly lives in `src/epcsaft/native/epcsaft_fugcoef.cpp` via `fugacity_coefficient_result_cpp(...)`.
- Residual chemical potentials and composition derivative terms live in `src/epcsaft/native/epcsaft_mu.cpp` and `src/epcsaft/native/epcsaft_ares.cpp`.
- Existing Ceres ownership patterns live in `src/epcsaft/native/epcsaft_regression.cpp`, where cost functions call explicit residual/Jacobian evaluators and report `cppad_implicit`.

## Existing Derivative Surfaces

- `composition_derivative_residual_helmholtz_result_cpp(...)` returns contribution-level composition derivative terms and backend labels.
- `association_density_response_cppad_cpp(...)` records an association-density implicit derivative helper for association internals.
- `cppad_eos_contribution_derivatives_cpp(...)` exposes CppAD contribution derivatives for selected EOS pieces.
- `cppad_pressure_density_derivative_cpp(...)` exists as a pressure-density derivative API, but it is not a complete composition-coupled density closure sensitivity for phase-state fugacity residuals.
- Regression Ceres paths already demonstrate the desired pattern: a dedicated residual/Jacobian evaluator feeds a `ceres::CostFunction`.

## Missing Exported Surface

No current native API returns:

- `d rho / d x` for the selected liquid density root,
- `d ln(phi_i) / d x_j` at the solved density with implicit density dependence,
- chain-rule derivatives from transformed electrolyte variables to electrolyte residual rows,
- a row-major Jacobian for `_evaluate_electrolyte_lle_residual_native`.

Without that surface, T017 cannot truthfully set `jacobian_available = true` or use Ceres as the accepted production electrolyte LLE solver.

## Recommended Next Boundary

Do not resume T017 as a direct Ceres implementation. Create a Judge decision task that either:

- creates a dedicated derivative-foundation child goal for phase-state fugacity sensitivities, or
- expands the Worker boundary to include a full phase-state sensitivity implementation across `epcsaft_density.cpp`, `epcsaft_fugcoef.cpp`, `epcsaft_mu.cpp`, `epcsaft_ares.cpp`, `epcsaft_core_internal.h`, and equilibrium tests.

The smallest safe implementation path is likely a new native phase-state sensitivity result type with tests that compare analytic or CppAD-implicit rows against exact known invariants, not against a derivative-approximation oracle.
