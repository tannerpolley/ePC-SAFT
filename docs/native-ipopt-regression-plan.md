# Native IPOPT Regression Toolkit for ePC-SAFT

## Summary

This document is a handoff plan for replacing the current Python/SciPy regression workflow with a native C++ regression engine built on IPOPT, while keeping a thin Python API wrapper for user-facing calls and results.

The first implementation scope is intentionally narrow:

- pure-neutral fitting only
- fitted parameters limited to \(m\), \(s\), and \(e\)
- weighted least-squares objective over liquid-density and pure-VLE fugacity-balance data
- exact first derivatives in production
- staged formulation where density remains a nested native solve in v1, with explicit density NLP variables deferred to a later phase

The goal is to make ePC-SAFT own its regression toolkit directly rather than trying to force the runtime into `parmest`. Useful `parmest` ideas should still be borrowed where they improve the workflow, especially typed experiment organization, structured results, and future extensibility for uncertainty analysis.

## Architecture

### Native Regression Subsystem

Add a dedicated native regression subsystem with the following responsibilities:

- typed regression record and residual-family data structures
- normalized regression problem/config structs
- residual evaluator that calls the native EOS runtime
- IPOPT `Ipopt::TNLP` adapter
- native solver/result payload

This subsystem should live separately from the core EOS property runtime conceptually, even if it reuses the existing native mixture and state code heavily.

### Optional IPOPT-Linked Extension

IPOPT support should be compiled as an optional feature. The core EOS extension must remain buildable and importable without IPOPT installed.

The recommended structure is:

- keep the current core `epcsaft.epcsaft` extension independent of IPOPT
- add a separate optional regression-linked extension or native bridge that links against IPOPT
- expose a runtime capability check so regression calls fail clearly when the package was built without IPOPT support

This keeps ordinary installs lightweight and avoids making the baseline thermodynamic runtime depend on solver toolchain availability.

### Thin Python Wrapper

Keep the current public Python regression surface, but make it a thin wrapper over the native engine:

- `fit_pure_neutral(...)`
- `FitProblem`
- `FitResult`
- `load_regression_records(...)`
- `write_fit_result(...)`

The Python layer should remain responsible for:

- user-facing argument normalization
- CSV/tabular record loading
- result-object rendering
- dataset write-back helpers

The Python layer should not remain the owner of the optimization loop.

## Optimization Formulation

### v1 Decision Variables

The v1 IPOPT model should optimize only the fitted pure-component parameters:

\[
\theta = [m, s, e]
\]

with simple bound constraints.

### Objective

Use weighted least squares over typed residual families:

\[
\min_{\theta} \; \frac{1}{2} \sum_i w_i r_i(\theta)^2
\]

where the v1 residual families are:

- density residuals
- pure-VLE fugacity-balance residuals

The weighting and residual-family bookkeeping should be explicit in the native problem representation, not implied by ad hoc array concatenation.

### Density Handling

Use the staged formulation selected for this effort:

- v1 keeps density as a nested native solve using the existing root-selection and stability logic
- a later v2 may promote density to an explicit NLP variable with closure constraints

This means IPOPT sees only parameter variables in v1, while the density closure remains internal to the EOS runtime.

### Failure Handling

Do not preserve the current large-penalty style used in the Python regression workflow. Invalid thermodynamic evaluations should instead be treated as solver-evaluation failures, with the implementation relying on:

- conservative parameter bounds
- reasonable initial guesses
- record pre-screening
- deterministic multistart logic outside IPOPT

This avoids hiding invalid regions behind arbitrary penalty residuals.

## Derivatives

### First Derivatives

Production v1 should use exact first derivatives. Do not rely on finite-difference gradients as the main implementation path.

The current native code already contains autodiff infrastructure, but it is mainly oriented around existing composition/dielectric derivative paths rather than direct regression sensitivities. The regression implementation should extend that infrastructure to support parameter sensitivities for:

- \(m\)
- \(s\)
- \(e\)

### Density Sensitivities

For density-based residuals, compute \(\mathrm{d}\rho / \mathrm{d}\theta\) through implicit differentiation of the pressure-density closure after the nested density solve converges.

This keeps the v1 formulation reduced-space while still giving exact first derivatives to IPOPT.

### Pure-VLE Sensitivities

For the pure-VLE fugacity-balance residuals, compute exact derivatives by combining:

- native autodiff for direct property sensitivities
- chain-rule terms through the liquid and vapor density closures

### Second Derivatives

Use IPOPT with limited-memory Hessian approximation in v1:

- `hessian_approximation=limited-memory`

Exact Hessians are explicitly deferred.

## Parmest-Inspired Ideas To Keep

Do not use `parmest` directly, but retain the useful design ideas that fit a native ePC-SAFT workflow:

- typed experiment/record organization
- explicit residual-family grouping
- structured fit results with family metrics and solver statistics
- a clean future seam for covariance, bootstrap, and leave-out workflows

These should be treated as architecture cues, not as a requirement to model the runtime as a Pyomo experiment system.

## Build and Packaging

IPOPT support is optional and should be documented that way.

### Supported Build Path

The primary supported build path should be a Conda environment that already provides IPOPT headers and libraries. This is the most practical route for the current repo and developer environment.

### Non-IPOPT Builds

The package must still:

- build cleanly without IPOPT
- import cleanly without IPOPT
- expose a clear actionable error when a regression call requires native IPOPT support but the installed build does not include it

### Build-System Implications

Any native-source split required by the new regression layer must also update the repo’s editable-build tracking, especially:

- `setup.py`
- `scripts/build_epcsaft.py`

so the tracked translation units and rebuild-stamp logic remain correct.

## Test Plan

### Build and Import

- build without IPOPT installed
- build with IPOPT installed
- verify non-regression imports work in both cases
- verify regression entrypoints raise a clear error in non-IPOPT builds

### Derivative Verification

- compare native exact gradients against finite differences on representative pure-neutral records
- add a debug-only small-case path that enables IPOPT’s derivative checker

### End-to-End Regression

- run the methane/ethane/propane pure-neutral benchmark through the native IPOPT path
- require parity with the current expected fit tolerances and residual metrics
- verify deterministic multistart selection

### API Compatibility

- keep `FitProblem` and `FitResult` close enough to the current public contract that existing downstream usage remains straightforward
- keep `write_fit_result(...)` working as a Python-side dataset persistence helper
- keep `fit_pure_ion(...)` and `fit_binary_pair(...)` as explicit deferred surfaces

## Assumptions

- No standalone CLI is included in v1.
- `write_fit_result(...)` remains Python-side even after optimization moves native.
- SciPy can remain in the project temporarily during migration and be removed only after the native regression cutover is complete.
- Electrolyte fitting, binary fitting, covariance estimation, bootstrap workflows, leave-out workflows, and explicit density-constraint NLP formulations are deferred.
- The handoff target is another implementation-capable thread/worktree, so this document is written to be self-contained without depending on chat context.

## Acceptance Criteria

The first implementation phase is complete when all of the following are true:

- ePC-SAFT has a native IPOPT-backed regression engine for pure-neutral fitting
- the public Python API remains a thin wrapper rather than the owner of the optimization loop
- v1 uses exact first derivatives
- v1 keeps density as a nested native solve
- IPOPT remains an optional build dependency
- the methane/ethane/propane benchmark passes through the new path with expected parity
