# Native IPOPT Regression Toolkit for ePC-SAFT

## Summary

ePC-SAFT now ships a native C++ IPOPT-backed regression path for the public `fit_pure_neutral(...)` workflow.

The package no longer uses a Python/SciPy optimization loop for the supported regression surface. Python remains responsible for record loading, argument normalization, result packaging, and dataset write-back, while the optimization model and solver callbacks live in the native runtime.

## Implemented v1 Scope

The first supported regression phase is intentionally narrow:

- one nonassociating neutral component only
- fitted parameters limited to \(m\), \(s\), and \(e\)
- weighted least-squares objective over two residual families:
  - liquid-density residuals
  - pure-VLE fugacity-balance residuals
- exact first derivatives
- density retained as a nested native solve
- IPOPT configured with limited-memory Hessian approximation

Ion regression, binary regression, associating pure-neutral regression, covariance estimation, bootstrap workflows, and explicit density-constrained NLP formulations remain deferred.

## Native Ownership

The native regression layer lives in `src/epcsaft/native/epcsaft_regression.cpp` and is linked directly into the main `epcsaft.epcsaft` extension.

Core responsibilities of that layer:

- typed density and pure-VLE record structs
- deterministic multistart generation
- IPOPT `Ipopt::TNLP` implementation
- exact reduced-space objective/gradient evaluation
- native result payloads for fitted values, family metrics, and solver statistics

The public Python wrapper in `src/epcsaft/regression.py` now:

- normalizes flat records
- resolves fixed pure-component metadata
- prepares bounds and initial guesses
- marshals payloads into the Cython/native seam
- repackages the native solver result into `FitResult`

## Exact-Derivative Strategy

The v1 implementation uses exact first derivatives.

Parameter sensitivities for \(m\), \(s\), and \(e\) are evaluated through repeated forward-mode autodiff passes in the native property model. Density-coupled terms use implicit differentiation after the nested density closure converges:

\[
\frac{\mathrm{d}\rho}{\mathrm{d}\theta}
=
-\frac{\partial p / \partial \theta}{\partial p / \partial \rho}
\]

That exact \(\partial p / \partial \rho\) path is evaluated in the native regression/property seam rather than reusing the finite-difference validator logic from the density root checker.

## Failure Handling

The native IPOPT callbacks do not synthesize large penalty residuals. If the EOS or density closure cannot be evaluated safely at a trial point, the native callback reports an evaluation failure to IPOPT.

This keeps invalid thermodynamic regions out of the objective algebra and makes the reduced-space NLP behavior more defensible.

## Build and Packaging

IPOPT is now a required build dependency for the package.

The supported developer path is the active Conda environment, with IPOPT headers and libraries discovered from that environment first. On Windows, the current build resolves IPOPT from locations under `%CONDA_PREFIX%\\Library\\include` and `%CONDA_PREFIX%\\Library\\lib`.

The package build now fails early with an actionable message when IPOPT headers or libraries are missing.

## Public API

The public regression surface remains:

- `fit_pure_neutral(...)`
- `FitBounds`
- `FitProblem`
- `FitResult`
- `load_regression_records(...)`
- `write_fit_result(...)`

`FitResult` keeps the existing compatibility fields such as `cost`, `status`, `message`, and `nfev`, and now also reports `backend="ipopt_native"`.

## Validation

The current cutover is validated by:

- editable build/import of the main extension with IPOPT linked
- hydrocarbon methane/ethane/propane benchmark parity checks
- regression API contract tests
- native exact-gradient validation against finite differences on a representative pure-neutral case

## Deferred Work

The main deferred extensions after v1 are:

- associating pure-neutral regression
- ion and binary regression
- exact Hessians
- covariance / bootstrap / uncertainty workflows
- explicit density-variable NLP formulations
- broader regression record families beyond liquid density and pure-VLE fugacity balance
