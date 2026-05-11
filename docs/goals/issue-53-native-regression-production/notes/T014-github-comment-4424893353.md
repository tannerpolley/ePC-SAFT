# T014 GitHub Comment 4424893353 Intake

Source: https://github.com/tannerpolley/ePC-SAFT/issues/53#issuecomment-4424893353

## Decision

The current branch must not be merged or described as completing issue #53 unless it is clearly framed as a partial Phase 0 / native residual-record contract slice.

## Current Branch Framing

The branch has useful groundwork:

- native result/status/schema contracts
- fixed-shape residual-record evaluation
- native residual-record solve boundary
- finite-difference production gate at that boundary
- benchmark scaffolding

It does not yet satisfy #53 because it does not provide a full native Ceres/CppAD thermodynamic parameter-iteration loop for mixed pressure/speciation rows.

## Required Next Tranche

1. Package-wide CppAD native AD layer, not a regression-only adapter.
2. Finite difference debug-only everywhere, gated by `EPCSAFT_ALLOW_FINITE_DIFFERENCE_DEBUG=1`.
3. Native C++ thermodynamic row evaluator for at least `ReactiveSpeciation` and `ReactiveElectrolyteBubble`.
4. Ceres must execute the production solve loop through a native dynamic cost-function path.
5. Use analytic/CppAD plus implicit sensitivities for nested speciation and bubble solves.
6. Capabilities remain honest until `native_hot_loop=true` is actually true.
7. Tests must prove real thermodynamic parameter movement, objective decrease, non-finite-difference derivatives, and no Python production objective.
8. Benchmarks must report `optimizer_backend=ceres`, `native_hot_loop=true`, `python_objective_used=false`, and `derivative_backend != finite_difference`.

## Misfire To Avoid

Do not complete the board by publishing a residual-record-only branch. That branch can only be a partial PR unless the next tranche is implemented.
