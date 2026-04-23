# Archived IPOPT Regression Note

This note is retained only as historical context.

## Current state

ePC-SAFT no longer ships an IPOPT-based regression backend.

The supported regression workflow is now:

- native C++ least-squares ownership inside the main `epcsaft.epcsaft` extension
- public Python surface kept at `fit_pure_neutral(...)`, `FitBounds`, `FitResult`, `load_regression_records(...)`, and `write_fit_result(...)`
- no IPOPT build dependency for the package

## Why IPOPT was removed

Two native IPOPT formulations were implemented and benchmarked:

1. reduced-space callback IPOPT
2. explicit-state IPOPT with pressure-closure constraints and square initialization

Neither formulation beat the native least-squares workflow on the current pure-neutral benchmark surface.

Observed failure modes included:

- excessive iteration count and solve time relative to least-squares
- strong basin sensitivity
- explicit-state drift away from good square-initialized feasible points
- a local Pyomo/`parmest` comparison that also did not outperform least-squares on the same methane case

Because the package already had a fast native least-squares workflow with better practical behavior, the IPOPT implementation was removed rather than kept as dead or misleading complexity.

## What would justify reintroduction

Any future IPOPT reintroduction should clear a higher bar than simple parity.

It should demonstrate at least one of:

- materially better fit quality on the benchmark suite
- materially better runtime on the benchmark suite
- support for a constrained regression surface that least-squares cannot handle cleanly

Until that bar is met, least-squares remains the authoritative regression engine for the package.
