# T017 Blocker: Production Jacobian Requires Phase-State Sensitivities

Date: 2026-05-15

## Environment Gate

`uv run python scripts/dev/build_epcsaft.py --profile full` passed. Runtime build info reports both Ceres and CppAD compiled and available.

## Blocker

T017 cannot honestly implement the production Jacobian and Ceres accepted electrolyte LLE solve yet because the residual Jacobian requires sensitivities of the phase-state fugacity residuals with respect to transformed electrolyte LLE variables.

The current residual evaluator from T016 can return finite residuals, compositions, densities, phase fraction, and diagnostics. It still reports `jacobian_available = false` because the required derivative chain is not implemented.

## Missing Derivative Chain

The required Jacobian path must cover:

- transformed variables to organic formula composition,
- transformed variables plus feed formula balance to aqueous formula composition,
- formula composition to explicit public species composition,
- composition to density closure,
- density and composition to fugacity coefficients,
- neutral and salt-pair equilibrium residuals,
- material-balance residual rows if retained in the residual vector.

Existing native CppAD and analytic-implicit derivative surfaces cover several EOS and regression pieces, but not this complete phase-equilibrium residual Jacobian as an exported native surface.

## Decision

Stop T017 source edits and create a read-only Scout task to map the exact phase-state sensitivity owner files and the smallest implementable derivative path. Do not add a placeholder Jacobian or label the residual-only evaluator as Ceres production.
