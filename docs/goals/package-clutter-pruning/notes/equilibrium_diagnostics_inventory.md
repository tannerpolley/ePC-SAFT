# Equilibrium Diagnostics Inventory

Phase 3 source scope from issue #120:

- `src/epcsaft/equilibrium_core/confidence.py`
- `src/epcsaft/equilibrium_core/thermo_diagnostics.py`

## Classification

| Original runtime path | Contents | Classification | New owner | Runtime package status |
| --- | --- | --- | --- | --- |
| `src/epcsaft/equilibrium_core/confidence.py` | Khudaida electrolyte LLE confidence suite, benchmark case loading, report/CSV/plot writers, CLI entry point. | Benchmark confidence score and analysis/report helper. | `scripts/validation/equilibrium_core/confidence.py` | Thin compatibility shim only. |
| `src/epcsaft/equilibrium_core/thermo_diagnostics.py` | Fixed-phase Khudaida validation diagnostics, cached tie-line comparisons, digitized-paper comparisons, solver-gate helpers. | Validation/report helper, not production result-object diagnostics. | `scripts/validation/equilibrium_core/thermo_diagnostics.py` | Thin compatibility shim only. |

## Runtime Diagnostics Kept

Production solve-result diagnostics remain on equilibrium results and state diagnostics payloads. This phase did not change solver-result fields such as convergence state, residual norms, material balance, charge balance, phase distance, iterations, messages, activity/fugacity payloads, or public result object behavior.

## Import Updates

- `scripts/validation/validate_electrolyte_lle_confidence.py` imports the relocated confidence implementation directly.
- `tests/workflows/validation/equilibrium_core/*` import relocated validation helpers directly.
- Package plot-smoke analysis tests import relocated validation helpers directly.
- Khudaida paper-validation scripts import relocated validation helpers directly.

## Test Boundary

The relocated confidence and Khudaida thermodynamic diagnostics tests now live under
`tests/workflows/validation/equilibrium_core` because they validate script-owned
reports, fixture matrices, and paper-validation helper behavior rather than core
equilibrium solver contracts.

Default workflow tests keep schema, fixture, charge-neutrality, cached-matrix, and
report-helper checks that are available in this checkout. Hard Khudaida native
validation cases require `EPCSAFT_KHUDAIDA_VALIDATION=1`, full confidence report
generation requires `EPCSAFT_EQUILIBRIUM_CONFIDENCE=1`, and generated digitized
CSV inputs require `EPCSAFT_KHUDAIDA_DIGITIZED_DATA=1`.

## Compatibility Decision

`src/epcsaft/equilibrium_core/confidence.py` and `src/epcsaft/equilibrium_core/thermo_diagnostics.py` remain only as import-forwarding compatibility shims. They contain no validation implementation logic.

## Acceptance Notes

The runtime package no longer owns benchmark confidence scoring, validation reports, or Khudaida analysis helpers. Scientific equations, solver residuals, and result diagnostics were not changed.
