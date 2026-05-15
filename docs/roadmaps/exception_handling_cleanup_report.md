# Exception Handling Cleanup Report

This report is the Phase 7 deliverable for issue #120. It classifies broad
exception handling in `src/epcsaft`, narrows clear over-catches, and records why
the remaining broad catches are intentional result-contract boundaries.

## Scope

- Source scan: `rg "except Exception|except BaseException|except:" src/epcsaft`
- Initial source count from Phase 0 intake: 23 broad catches.
- Current source count after cleanup: 2 broad catches.
- No EOS equations, solver math, native residuals, parameter values, or public
  result schemas were changed.

## Changes Made

| Area | Original sites | New boundary | Rationale |
| --- | ---: | --- | --- |
| State construction and native derivative probes in `src/epcsaft/epcsaft.py` | 8 | `_NATIVE_CALL_ERRORS` and `_DERIVATIVE_VALUE_ERRORS` | Keeps native and derivative capability failures mapped into existing public result contracts without swallowing arbitrary Python bugs. |
| Neutral bubble/dew and phase-state evaluation in `src/epcsaft/equilibrium.py` | 3 | `_SOLVER_EVALUATION_ERRORS` or explicit state-construction errors | Preserves candidate-failure recording while preventing unrelated programming errors from being treated as failed thermodynamic candidates. |
| Activity-fixed-point and scalar-search speciation helpers in `src/epcsaft/reactive_speciation.py` | 3 narrowed, 1 retained | `_SPECIATION_EVALUATION_ERRORS` and `_COMPOSITION_NORMALIZATION_ERRORS` | Candidate and normalization failures remain structured, but unexpected bugs are no longer silently skipped in scalar search or initial failure shaping. |
| Runtime metadata import probes in `src/epcsaft/runtime.py` | 3 | `(ImportError, OSError)` | Runtime metadata still handles missing or unloadable native extensions, while unrelated metadata bugs surface loudly. |
| CLI import probe in `src/epcsaft/__main__.py` | 1 | `(ImportError, OSError)` | The CLI still reports import/load failures cleanly, but does not hide arbitrary execution errors. |
| Optional IPOPT adapter in `src/epcsaft/_optional_backends/ipopt.py` | 3 | import/load, status-code, and native-seed typed boundaries | Optional dependency import diagnostics, IPOPT status parsing, and native seed failure recording now have specific exception boundaries. |

## Remaining Broad Catches

| Site | Classification | Why it remains broad |
| --- | --- | --- |
| `src/epcsaft/reactive_regression.py` row evaluation | Intentional batch-result boundary | Each row in a regression objective must become either a successful row result or a failed row result with penalty residuals. The catch intentionally excludes `BaseException`, but captures row-level user data, mixture factory, speciation, bubble, and result-conversion failures so one bad row does not abort objective accounting. |
| `src/epcsaft/reactive_speciation.py` sweep point evaluation | Intentional `error_mode="result"` boundary | In result mode, each sweep point must return a structured failure object instead of aborting the whole sweep. Strict mode still re-raises before conversion, so solver failures are not hidden when callers request strict behavior. |

## Sites Narrowed Or Removed From The Broad Scan

- `src/epcsaft/epcsaft.py`: state construction, native derivative support probes,
  activity derivative value assembly, and relative-permittivity derivative value
  assembly.
- `src/epcsaft/equilibrium.py`: phase state construction and neutral
  bubble/dew candidate evaluation.
- `src/epcsaft/reactive_speciation.py`: activity-fixed-point state evaluation,
  scalar binary activity grid evaluation, and structured failure initial
  composition normalization.
- `src/epcsaft/runtime.py`: native extension path, CppAD metadata, and Ceres
  metadata import probes.
- `src/epcsaft/__main__.py`: package/native import failure reporting.
- `src/epcsaft/_optional_backends/ipopt.py`: optional dependency import
  diagnostics, status parsing, native seed failure recording, and stale moved
  `_core` imports.

## Validation

- `uv run python -m compileall -q src/epcsaft`: pass.
- `uv run python run_pytest.py tests/api -q`: pass, 165 passed and 1 skipped.
- `uv run python run_pytest.py tests/native/runtime -q`: pass, 21 passed.

## Follow-Up Boundary

Do not convert the two remaining broad catches without a product decision on
batch/sweep failure semantics. They are not generic dodge paths; they are public
result-mode boundaries.
