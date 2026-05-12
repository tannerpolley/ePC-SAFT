T027 receipt: native thermodynamic wrapper routing
===================================================

Decision
--------

`fit_reactive_electrolyte_parameters(..., backend="native")` now attempts the
production native thermodynamic regression route for the supported slice:

- reactive speciation rows only;
- ideal mole-fraction reaction standard states;
- linear speciation targets;
- reaction `logK` / reaction equilibrium constant parameters.

Unsupported high-level cases remain on the existing explicit native residual
record boundary. `backend="python_compat"` remains the only Python optimizer
loop and is still labeled non-production/debug compatibility in diagnostics.

Implementation
--------------

- Added a conservative serializer from `ReactiveElectrolyteRegressionContext`
  to the native `fit_native_thermo_regression(...)` request payload.
- Added normalized `ReactiveRegressionFitResult` wrapping for the native
  thermodynamic fit result, with diagnostics `backend="native_thermo"`,
  `python_optimizer=False`, and `Backend_unavailable_jacobian=False`.
- Added route coverage for a supported speciation/logK batch.
- Updated docs to distinguish the production native thermodynamic slice from
  the residual-record compatibility boundary and Python compatibility loop.

Validation
----------

- `uv run python run_pytest.py tests/api/test_reactive_regression.py tests/workflows/test_benchmark_reactive_regression.py -q`
  - pass: 20 passed in 74.08s
- `uv run python scripts/benchmark_native_ceres_thermo_regression.py --warmup 1 --repeat 3`
  - pass: `reactive_speciation_logk_implicit`, median 0.347 ms,
    status `backend_unavailable`, backend `backend_unavailable`, derivative
    `implicit`, fixed native thermodynamic objective reported
- `uv run ruff check src/epcsaft/reactive_regression.py tests/api/test_reactive_regression.py`
  - pass
- `uv run black --check src/epcsaft/reactive_regression.py tests/api/test_reactive_regression.py`
  - pass

Remaining gap
-------------

This task does not close issue #53 by itself. Earlier receipts still show that
production Ceres sensitivities for Born-SSM+DS `d_born`/`f_solv` and reactive
electrolyte bubble pressure rows are not implemented; they report
`backend_unavailable` rather than using backend-unavailable shortcuts.

