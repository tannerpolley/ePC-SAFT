# T011 Reactive Phase Residual Surface

## Result

Implemented the issue #117 Stages 1-3 residual-surface foundation for coupled reactive phase equilibrium.

## What Changed

- Added `ReactivePhaseResidualEvaluationNative` to the native equilibrium contract.
- Added `_evaluate_reactive_phase_equilibrium_residual_native(...)` through the pybind layer.
- Added `src/epcsaft/native/equilibrium/reactive_phase_equilibrium_problem.cpp` with a two-liquid coupled residual evaluator.
- Added native tests for neutral reactive LLE and reactive electrolyte LLE residual payloads.
- Updated the GoalBuddy allowed-file boundary to include `src/epcsaft/bindings.cpp`, because exposing a native residual surface requires the pybind entrypoint.

## Residual Blocks

The evaluator unpacks one `log_phase_species_amounts` vector into two nonnegative liquid phase amount/composition states. The residual payload includes:

- element or balance-row residuals,
- reaction-equilibrium residuals for both liquid phases,
- neutral interphase equilibrium residuals,
- electroneutral ion-combination residuals when charged species exist,
- per-phase charge residuals,
- phase distance and residual norms.

Reaction and phase residuals are evaluated from the same unpacked native state. The evaluator intentionally reports `solver_backend = residual_surface_only`, `ceres_accepted_solve = false`, and `jacobian_available = false`; the production Jacobian and accepted Ceres solve remain queued for T012.

## Validation

- `uv run python scripts/dev/build_epcsaft.py`: pass, fast profile.
- `uv run python run_pytest.py tests/native/equilibrium/test_reactive_phase_equilibrium_residual_surface.py -q`: pass, 2 tests.
- `uv run python run_pytest.py tests/native/equilibrium -q`: initially failed after the fast rebuild because Ceres was disabled for existing #116 tests.
- `uv run python scripts/dev/build_epcsaft.py --profile full`: pass.
- `uv run python run_pytest.py tests/native/equilibrium -q`: pass, 24 tests.

## Next

T012 should add the coupled residual Jacobian, solved-state sensitivity diagnostics, and the accepted Ceres trust-region route. The T011 residual surface is not sufficient for issue #117 completion by itself.
