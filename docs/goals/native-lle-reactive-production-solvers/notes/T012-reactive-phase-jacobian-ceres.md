# T012 Reactive Phase Jacobian And Ceres Route

## Result

Implemented the issue #117 Stages 4-6 native Jacobian and Ceres solve foundation for coupled reactive phase equilibrium.

## What Changed

- Added a coupled residual Jacobian for `log_phase_species_amounts`.
- The Jacobian uses analytic chain-rule blocks around the existing CppAD/implicit phase-state fugacity composition sensitivities.
- Added `_solve_reactive_phase_equilibrium_native(...)` through the pybind layer.
- Added a native Ceres trust-region residual solve with diagnostics for solver backend, method, trust-region strategy, linear solver, termination, costs, Jacobian backend, derivative backend, and solved-state sensitivity backend.
- Added native tests for the coupled residual Jacobian and Ceres route.

## Solver Boundary

This slice proves the native coupled residual/Jacobian/Ceres foundation directly. It does not yet make `ReactivePhaseEquilibriumProblem.solve(...)`, `kind="reactive_lle"`, or `kind="reactive_electrolyte_lle"` use the coupled production route. That public API work remains T013.

## Validation

- `uv run python scripts/dev/build_epcsaft.py --profile full`: pass during iteration.
- `uv run python run_pytest.py tests/native/equilibrium/test_reactive_phase_equilibrium_ceres_solver.py tests/native/equilibrium/test_reactive_phase_equilibrium_residual_jacobian.py -q`: pass, 2 tests during iteration.
- `uv run python scripts/dev/build_epcsaft.py --clean --enable-ceres --enable-cppad`: pass.
- `uv run python run_pytest.py tests/native/equilibrium/test_reactive_phase_equilibrium_ceres_solver.py tests/native/equilibrium/test_reactive_phase_equilibrium_residual_jacobian.py -q`: pass, 2 tests after exact clean build.
- `uv run python run_pytest.py tests/native/equilibrium -q`: pass, 26 tests.

## Next

T013 should route the generic Python production APIs to the native coupled solve, add neutral and electrolyte reactive LLE benchmark tests, and keep the staged helper explicitly separate.
