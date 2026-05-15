# T013 Public API And Benchmarks Receipt

## Result

Done. Production reactive LLE routes now call the native coupled reactive phase-equilibrium solver instead of the explicit staged helper.

## Changes

- Added a Python adapter for `_solve_reactive_phase_equilibrium_native(...)` that builds the coupled request from public species, balance, reaction, feed, initial-phase, and option inputs.
- Routed `ReactivePhaseEquilibriumProblem.solve(...)`, `kind="reactive_lle"`, and `kind="reactive_electrolyte_lle"` to the native coupled production route.
- Kept `kind="reactive_staged"` and `reactive_staged_equilibrium(...)` as the explicit staged compatibility route.
- Added neutral reactive LLE and reactive electrolyte LLE benchmark tests with checks for reaction, phase, ionic, material, element, charge, phase-distance, composition, phase-amount, reaction-extent, solver, Jacobian, and derivative diagnostics.
- Updated the prior staged reactive LLE regression so it asserts only the explicit staged route.

## Validation

- `uv run python run_pytest.py tests/equilibrium/reactive/test_reactive_lle.py tests/equilibrium/reactive/test_reactive_lle_coupled_solver.py tests/equilibrium/reactive/test_reactive_electrolyte_lle_coupled_solver.py tests/api/reactive/test_reactive_phase_equilibrium_problem_routes_native.py tests/api/reactive/test_staged_reactive_route_not_production.py -q`: pass, 8 tests.
- `uv run python run_pytest.py tests/equilibrium/reactive -q`: pass, 4 tests.
- `uv run python run_pytest.py tests/api/reactive -q`: pass, 65 passed and 1 skipped for optional `cyipopt`.

## Next

T014 owns the clean Ceres+CppAD build, native/reactive/API validation ladder, docs validation, diff check, and route-guard audit.
