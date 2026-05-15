# T006 Python Route and Benchmark Coverage

## Result

Done. The electrolyte LLE Python route now exposes the native Ceres residual solve and production derivative diagnostics, and the public benchmark/API tests assert the issue #116 production route instead of the old predictive route labels.

## Source Evidence

- `src/epcsaft/native/epcsaft_equilibrium.cpp` reports separated neutral fugacity, ionic equilibrium, material balance, phase charge, scaled solver, and unscaled solver residual norms from the accepted solved state.
- `src/epcsaft/equilibrium.py` failure diagnostics now name the Ceres residual solve and `cppad_implicit` derivative route instead of the old transformed Newton label.
- `tests/equilibrium/electrolyte/test_distributed_ion_lle_production_solver.py` covers the mixed-salt distributed-ion Ascani-style fixture.
- `tests/equilibrium/electrolyte/test_salting_out_lle_benchmark.py` covers a repo-contained quaternary water/butanol/NaCl salting-out fixture.
- `tests/api/equilibrium/test_electrolyte_lle_problem_production_route.py` proves the generic `ElectrolyteLLEProblem` path routes to the native Ceres production solver.

## Validation

- `uv run python run_pytest.py tests/equilibrium/electrolyte/test_distributed_ion_lle_production_solver.py tests/equilibrium/electrolyte/test_salting_out_lle_benchmark.py tests/api/equilibrium/test_electrolyte_lle_problem_production_route.py -q`: pass, 3 tests
- `uv run python run_pytest.py tests/equilibrium/electrolyte/test_electrolyte_lle_solver_contracts.py -q`: pass, 10 tests
- `uv run python run_pytest.py tests/equilibrium/electrolyte -q`: pass, 37 passed, 7 skipped
- `uv run python run_pytest.py tests/api/equilibrium/test_electrolyte_lle_problem_production_route.py -q`: pass, 1 test
- `uv run python run_pytest.py tests/native/equilibrium -q`: pass, 22 tests
- `git diff --check`: pass

## Guard Search

The route-guard search still finds historical intake notes, #117 staged-route work, neutral LLE tests, IPOPT seed naming, and explicit tests asserting the old route is absent from accepted electrolyte diagnostics. It no longer identifies an accepted electrolyte LLE production result using the old route labels.

## Next Task

Run T007 issue #116 Stage 10 validation and targeted guard review.
