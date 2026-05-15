# T026 Association Sensitivity Repair

## Result

Done. The active association phase-state fugacity sensitivity path now matches projected solved-state perturbation behavior, and the accepted electrolyte LLE Ceres solve uses the transformed residual Jacobian with `cppad_implicit` provenance.

## Source Evidence

- `src/epcsaft/native/epcsaft_ares.cpp` now keeps `thermo.den` on the CppAD tape when building association composition-response zeta derivatives. This preserves the density leg needed by pressure-state association fugacity sensitivities.
- `src/epcsaft/native/epcsaft_equilibrium.cpp` removed the local residual slope Ceres Jacobian path and copies `electrolyte_residual_jacobian_row_major(...)` into the Ceres residual block.
- Accepted electrolyte LLE diagnostics now report `jacobian_backend = cppad_implicit` and `derivative_backend = cppad_implicit`.

## Test Evidence

- `tests/native/cppad/test_phase_state_sensitivities.py` now checks projected pressure-state perturbations for both the neutral and active association/ionic fixtures.
- `tests/native/equilibrium/test_electrolyte_lle_ceres_solver.py` now fails if accepted diagnostics contain `local_residual_slope`.

## Validation

- `uv run python scripts/dev/build_epcsaft.py --profile full`: pass
- `uv run python run_pytest.py tests/native/cppad/test_phase_state_sensitivities.py tests/native/equilibrium/test_electrolyte_lle_ceres_solver.py tests/native/equilibrium/test_electrolyte_lle_residual_jacobian.py -q`: pass, 5 tests
- `uv run python run_pytest.py tests/native/equilibrium tests/native/cppad -q`: pass, 52 tests
- `git diff --check`: pass

## Next Task

Continue issue #116 with T006: TPD/g-hat role guards, generic Python route proof, and distributed-ion plus salting-out benchmark coverage.
