# T022 Completion Audit

## Decision

Rejected for full issue #53 completion.

The branch now contains a real native Ceres thermodynamic regression slice and a backend-unavailable debug gate, but it still does not satisfy the full GitHub issue #53 / comment 4424893353 acceptance criteria.

## What Is Now Real

- Native C++ owns a supported Ceres hot loop for reactive-speciation rows with reaction logK parameters.
- That supported slice uses analytic/implicit derivatives and reports no backend-unavailable use.
- Python serializes the supported native thermodynamic fit request once; C++ applies parameter vectors and owns row evaluation during Ceres iterations.
- Unsupported derivative combinations return ``backend_unavailable`` rather than silent backend-unavailable fallback.
- Explicit backend-unavailable chemical-equilibrium/electrolyte-LLE diagnostics require ``EPCSAFT_ALLOW_DERIVATIVE_BACKEND_DEBUG=1``.

## Blockers To Closing Issue #53

1. Reactive electrolyte bubble rows are native-evaluated but not yet differentiable for Ceres.
   Evidence: `src/epcsaft/native/regression/thermo_regression.cpp` supports `reactive_electrolyte_bubble` in the row evaluator, but `thermo_derivative_supported(...)` rejects row modes other than `reactive_speciation`.

2. Born/SSM+DS and `k_ij` parameter kinds can be applied to native mixtures, but production Ceres derivatives are not implemented for them.
   Evidence: `apply_native_thermo_parameters(...)` handles `born_radius`/`born_diameter` and `binary_interaction`/`k_ij`, while `thermo_derivative_supported(...)` currently accepts only `reaction_equilibrium_constant` / `log_equilibrium_constant`.

3. Activity-coupled and concentration-coupled reactive speciation still lack production analytic/autodiff/implicit derivatives.
   Evidence: native chemical equilibrium now returns `backend_unavailable` for these cases unless the backend-unavailable debug gate is set.

4. The public high-level reactive regression helper still has a compatibility Python optimizer path.
   Evidence: `src/epcsaft/reactive_regression.py` still contains `Backend_unavailable_jacobian(...)`, `backend="python_compat"`, and compatibility diagnostics. This is acceptable as debug/compatibility but not enough to say all production fitting has moved to the new native thermodynamic Ceres path.

5. The default local build has Ceres disabled.
   Evidence: `scripts/doctor.py` reports `native_dependency_ceres: enabled=OFF`; the default benchmark reports `backend_unavailable`. T019 did prove the Ceres-enabled temp build can run the supported native hot loop with objective decrease.

## Validation Evidence

- `uv run python scripts/build_epcsaft.py --build-only --parallel 10`: passed.
- `uv run python scripts/validate_project.py quick`: passed, `23 passed in 15.64s`.
- `uv run python scripts/validate_project.py docs`: passed.
- `uv run python run_pytest.py tests/api/test_Backend_unavailable_debug_policy.py tests/api/test_reactive_speciation.py tests/api/test_runtime.py tests/native/test_chemical_equilibrium_native.py tests/native/test_equilibrium_native_contracts.py tests/equilibrium/test_lle.py -q`: passed, `118 passed, 1 skipped`.
- `uv run python run_pytest.py tests/api/test_runtime.py tests/workflows/test_benchmark_native_regression.py tests/workflows/test_benchmark_reactive_regression.py -q`: passed, `52 passed`.
- `uv run python scripts/benchmark_native_regression.py --case native_reactive_born_kij_tiny --warmup 0 --repeat 1`: passed.
- `uv run python scripts/benchmark_native_ceres_thermo_regression.py --warmup 1 --repeat 3`: passed on default build with honest `backend_unavailable`.

## Required Next Tasks

Continue with a new implementation tranche before claiming issue #53 complete:

1. Add native implicit/analytic parameter sensitivities for reactive electrolyte bubble pressure targets.
2. Add Ceres derivative support for Born/SSM+DS and `k_ij` parameters, or route them through CppAD-templated residuals.
3. Add activity/concentration standard-state derivatives for reactive speciation, or keep them unsupported and explicitly out of the production slice.
4. Wire high-level production reactive regression calls to the native thermodynamic Ceres fit where the supported slice applies.
5. Add a Ceres-enabled CI/build check or documented release verification command that proves `optimizer_backend=ceres`, `native_hot_loop=true`, `python_objective_used=false`, and objective decrease.

