# Native Ipopt Test Audit

Date: 2026-05-16
Branch: `codex/native-ipopt-derivative-gates`
Plan: `docs/superpowers/plans/2026-05-16-native-ipopt-derivative-gates.md`

## Current Evidence

- Test inventory command: `uv run python -m pytest tests --collect-only -q`
  - Result: `531 tests collected in 2.22s`
  - Tracked Python test files: 144
- Wrapper inventory command: `uv run python run_pytest.py --all --collect-only -q`
  - Result: `531 tests collected in 2.50s`; wrapper wall time `3.518s`.
- Full duration command: `uv run python run_pytest.py --all -q --durations=30`
  - Result: failed; `496 passed, 8 failed, 27 skipped` in `239.64s`; wrapper wall time `240.664s`.
- Quick validation command: `uv run python scripts/dev/validate_project.py quick`
  - Result after the text gate landed: `32 passed in 20.48s`
  - This is comfortably under the 10 minute quick-gate target.
- Docs validation command: `uv run python scripts/dev/validate_project.py docs`
  - Result: Sphinx HTML build passed.
- Text gate command: `uv run python scripts/dev/check_text_gates.py`
  - Result: passed.
- Test hygiene command: `uv run ruff check tests --select F401,F841,ARG001,ARG002`
  - Initial result: 49 test-only issues.
  - Cleanup result: all checks passed.

## Current Test Slice Map

- `generic` / default quick slice: fast public API, representative runtime, selected regression, selected equilibrium, repo workflow tests.
- `confidence`: quick slice plus a few native confidence checks.
- `equilibrium-confidence`: Khudaida electrolyte validation confidence checks.
- `equilibrium-api`: public equilibrium and reactive API checks.
- `runtime`, `api`, `native`: focused subsystem slices.
- `profile` and `profile-full`: opt-in runtime profiling tests.
- `all`: full historical test tree.

## Slow Or Opt-In Areas

- Slowest observed tests from the full duration run:
  - `23.21s` `tests/equilibrium/electrolyte/test_electrolyte_lle_smokes.py::test_electrolyte_lle_direct_feed_reports_production_solver_derivatives`
  - `20.99s` `tests/api/reactive/test_staged_reactive_route_not_production.py::test_explicit_staged_kind_remains_separate_from_production_reactive_lle`
  - `17.82s` `tests/equilibrium/electrolyte/test_electrolyte_lle_smokes.py::test_one_salt_smoke_reports_production_solver_derivatives`
  - `15.10s` `tests/workflows/paper_validation/test_rezaee_2026_paper_validation.py::test_rezaee_source_backed_paper_validation_generates_pre_surrogate_rows`
  - `10.06s` `tests/workflows/benchmarks/test_benchmark_reactive_regression.py::test_reactive_regression_benchmark_schema_for_one_case`
  - `9.45s` `tests/equilibrium/electrolyte/test_electrolyte_lle_solver_contracts.py::test_electrolyte_lle_accepts_legacy_option_dict_before_production_failure`
  - `9.33s` `tests/equilibrium/electrolyte/test_electrolyte_lle_solver_contracts.py::test_electrolyte_lle_solver_failure_reports_production_derivatives`
  - `9.00s` `tests/equilibrium/electrolyte/test_electrolyte_lle_solver_contracts.py::test_experimental_coupled_density_lle_option_is_reported_without_changing_default_gate`
  - `8.53s` `tests/equilibrium/electrolyte/test_electrolyte_lle_solver_contracts.py::test_electrolyte_lle_best_effort_reports_production_derivatives`
  - `5.35s` `tests/workflows/benchmarks/test_benchmark_literature_suite.py::test_literature_benchmark_script_runs_and_writes_json`
- `tests/profile/**` is opt-in through `EPCSAFT_RUN_PERF`.
- MEA regression literature tests are opt-in through `EPCSAFT_RUN_MEA_TABLE2_REGRESSION` or `EPCSAFT_RUN_MEA_REGRESSION`.
- Khudaida matrix and digitized-data validation tests use explicit skip gates for unavailable local validation assets.
- Hubach hard-case electrolyte LLE tests use opt-in legacy-candidate gating.
- Full-suite duration evidence is being collected separately; keep this audit updated with the slowest individual tests before Task 2 is considered complete.

## Unnecessary Or Weak Coverage Found

- Test import and unused mock-argument bloat was mechanical and removed.
- Several tests still assert legacy missing-status behavior instead of implemented derivative coverage.
- Several tests still protect fallback diagnostics, debug skip controls, or legacy custom solver routes.
- Some Ceres and Ipopt tests currently skip when optional backends are absent; this conflicts with the new required native dependency direction and must be tightened after build gates are updated.
- `tests/native/equilibrium/test_reactive_phase_equilibrium_residual_jacobian.py` still uses source-perturbation comparison language and should be rewritten around analytical or CppAD Jacobian evidence.

Full-suite failures from the duration run that align with planned cleanup or next implementation:

- `tests/native/contracts/test_equilibrium_native_contracts.py::test_public_electrolyte_lle_reports_unavailable_solver_derivatives`
- `tests/native/contracts/test_equilibrium_native_contracts.py::test_native_electrolyte_lle_residual_evaluator_reports_unavailable_derivatives`
- `tests/native/contracts/test_equilibrium_native_contracts.py::test_native_electrolyte_lle_residual_evaluator_rejects_auto_without_autodiff`
- `tests/native/contracts/test_equilibrium_native_contracts.py::test_package_runtime_has_no_external_optimizer_dependency_or_imports`
- `tests/native/contracts/test_equilibrium_native_contracts.py::test_public_equilibrium_does_not_expose_python_backend_tokens`
- `tests/native/cppad/test_cppad_lle_derivatives.py::test_associating_neutral_lle_solves_without_numerical_derivative_derivatives`
- `tests/native/equilibrium/test_reactive_phase_equilibrium_ceres_solver.py::test_reactive_phase_native_ceres_solver_rejects_unusable_solution`
- `tests/workflows/repo/test_dependency_issue_triage.py::test_issue_number_and_url_resolve_to_same_gh_issue_view_call`

The last two failures were caused or exposed by this cleanup slice and were fixed before commit. The remaining failures are gate-aligned evidence for the upcoming dependency and solver-route overhaul.

## Dependency And Solver Gate Gaps

- `pyproject.toml` still includes SciPy in the test dependency group.
- `analyses/paper_validation/application/2026_rezaee/scripts/rezaee_reactive_equilibrium_fit.py` still imports `scipy.optimize.least_squares`.
- `src/epcsaft/_optional_backends/ipopt.py` still owns a Python `cyipopt` residual-minimization route.
- Public equilibrium options still expose `newton` and `ipopt` through Python-side routing.
- Native and Python equilibrium paths still contain custom bracketing, bisection, golden-section, retry, and fallback behavior.
- Regression still exposes native least-squares compatibility paths that must be replaced by Ceres-only production routes.

## Required Cleanup Still Open

- Delete or rewrite tests that only protect legacy missing-status behavior.
- Add passing tracked gates for no SciPy package/dev/test dependency after SciPy is removed.
- Add passing tracked gates for no Python production solve loop after native Ipopt routes exist.
- Add passing tracked gates for no Eigen nonlinear optimizer route while still allowing Eigen linear algebra.
- Move any slow scientific matrix coverage that is not already opt-in out of the quick gate.
- Replace optional-backend skip behavior with required-backend validation once Ceres, CppAD, and native Ipopt build gates are in place.

## Task 2 Status

Completed in the first Task 2 slice:

- Collected test inventory.
- Collected full-suite duration evidence and slowest-test list.
- Confirmed current quick-gate duration.
- Removed mechanical unused-import and unused-argument test bloat.
- Rewrote one stale derivative-free LLE assertion to require the current Ceres plus CppAD-implicit route.
- Added this tracked audit artifact.

Not complete yet:

- Deleting or rewriting legacy status-only tests.
- New strict dependency and solver ownership gates that require Task 3 and native Ipopt implementation to pass.
