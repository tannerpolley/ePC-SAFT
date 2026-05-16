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
- Full-suite rerun command after the cleanup slices: `uv run python run_pytest.py --all -q`
  - Result: `524 passed, 21 skipped` in `164.09s`; wrapper wall time `165.4s`.
- Quick validation command: `uv run python scripts/dev/validate_project.py quick`
  - Result after the native dependency gate landed: `35 passed in 21.53s`
  - This is comfortably under the 10 minute quick-gate target.
- Ceres/CppAD validation command: `uv run python scripts/dev/validate_project.py ceres-cppad`
  - Result after the native dependency gate landed: `4 passed in 2.97s`; wrapper completed after an incremental full-profile native build.
- Package boundary command: `uv run python scripts/dev/build_dist.py`
  - Result after the Ceres dependency gate landed: wheel smoke passed and the built wheel had no vendored Ceres development artifacts.
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
- NumPy testing imports were removed from tests because the lazy `np.testing` import path can stall Windows validation before the thermodynamic assertion runs.
- Several tests still assert legacy missing-status behavior instead of implemented derivative coverage.
- Several tests still protect fallback diagnostics, debug skip controls, or legacy custom solver routes.
- Ceres and CppAD are now required by the local dev script, package build backend, and CMake configure gate. Ipopt remains a system-dependency opt-in until the adapter is implemented.
- Native derivative tests now avoid perturbation-derived oracle checks; the reactive-phase Jacobian test is rewritten around CppAD backend, shape, finite-value, and analytical row-contract evidence.
- The public implicit-sensitivity placeholder helper was removed. Reactive speciation tests now require real analytical/CppAD-backed implicit solve payloads and raising behavior for unsupported implicit derivative requests.
- The legacy reactive-speciation `jacobian_backend="autodiff"` selector was removed; option tests now reject it and explicit derivative requests are limited to analytical/CppAD semantics.
- The public `EquilibriumOptions.jacobian_backend` selector now rejects the same legacy backend spelling, keeping phase-equilibrium derivative requests aligned to analytical/CppAD semantics.
- Stale native electrolyte LLE contract tests were rewritten from missing-derivative assertions to the implemented Ceres plus CppAD-implicit residual/Jacobian surface.
- The reactive phase diagnostic extent helper no longer uses NumPy's least-squares convenience path; it uses a direct library linear solve for the small stoichiometric normal system.
- No-reaction, failed, and best-effort reactive speciation paths now omit implicit-sensitivity payloads when the native route has no reaction-constant sensitivity matrix instead of returning placeholders or raising during diagnostic normalization.
- The MIAC electrolyte fixture check now uses a strict approximate comparison instead of exact binary float equality.
- A repo workflow gate now scans native C++ sources for Ceres numeric-diff APIs and legacy perturbation route tokens.

The failure list from the initial full-duration run has been retired. Each listed node now passes individually after the dependency, contract, and derivative-surface cleanup slices:

- `tests/native/contracts/test_equilibrium_native_contracts.py::test_package_runtime_has_no_external_optimizer_dependency_or_imports`
- `tests/native/cppad/test_cppad_lle_derivatives.py::test_associating_neutral_lle_solves_with_cppad_implicit_derivatives`
- `tests/native/equilibrium/test_reactive_phase_equilibrium_ceres_solver.py::test_reactive_phase_native_ceres_solver_rejects_unusable_solution`
- `tests/workflows/repo/test_dependency_issue_triage.py::test_issue_number_and_url_resolve_to_same_gh_issue_view_call`

## Dependency And Solver Gate Gaps

- The first Task 3 slice removed the legacy numerical package from the test dependency group.
- The first Task 3 slice deleted the legacy Rezaee fitting script and its generated fit outputs.
- The second Task 3 slice removed the Python IPOPT adapter, added native system Ipopt discovery, and added doctor/build status reporting. The native Ipopt adapter is present; public route wiring is currently limited to explicit homogeneous ideal reactive speciation.
- The third Task 3 slice made Ceres and CppAD mandatory native dependencies for dev-script, package-backend, and CMake builds, excluded vendored Ceres install rules from package artifacts, and validated the actual local extension with Ceres enabled.
- Public `EquilibriumOptions` now accepts only `auto` and explicit `ipopt`; `ReactiveSpeciationOptions` has the same public solver selector shape.
- Native equilibrium paths still contain custom bracketing, bisection, golden-section, retry, and fallback behavior, but the Python neutral bubble/dew scalar route has been removed from accepted public execution.
- Regression public pure-neutral fitting now defaults to Ceres and rejects the old native least-squares backend. The private generic Eigen least-squares binding and perturbation-Jacobian route are removed; generic supported production fits go through Ceres, while associating/MEA benchmark helpers are residual scorers until native Ceres derivative coverage exists.

## Required Cleanup Still Open

- Delete or rewrite tests that only protect legacy missing-status behavior.
- Add passing tracked gates for no legacy numerical package/dev/test dependency after dependency cleanup. Done in the Task 3 dependency slices.
- Add passing tracked gates for no Python production solve loop after native Ipopt routes exist.
- Add passing tracked gates for no Eigen nonlinear optimizer route while still allowing Eigen linear algebra.
- Move any slow scientific matrix coverage that is not already opt-in out of the quick gate.
- Continue replacing optional-backend skip behavior with required-backend validation as solver routes move to Ipopt and regression routes become Ceres-only.

## Task 2 Status

Completed in the first Task 2 slice:

- Collected test inventory.
- Collected full-suite duration evidence and slowest-test list.
- Confirmed current quick-gate duration.
- Removed mechanical unused-import and unused-argument test bloat.
- Rewrote one stale derivative-free LLE assertion to require the current Ceres plus CppAD-implicit route.
- Added this tracked audit artifact.
- Removed NumPy testing import usage from the test tree to keep Windows validation deterministic.

Not complete yet:

- Deleting or rewriting legacy status-only tests.
- New strict dependency and solver ownership gates that require Task 3 and native Ipopt implementation to pass.
