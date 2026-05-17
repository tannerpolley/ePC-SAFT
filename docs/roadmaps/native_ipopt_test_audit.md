# Native Ipopt Test Audit

Date: 2026-05-16
Branch: `codex/native-ipopt-derivative-gates`
Plan: `docs/superpowers/plans/2026-05-16-native-ipopt-derivative-gates.md`

## Current Evidence

- Test inventory command: `uv run python -m pytest tests --collect-only -q`
  - Result: `531 tests collected in 2.22s`
  - Tracked Python test files: 144
- Wrapper inventory command: `uv run python run_pytest.py --all --collect-only -q`
  - Latest result after duplicate staged reactive LLE prune: `527 tests collected in 1.44s`; wrapper wall time `3.1s`.
- Full duration command: `uv run python run_pytest.py --all -q --durations=30`
  - Result: failed; `496 passed, 8 failed, 27 skipped` in `239.64s`; wrapper wall time `240.664s`.
- Full-suite rerun command after the cleanup slices: `uv run python run_pytest.py --all -q`
  - Result: `524 passed, 21 skipped` in `164.09s`; wrapper wall time `165.4s`.
- Quick validation command: `uv run python scripts/dev/validate_project.py quick`
  - Result after the native dependency gate landed: `35 passed in 21.53s`
  - This is comfortably under the 10 minute quick-gate target.
  - Latest result after the fixed-temperature pressure-route and workflow-gate slices: `41 passed in 4.21s`;
    wrapper command completed in `6.6s`.
  - Latest result after the CppAD capability and staged-test cleanup slices: `40 passed in 4.44s`;
    wrapper command completed in `7.1s`.
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
  - `23.21s` `tests/equilibrium/electrolyte/test_electrolyte_lle_smokes.py::test_electrolyte_lle_direct_feed_reports_current_ceres_derivatives`
  - `20.99s` `tests/api/reactive/test_staged_reactive_route_not_production.py::test_explicit_staged_kind_remains_separate_from_production_reactive_lle`
  - `17.82s` `tests/equilibrium/electrolyte/test_electrolyte_lle_smokes.py::test_one_salt_smoke_reports_current_ceres_derivatives`
  - `15.10s` `tests/workflows/paper_validation/test_rezaee_2026_paper_validation.py::test_rezaee_source_backed_paper_validation_generates_pre_surrogate_rows`
  - `10.06s` `tests/workflows/benchmarks/test_benchmark_reactive_regression.py::test_reactive_regression_benchmark_schema_for_one_case`
- `tests/profile/**` is opt-in through `EPCSAFT_RUN_PERF`.
- MEA regression literature tests are opt-in through `EPCSAFT_RUN_MEA_TABLE2_REGRESSION` or `EPCSAFT_RUN_MEA_REGRESSION`.
- Khudaida matrix and digitized-data validation tests use explicit skip gates for unavailable local validation assets.

## Unnecessary Or Weak Coverage Found

- Test import and unused mock-argument bloat was mechanical and removed.
- Hubach electrolyte LLE status-only continuation tests were deleted, and removed-option checks now run as fast validation-path tests instead of opt-in skips.
- NumPy testing imports were removed from tests because the lazy `np.testing` import path can stall Windows validation before the thermodynamic assertion runs.
- The reactive-speciation runtime capability no longer exposes a derivative-gap status string; tests use the structured
  supported and route-gated standard-state lists instead.
- Remaining weak coverage is concentrated in opt-in diagnostics, debug skip controls, and duplicate route-gate docs;
  continue deleting tests that only protect status/debug surfaces.
- Ceres and CppAD are now required by the local dev script, package build backend, and CMake configure gate. Ipopt remains a system-dependency opt-in until the adapter is implemented.
- Native derivative tests now avoid shifted-source oracle checks; the reactive-phase Jacobian test is rewritten around CppAD backend, shape, finite-value, and analytical row-contract evidence.
- The electrolyte LLE confidence report no longer writes parameter-shift sensitivity CSV or plot artifacts; it keeps benchmark, continuation, oracle, stress, and residual/error outputs.
- The public implicit-sensitivity placeholder helper was removed. Reactive speciation tests now require real analytical/CppAD-backed implicit solve payloads and raising behavior for unsupported implicit derivative requests.
- The legacy reactive-speciation `jacobian_backend="autodiff"` selector was removed; option tests now reject it and explicit derivative requests are limited to analytical/CppAD semantics.
- The public `EquilibriumOptions.jacobian_backend` selector now rejects the same legacy backend spelling, keeping phase-equilibrium derivative requests aligned to analytical/CppAD semantics.
- Stale native electrolyte LLE and reactive-phase tests that protected accepted Ceres equilibrium solves were deleted. Residual-surface derivative tests remain as private diagnostic coverage until native Ipopt route builders own production equilibrium.
- Old PR #126 and issue-specific Ceres-equilibrium handoff documents were removed from active docs, and literature benchmark metadata now points blocked Ascani LLE/reactive-phase cases at the native Ipopt gate plan.
- Public route-pending errors now describe the current Ipopt ownership requirement without naming retired solver routes as compatibility context.
- Runtime reactive-regression capability labels now describe residual-evaluation contexts instead of naming Python orchestration as a solver backend.
- The unreferenced tracked LaTeX backup `docs/latex/equations_old.tex` was deleted; `docs/latex/equations.tex` remains the equation source of truth.
- Completed, unreferenced JetBrains cleanup plan artifacts were removed now that the repo-owned script exists.
- Remaining unreferenced stale handoff/planning artifacts under `docs/handoffs/` and `docs/plans/` were removed so the native Ipopt gate plan is the active implementation handoff.
- The public reactive-electrolyte fit helper now validates fit requests and raises until native Ceres owns that optimizer with exact derivatives; the retained public path is objective/residual evaluation, not a manufactured fit result.
- The unused private `_solve_equilibrium_native` pybind wrapper and Python payload adapter were removed; public equilibrium routes remain gated to native Ipopt builders, and private residual-surface bindings stay separate.
- The `_core` pybind module now disables pybind11 release extras and MSVC optimization for the large binding translation unit; native thermodynamic objects remain Release-optimized.
- The reactive phase diagnostic extent helper no longer uses NumPy's least-squares convenience path; it uses a direct library linear solve for the small stoichiometric normal system.
- No-reaction, failed, and best-effort reactive speciation paths now omit implicit-sensitivity payloads when the native route has no reaction-constant sensitivity matrix instead of returning placeholders or raising during diagnostic normalization.
- The MIAC electrolyte fixture check now uses a strict approximate comparison instead of exact binary float equality.
- A repo workflow gate now scans native C++ sources for Ceres non-exact derivative APIs and legacy shifted-source route tokens.
- A repo workflow gate now scans public Python solver surfaces for external optimizer/root-loop calls, so Python facades
  remain request/result surfaces instead of owning production solve algorithms.
- Duplicate distributed-ion electrolyte LLE route-pending checks were removed; mixed-salt public route-gate coverage
  remains in the route-specific solver-contract test file.
- The retired tracked `docs/goals/**` GoalBuddy boards were removed from active docs. The only referenced benchmark
  evidence was condensed into `docs/roadmaps/native_associating_binary_ceres_benchmark.md`, and the text gate now blocks
  the stale derivative-proof phrases that appeared only in historical board state.
- Duplicate public route-pending checks were pruned from the native equilibrium contract file. Route-specific public
  tests still own TP flash, stability, and electrolyte LLE route-gate behavior, while the native contract file keeps the
  native residual and dependency-boundary checks.
- The standalone salting-out LLE route-pending benchmark test was removed because the smoke tests already cover strict
  electrolyte LLE route gating and the solver-contract tests cover mixed-salt route gates. Its typed problem
  fixture now lives with the `ElectrolyteLLEProblem` problem-object test.
- A duplicate electrolyte LLE explicit-Ipopt route-pending test was removed. The retained native-route request test
  already sets `solver_backend="ipopt"` and verifies the exact route payload before the local no-Ipopt gate.
- The no-op public `EquilibriumOptions.include_phase_diagnostics` switch was removed after confirming native route
  payloads do not read it. Current option acceptance stays covered by exact field-set and unknown-key tests.
- A duplicate TP-flash route-pending test for the retired stability-precheck option was removed. The route-specific
  TP-flash test still verifies the native route payload before the local no-Ipopt gate.
- The no-reaction mixed reactive-regression objective status test was removed. Target-family accounting is covered in
  diagnostics tests, and retained regression setup tests now exercise a real nonideal native derivative-block gate.
- `equilibrium_curve(...)` no longer carries accepted phase splits forward as Python-level seeds; curve points use
  route-owned canonical initial points, and user-supplied phase seeds are rejected before route execution.
- Supported generic native Ceres regression routes no longer treat `max_nfev=1` as an initial-residual shortcut. They
  run through Ceres, reject nonpositive evaluation limits before native dispatch, and the text gate blocks the retired
  initial-evaluation optimizer message.
- The accepted Ceres regression surfaces were checked against the autodiff gate. Current production residuals depend on
  implicit density/EOS derivatives or Born/activity derivative helpers and remain `cppad_implicit`; no direct-template
  production Ceres residual is active yet. Native Ceres option setup now uses one shared helper.
- The older reactive-electrolyte bubble residual wrapper and its public result payload were removed. The retained
  diagnostic path is the native-first `ReactiveElectrolyteRegressionContext` plus
  `evaluate_reactive_regression_objective(...)` surface.
- Reactive-regression benchmarks no longer keep an in-process legacy timing baseline against the removed wrapper. The
  retained benchmark comparison path is optional external JSON baseline reporting against supported cases.
- Electrolyte/property runtime options now use `cppad` for explicit CppAD derivative requests, and the old generic
  AD option spelling is rejected outside Ceres-owned regression contexts.
- Equilibrium cookbook and downstream local-install docs now distinguish Ipopt-enabled native routes from route-gated
  stability, bubble/dew temperature, and reactive-electrolyte bubble routes.
- The public activity-coefficient contribution decomposition path now raises the package's typed `InputError` instead
  of a generic Python unsupported-operation exception, and the regression derivative table now describes neutral LLE as
  a native Ipopt route when compiled rather than a pending route. Runtime capabilities omit the unavailable
  activity-coefficient decomposition flag.
- The strict text gate now blocks generic "not implemented" wording in executable source, tests, and scripts, so
  unsupported derivative or route paths must state the positive required backend/formulation instead.
- Reactive-speciation diagnostics omit association-coupling metadata until a native route returns real coupled
  association variables; the executable text gate blocks the old association solver-status field name.
- Native Ceres regression tests now assert required Ceres/CppAD build support instead of skipping when those required
  dependencies are absent. The executable text gate blocks the retired optional-Ceres skip wording.
- Runtime build-contract and CppAD smoke tests now require enabled Ceres/CppAD native dependencies instead of accepting
  disabled or unconfigured states as valid local test outcomes.
- Native CppAD derivative helpers no longer keep disabled-payload branches, and native Ceres regression no longer keeps
  unreachable disabled-backend exception branches. The source text gate now blocks those retired strings case-insensitively.
- Runtime build metadata now uses `native_dependencies` instead of the retired optional-dependency label. Ceres and
  CppAD entries are marked as required native dependencies, and source/test/script gates block the old dependency-status
  wording.
- The stale staged TP-reactive alias guard and duplicate route-pending test were removed. Explicit
  `phase_kind` staged workflow coverage remains, and the executable text gate blocks the old alias token.
- Duplicate equilibrium capability tests that only repeated runtime metadata route-status checks were deleted. The
  retained capability tests cover derivative policy, reactive-speciation standard-state gates, and reactive
  phase-equilibrium reaction-scope metadata.
- Generic native residual-score benchmark paths now report only the native residual evaluator backend; optimizer,
  derivative, and Jacobian backend fields stay empty because no optimization or derivative evaluation is performed.
- Generic staged reactive workflow diagnostics no longer include a benchmark-specific attempt field or negative
  neutral-route status. The retained diagnostics cover staged method metadata, residual audits, derivative policy, and
  phase-route diagnostics.
- Native Ipopt route gates now use the positive `ipopt_dependency_required` status for local no-Ipopt builds, and the
  executable text gate blocks the retired build-requirement status literal in active source, tests, and scripts.
- Public Ipopt capabilities and chemical-equilibrium diagnostics no longer publish negative Hessian availability flags;
  Ipopt Hessian behavior stays internal to the adapter while public contracts require exact gradients and Jacobians.
- The ideal reactive-speciation Ipopt route now reports the positive analytical derivative backend without no-op activity
  derivative policy or Hessian-mode diagnostics.
- Native Ipopt smoke and route result payloads no longer expose Hessian-strategy or exact-Hessian-required fields; those
  choices stay inside the adapter while result payloads report exact gradient/Jacobian requirements.
- Runtime reactive phase-equilibrium capabilities no longer emit an empty negative reaction-scope list. The retained
  contract is the positive `supported_reaction_scopes` list plus cross-phase quotient metadata, and the executable text
  gate blocks the retired empty-field label.
- Runtime Ipopt dependency probing now reports `native_extension_missing` or `ipopt_probe_missing` instead of a generic
  configuration-status fallback, and the executable text gate blocks the retired fallback label.
- Equilibrium capabilities now list implemented native Ipopt routes only. Unimplemented stability, reactive bubble, and
  reactive phase-equilibrium routes remain public API contracts that fail loudly, but they are no longer advertised as
  available-false capability rows. Implemented but locally uncompiled Ipopt routes use `ipopt_dependency_required`.
- Derivative coverage matrices now represent out-of-scope rows through the existing classification/backend labels
  instead of carrying a redundant negative applicability column. Active source, tests, and scripts now block the retired
  column label.
- Electrolyte bubble options no longer expose public pressure or vapor-composition seed controls. Reactive electrolyte
  bubble sweeps no longer expose continuation controls, and reactive-regression residual contexts use each row's explicit
  initial composition/defaults without converting row-result composition, pressure, or vapor data into later route seeds.
  Bubble-pressure route initialization remains owned by the native Ipopt route builder.
- Public electrolyte LLE no longer accepts user phase seeds or exports electrolyte LLE seed-helper wrappers. The route
  builder owns the canonical initial point, while private native residual-surface tests keep explicit evaluator phase
  payloads where needed for diagnostic Jacobian coverage.
- Public neutral LLE and reactive LLE facades no longer accept user phase seeds. Direct LLE, typed LLE, reactive phase,
  and staged reactive LLE calls all use route-owned canonical initial points; public validation-script paths with
  explicit phase initializers were removed.
- Retained analysis scripts no longer pass public equilibrium phase initializers or manual solve controls. The repo
  workflow tests block those keywords while allowing private native residual/Jacobian evaluator phase data.
- Duplicate electrolyte and Hubach tests that individually protected old option names were deleted. The core LLE tests
  now assert the current `EquilibriumOptions` field set and one generic unknown-option rejection path instead of
  preserving one assertion per retired compatibility key.
- Reactive bubble seed-control and LLE problem-object surface checks were folded into exact field-set coverage.
- A duplicate staged reactive LLE test was folded into the stronger reaction-coordinate and split-diagnostic test.
- Reactive-regression row/options/result surface checks now use one compact current-surface assertion instead of three
  separate removed-field tests.
- Duplicate neutral LLE explicit-Ipopt and stability-option route-gate tests were removed. The quick gate keeps one
  representative native route-request test per active public equilibrium route, while field-set and invalid-option
  tests keep current option-surface coverage.
- Staged reactive-equilibrium diagnostics now report direct counts, labels, residuals, and minimum composition instead
  of redundant pass/reported/split status strings.
- Python equilibrium facades no longer fabricate rejected, empty-start, or accepted route-status defaults when a
  native route omits that metadata; diagnostics now preserve only non-empty route fields returned by native code.
- Reactive-regression objective-only summary payloads no longer add fit/covariance status placeholders;
  summaries keep objective metrics plus fit metadata only when a native fit result exists.
- Reactive-speciation diagnostics no longer synthesize solved density, bubble-pressure, association-coupling, or
  best-state metadata for ideal native routes; those fields now appear only when a native route returns real evidence.
- Reactive-speciation sweep input-validation failures no longer fabricate backend or selected-solver labels; those
  structured failure payloads report failure context only, and the executable text gate blocks the retired fake label.
- `EquilibriumOptions.timeout_seconds` is now wired through public neutral/electrolyte equilibrium facades into the
  native Ipopt adapter as a wall-clock option instead of being a normalized-but-unused public control.

The failure list from the initial full-duration run has been retired. Each listed node now passes individually after the dependency, contract, and derivative-surface cleanup slices:

- `tests/native/contracts/test_equilibrium_native_contracts.py::test_package_runtime_has_no_external_optimizer_dependency_or_imports`
- deleted obsolete accepted-solve coverage in `tests/native/cppad/test_cppad_lle_derivatives.py`
- deleted obsolete accepted-solve coverage in `tests/native/equilibrium/test_reactive_phase_equilibrium_ceres_solver.py`
- `tests/workflows/repo/test_dependency_issue_triage.py::test_issue_number_and_url_resolve_to_same_gh_issue_view_call`

## Dependency And Solver Gate Gaps

- The first Task 3 slice removed the legacy numerical package from the test dependency group.
- The first Task 3 slice deleted the legacy Rezaee fitting script and its generated fit outputs.
- The second Task 3 slice removed the Python IPOPT adapter, added native system Ipopt discovery, and added doctor/build status reporting. The native Ipopt adapter is present; public route wiring currently covers homogeneous ideal reactive speciation, neutral TP/LLE/bubble/dew pressure routes, electrolyte LLE, and fixed-liquid electrolyte bubble pressure when Ipopt is compiled.
- The third Task 3 slice made Ceres and CppAD mandatory native dependencies for dev-script, package-backend, and CMake builds, excluded vendored Ceres install rules from package artifacts, and validated the actual local extension with Ceres enabled.
- Public `EquilibriumOptions` now accepts only `auto` and explicit `ipopt`; `ReactiveSpeciationOptions` has the same public solver selector shape.
- Native homogeneous chemical-equilibrium solves now dispatch only to the explicit Ipopt ideal-speciation NLP. Remaining known custom bracket/root behavior is in density closure diagnostics, which stay evidence-gated rather than accepted as a general equilibrium solver route.
- A workflow gate now enforces that custom scalar root/search solver tokens are confined to the density-closure exception
  files and cannot appear in new active `src/epcsaft` solver surfaces.
- Regression public pure-neutral fitting now defaults to Ceres and rejects the old native least-squares backend. The private generic Eigen least-squares binding, shifted-source Jacobian route, and repeated-start Ceres controls are removed; generic supported production fits go through one canonical Ceres start, while associating/MEA benchmark helpers are residual scorers until native Ceres derivative coverage exists.
- Runtime derivative capabilities now list implemented production coverage only. Open derivative blockers remain in the roadmap and state-level coverage matrices, not in `epcsaft.capabilities()` as capability rows.
- Reactive electrolyte batch regression is no longer described as a production optimizer in runtime capabilities; it is a diagnostic residual context until native Ceres owns that route.

## Required Cleanup Still Open

- Delete or rewrite tests that only protect duplicate route-gate or status/debug behavior.
- Add passing tracked gates for no legacy numerical package/dev/test dependency after dependency cleanup. Done in the Task 3 dependency slices.
- Add passing tracked gates for no Python production solve loop after native Ipopt routes exist. Done for public Python solver surfaces.
- Add passing tracked gates for no Eigen nonlinear optimizer route while still allowing Eigen linear algebra. Done for native regression sources, with Ceres non-exact derivative sources gated separately.
- Move any slow scientific matrix coverage that is not already opt-in out of the quick gate. Done for the current quick slice; continue checking when adding new validation tests.
- Continue replacing optional-backend skip behavior with required-backend validation as solver routes move to Ipopt and regression routes become Ceres-only.

## Task 2 Status

Completed in the first Task 2 slice:

- Collected test inventory.
- Collected full-suite duration evidence and slowest-test list.
- Confirmed current quick-gate duration.
- Removed mechanical unused-import and unused-argument test bloat.
- Rewrote one stale derivative-absent LLE assertion to require the current Ceres plus CppAD-implicit route.
- Added this tracked audit artifact.
- Removed NumPy testing import usage from the test tree to keep Windows validation deterministic.

Still ongoing:

- Continue pruning duplicate route-gate and status/debug tests when structured route or derivative coverage already
  exists.
- New strict dependency and solver ownership gates that require Task 3 and native Ipopt implementation to pass.
