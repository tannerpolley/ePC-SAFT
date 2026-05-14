# PR #56 salvage map for staged roadmap work

PR: https://github.com/tannerpolley/ePC-SAFT/pull/56

This is a reference-only disposition map for PR #56. Do not continue PR #56, merge it, base new work on it, or cherry-pick from it directly. Use it only as historical evidence while implementing the staged roadmap from fresh `main` branches.

## Safe to salvage

- The high-level problem decomposition is useful: split derivative cleanup, CppAD substrate work, Ceres build plumbing, regression API design, native regression validation, capability reporting, and package organization into separate roadmap slices rather than one broad PR.
- The PR reinforces that derivative status values need a package-wide contract. Salvage only the contract intent, not the implementation. The current roadmap should converge on explicit `backend_unavailable` behavior when analytic, CppAD, or implicit derivative coverage is absent.
- The PR identifies useful validation surfaces: runtime capability reporting, native regression benchmark smoke checks, CMake option reporting, and focused tests around native regression contracts. Rebuild those surfaces from current `main` under the corresponding roadmap issues.
- The PR shows that residual-record boundaries are a useful compatibility concept for transporting regression rows across Python/native layers. Salvage that as a boundary design idea only, not as a production regression claim.
- The PR's broad changed-file inventory is useful as a search map for future work. Treat it as an inspection index, not a source patch.

## Must not salvage

- `finite_difference` renamed to `unsupported_derivative`. Issue #58 requires finite-difference execution and terminology removal, not a semantic rename.
- `unsupported_derivative` as an executable backend. Unsupported derivative paths must fail or report `backend_unavailable`; they must not compute comparison Jacobians, fallback perturbations, or production derivatives.
- Tests that pass when the production backend is unavailable. Tests for Ceres, CppAD, and native regression must prove the supported backend path or assert a clear unavailable result in a non-production compatibility check.
- Capability claims with `native_hot_loop=True` before tests prove the native Ceres loop is actually built, selected, and executing the optimization hot loop.
- Residual-record solving labeled as production regression. Residual records can be a transport or compatibility boundary, but production regression requires real optimizer ownership, derivative policy, validation, and benchmark evidence.
- Broad unrelated goal-note churn. PR #56 added a large `docs/goals/issue-53-native-regression-production/` execution trail. Do not salvage this into staged roadmap branches unless a future issue explicitly needs a compact handoff excerpt.
- Broad production edits from the PR as a unit. The PR spans build configuration, docs, runtime reporting, Python APIs, C++ native kernels, benchmarks, tests, and historical notes. Future work must be one issue per branch from current `main`.
- Documentation that advertises narrow sketches as production-ready capability. Capability pages, README text, and benchmark output should stay conservative until the matching build and tests prove the claim.

## Maybe salvage after rewrite

- CMake Ceres/CppAD ideas: PR #56 sketches package-level `EPCSAFT_ENABLE_CERES`, `EPCSAFT_USE_SYSTEM_CERES`, `EPCSAFT_ENABLE_CPPAD`, and `EPCSAFT_USE_SYSTEM_CPPAD` options in `CMakeLists.txt`, `CMakePresets.json`, `scripts/build_epcsaft.py`, and `scripts/doctor.py`. Rebuild this under issue #68 with a clean build matrix and dependency policy.
- Runtime dependency-reporting ideas: `scripts/doctor.py`, `src/epcsaft/runtime.py`, and related tests show useful reporting concepts for native dependency availability. Rewrite under issue #68 so reporting distinguishes configured, compiled, imported, and actively used states.
- Ceres cost-function sketches: `src/epcsaft/native/regression/thermo_regression.cpp`, `src/epcsaft/native/regression/thermo_regression.h`, `src/epcsaft/benchmarks/native_ceres_thermo_regression.py`, and `tests/native/test_native_ceres_thermodynamic_regression.py` contain cost-function shape ideas. Rebuild under issues #65, #66, and #67 only after API shape, validation data, and derivative requirements are explicit.
- CppAD scalar scaffold sketches: `src/epcsaft/native/autodiff/ad_scalar.h`, `src/epcsaft/native/autodiff/ad_derivative_checks.*`, and the CppAD derivative test files show a possible substrate direction. Rewrite under issues #59, #60, #61, #62, and #63 with scalar-type ownership and equation templating planned before broad kernel edits.
- Native regression schema ideas: `src/epcsaft/native/regression/regression_types.*`, `src/epcsaft/native_regression.py`, and native regression tests suggest a typed result/schema boundary. Rebuild under issues #65, #66, and #67 with clear public API contracts and literature-backed validation.
- Benchmark harness ideas: `scripts/benchmark_native_ceres_thermo_regression.py`, `scripts/benchmark_native_regression.py`, and `tests/workflows/test_benchmark_native_regression.py` may be useful after the backend paths are real. Rewrite so skipped/unavailable backends cannot masquerade as passing production evidence.
- Documentation structure ideas in `docs/pages/parameter_regression.rst`, `docs/pages/diagnostics.rst`, and `docs/pages/package_architecture.rst` can inform future docs, but only after the implementation slice they describe is actually merged and validated.

## Files inspected

Inspection commands used from `codex/pr56-salvage-map`:

- `gh pr view 56 --repo tannerpolley/ePC-SAFT --json number,title,state,author,headRefName,baseRefName,body,files,commits,comments,reviews,url`
- `gh pr diff 56 --repo tannerpolley/ePC-SAFT --name-only`
- `git fetch origin pull/56/head:refs/remotes/origin/pr/56`
- `git diff --stat main...origin/pr/56`
- `git diff --name-status main...origin/pr/56`
- `git grep -n -i -E "finite[_ -]?difference|unsupported_derivative|native_hot_loop|backend_unavailable|residual-record|residual record|Ceres|CppAD|cost-function|cost function" origin/pr/56 -- CMakeLists.txt CMakePresets.json README.md docs/pages scripts src tests`

PR #56 changed 159 files with roughly 15.7k insertions and 647 deletions. The inspected changed-file groups were:

- Build and dependency plumbing: `CMakeLists.txt`, `CMakePresets.json`, `cmake/FindEigen3.cmake`, `scripts/build_epcsaft.py`, `scripts/doctor.py`.
- User and developer docs: `README.md`, `docs/pages/api_reference.rst`, `docs/pages/development_workflows.rst`, `docs/pages/diagnostics.rst`, `docs/pages/electrolyte_vle_reactive_workflow.rst`, `docs/pages/equilibrium_cookbook.rst`, `docs/pages/package_architecture.rst`, `docs/pages/parameter_regression.rst`, `docs/pages/user_options.rst`.
- Goal-note trail: `docs/goals/issue-53-native-regression-production/goal.md`, `docs/goals/issue-53-native-regression-production/state.yaml`, and notes `T001` through `T060`.
- Python runtime and public surfaces: `src/epcsaft/__init__.py`, `src/epcsaft/epcsaft.py`, `src/epcsaft/equilibrium.py`, `src/epcsaft/ipopt_backend.py`, `src/epcsaft/native_regression.py`, `src/epcsaft/parameters.py`, `src/epcsaft/reactive_regression.py`, `src/epcsaft/reactive_speciation.py`, `src/epcsaft/regression.py`, `src/epcsaft/runtime.py`.
- Benchmark modules and scripts: `src/epcsaft/benchmarks/*`, `scripts/benchmark_native_ceres_thermo_regression.py`, `scripts/benchmark_native_regression.py`, `scripts/benchmark_neutral_equilibrium.py`, `scripts/benchmark_reactive_regression.py`, `scripts/profile_regression_runtime.py`, `scripts/validate_hydrocarbon_regression.py`.
- Native bindings and kernels: `src/epcsaft/bindings.cpp`, `src/epcsaft/native/autodiff/*`, `src/epcsaft/native/contributions/*`, `src/epcsaft/native/epcsaft_*.cpp`, `src/epcsaft/native/epcsaft_*.h`, `src/epcsaft/native/regression/*`.
- Tests: `tests/api/*`, `tests/equilibrium/test_lle.py`, `tests/native/*`, `tests/regression/*`, `tests/workflows/*`, plus package plot smoke tests under `analyses/package_validation/package_plot_smokes/tests/plots/`.
- Data and maintenance scripts touched incidentally: `scripts/_env.py`, `scripts/_epcsaft_oop.py`, `scripts/data/*`, `scripts/plot_outputs.py`, `scripts/sync_equation_registry.py`, `run_pytest.py`, and one archival paper markdown file.

## Follow-up issues

- #58, Roadmap slice 1: eradicate finite-difference execution and terminology from production code. Use PR #56 only as a negative reference for why renaming `finite_difference` to `unsupported_derivative` is not enough.
- #59, Roadmap slice 2: add package-wide CppAD scalar and derivative substrate for the EOS harness. Rewrite any CppAD scalar scaffolding from current `main`; do not transplant PR #56 files.
- #60, Roadmap slice 3: templatize EOS contribution and property derivatives over scalar type. Use the PR only to identify affected kernels, then design the templating deliberately.
- #61, Roadmap slice 4: implement liquid-electrolyte SSM+DS Born derivatives for `d_born` and `f_solv`. Treat PR #56 Born derivative fragments as exploratory only.
- #62, Roadmap slice 5: wire CppAD/analytic derivatives into fugacity, activity, chemical-potential, pressure, density, and relative-permittivity APIs. Rebuild API routing with real derivative evidence.
- #63, Roadmap slice 6: rework equilibrium derivative and solver interfaces around CppAD, analytic derivatives, and implicit sensitivities. Do not salvage `unsupported_derivative` solver execution.
- #64, Roadmap slice 7: define staged reactive workflow boundaries using literature reaction constants first. Keep residual-record compatibility separate from production regression claims.
- #65, Roadmap slice 8: design easy pure-component and binary ePC-SAFT regression APIs. Use native regression schema sketches only after public API shape is settled.
- #66, Roadmap slice 9: implement Ceres regression backends for pure and binary ePC-SAFT parameter fitting. Rewrite Ceres cost functions and backend selection from scratch with tests that require the real backend.
- #67, Roadmap slice 10: add literature-backed regression validation using existing papers and analysis outputs. Convert any PR #56 benchmark ideas into validation only when backed by source data and pass/fail thresholds.
- #68, Roadmap slice 11: add Ceres/CppAD build matrix, CI gates, and honest capabilities reporting. This is the best home for CMake option, doctor, runtime reporting, and capability-claim cleanup.
- #69, Roadmap slice 12: decide and implement package organization for EOS, equilibrium, and regression modules. Use PR #56's sprawl as evidence that organization should be decided before another broad native-regression tranche.
