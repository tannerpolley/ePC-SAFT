# Module Ownership And Split Plan

This plan is the Phase 6 deliverable for issue #120. It documents ownership and
future split boundaries for the largest source modules without changing runtime
behavior.

## Non-Goals

- Do not change EOS equations, association internals, equilibrium residuals,
  regression objectives, parameter values, benchmark targets, or public
  scientific outputs in a split-only cleanup.
- Do not change pybind signatures or native request/result payload shapes without
  tests that exercise the exact public route.
- Do not close implementation issues with this document. A split plan is cleanup
  evidence, not proof of a production workflow.

## Size Snapshot

| Path | Current line count | Primary ownership |
| --- | ---: | --- |
| `src/epcsaft/regression.py` | 3356 | Python regression API, schema, input normalization, native-call orchestration, result formatting, CSV updates. |
| `src/epcsaft/native/epcsaft_equilibrium.cpp` | 2900 | Native equilibrium kernels, stability, LLE, electrolyte LLE, bubble/dew, residual evaluation. |
| `src/epcsaft/epcsaft.py` | 2754 | Public mixture/state facade, state/property wrappers, equilibrium dispatch, native regression shims. |
| `src/epcsaft/native/epcsaft_regression.cpp` | 2630 | Native regression residuals, native optimizer loops, Ceres cost functions, derivative-backed residual packing. |
| `src/epcsaft/equilibrium.py` | 1870 | Python equilibrium dataclasses, input normalization, public wrappers, native payload conversion. |
| `src/epcsaft/reactive_regression.py` | 1683 | Reactive electrolyte regression row/batch schema, objective evaluation, reporting helpers. |
| `src/epcsaft/parameters.py` | 1486 | Dataset resolution, parameter loading, schema migration, electrolyte parameter normalization. |
| `src/epcsaft/bindings.cpp` | 1517 | pybind11 boundary, native request conversion, result conversion, exposed native entrypoints. |
| `src/epcsaft/reactive_speciation.py` | 1435 | Reactive-speciation public API, reaction conventions, native and activity-fixed-point routes, sweep continuation. |

## Ownership Matrix

| Path | Runtime responsibilities | Validation/docs responsibilities | Safe future splits | Unsafe future splits | Tests that must move with a split |
| --- | --- | --- | --- | --- | --- |
| `src/epcsaft/regression.py` | Owns public regression problem/result objects, target-row parsing, Python-to-native payload preparation, native result normalization, and public fit helpers. | Owns regression API examples and package-facing regression contracts. | Move dataclasses and target-row validators to `regression_models.py`; move CSV/file writers to `regression_io.py`; move native payload builders to `regression_native_payloads.py`; keep compatibility imports in `regression.py`. | Moving target-family semantics or native solver dispatch without `tests/api/regression` and `tests/native/ceres`; changing residual scales, bounds, or parameter names during a split. | `tests/api/regression`, `tests/native/ceres`, `scripts/validation/validate_hydrocarbon_regression.py`, `scripts/benchmarks/profile_regression_runtime.py`. |
| `src/epcsaft/native/epcsaft_equilibrium.cpp` | Owns compiled phase-equilibrium algorithms, stability trials, electrolyte LLE residuals, bubble/dew calculations, fixed-liquid electrolyte bubble pressure, and native residual evaluators. | Owns native behavior that supports equilibrium docs and diagnostics. | Move neutral bubble/dew routines, stability helpers, and electrolyte bubble helpers into separately compiled native files with unchanged declarations in `epcsaft_equilibrium.h`. | Splitting transformed-basis electrolyte LLE, phase-label logic, density closure, or acceptance gates without native contract tests and API tests. | `tests/equilibrium`, `tests/native/contracts`, `tests/api/runtime`, `tests/api/package`. |
| `src/epcsaft/epcsaft.py` | Owns `ePCSAFTMixture`, `ePCSAFTState`, public state/property methods, high-level equilibrium dispatch, and legacy-compatible native regression shims. | Owns user-facing facade behavior and method-level examples. | Move pure state/property helper functions to `state_api.py`; move public equilibrium dispatch helpers to `mixture_equilibrium.py`; keep public methods and imports stable on `ePCSAFTMixture`. | Moving public methods or changing dispatch aliases without runtime and downstream smoke coverage. | `tests/api/runtime`, `tests/equilibrium`, `tests/api/package`, docs examples in `docs/pages/api_reference.rst` and cookbook pages. |
| `src/epcsaft/native/epcsaft_regression.cpp` | Owns native regression objective evaluation, residual/Jacobian packing, bounded transforms, candidate starts, Ceres cost functions, and native fit entrypoints. | Owns native proof behind regression capability metadata. | Move pure-neutral helpers, generic target evaluators, and Ceres cost-function classes into compiled sibling files with declarations kept in `epcsaft_electrolyte.h` until a dedicated regression header exists. | Moving objective math, transform definitions, residual ordering, or target-family indexing without exact native/API regression tests. | `tests/native/ceres`, `tests/api/regression`, `tests/native/contracts/test_ceres_cppad_build_contract.py`. |
| `src/epcsaft/equilibrium.py` | Owns Python equilibrium option/result dataclasses, feed normalization, salt/formula transforms, native result conversion, and public neutral/electrolyte wrapper functions. | Owns user-facing equilibrium API contracts and JSON-like result shapes. | Move dataclasses to `equilibrium_models.py`; move electrolyte feed/salt normalization to `equilibrium_core/electrolyte_inputs.py`; move native payload conversion to `equilibrium_native_payloads.py`. | Splitting result conversion or feed normalization without tests for charge neutrality, formula/explicit transforms, and public route aliases. | `tests/equilibrium`, `tests/api/runtime`, `tests/api/package`, docs cookbook examples. |
| `src/epcsaft/reactive_regression.py` | Owns reactive-electrolyte row/batch models, residual packing, objective/Jacobian wrappers, fit result shaping, and report writers. | Owns reactive regression docs and benchmark/report helper surfaces. | Move dataclasses to `reactive_regression_models.py`; move report writers to `reactive_regression_io.py`; move row residual packing to `reactive_regression_residuals.py`. | Moving continuation seed handling, residual-family names, or parameter application semantics without reactive API/regression tests. | `tests/api/reactive`, `tests/api/regression`, `scripts/benchmarks/helpers/reactive_regression.py`. |
| `src/epcsaft/parameters.py` | Owns dataset discovery, parameter bundle loading, electrolyte schema normalization, component matching, and molality/mole-fraction helpers. | Owns parameter docs and dataset bundle compatibility contracts. | Move dataset-path discovery to `parameter_sources.py`; move schema migration/default pruning to `parameter_migration.py`; move electrolyte normalization to `electrolyte_parameters.py`. | Moving default handling, pure/binary matrix conventions, or charge/electrostatic schema normalization without dataset validation tests. | `tests/api/package`, `tests/workflows/repo`, dataset validation tests, downstream integration smokes. |
| `src/epcsaft/bindings.cpp` | Owns pybind11 converters, native request parsing, result serialization, exception mapping, and `_core` entrypoint exposure. | Owns native/Python ABI stability for tests and downstream users. | Move converter helpers into header-only or compiled binding helper files only after CMake explicitly lists them; group regression, equilibrium, and derivative binding converters by concern. | Changing exported function names, request keys, exception classes, or result keys during a mechanical split. | `tests/api/runtime`, `tests/api/regression`, `tests/equilibrium`, `tests/native/contracts`, `scripts/dev/build_dist.py`. |
| `src/epcsaft/reactive_speciation.py` | Owns reaction conventions, reactive-speciation options/results, native speciation calls, activity fixed-point route, sweep continuation, and structured failure results. | Owns reactive-speciation API docs and validation examples. | Move reaction convention/dataclasses to `reactive_models.py`; move sweep/continuation helpers to `reactive_sweeps.py`; move activity fixed-point implementation to `reactive_activity_routes.py`. | Moving reaction-standard-state semantics, mass/charge residual gates, or native payload construction without reactive speciation tests. | `tests/api/reactive`, `tests/equilibrium/reactive`, docs reactive workflow pages. |

## Safe Split Order

1. Extract pure Python dataclasses and serializers first, with compatibility
   imports left in the original modules.
2. Extract report writers and CSV/file helpers next because they have narrow
   side effects and clear tests.
3. Extract native binding converters only after CMake owns the new files
   explicitly and `build_dist.py` passes.
4. Extract native algorithms last, one compiled source file at a time, with the
   public declarations and pybind entrypoints unchanged.

## Required Gates For Future Split PRs

Every future split PR must include:

- A before/after import scan for public paths.
- Focused tests named in the ownership matrix.
- `uv run python scripts/dev/build_epcsaft.py` when native files or bindings
  move.
- `uv run python scripts/dev/build_dist.py` when package source ownership or
  CMake ownership changes.
- `uv run python scripts/dev/validate_project.py quick` before handoff.

## Explicitly Deferred

The following are not cleanup-only work and require separate implementation
issues:

- Adding new derivative coverage, residual families, or optimizer behavior.
- Reworking association, density, equilibrium, or regression math.
- Changing public route aliases or capability claims.
- Retiring compatibility imports for currently documented public modules.
