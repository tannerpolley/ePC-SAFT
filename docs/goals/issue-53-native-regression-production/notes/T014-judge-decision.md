# T014 Judge Decision

Result: `approved_with_constraints`

## Decision

Continue on `codex/issue-53-native-regression-production`. Do not fork a new `codex/issue-53-native-ceres-thermo-regression` branch. The existing branch already contains the Phase 0 groundwork and should carry the next tranche.

Do not open a completion PR or close issue #53 from the current residual-record-only state. If an early review PR is needed, it must be draft-only and explicitly framed as Phase 0 native residual-record groundwork.

## Approved Next Tranche

Implement `Phase 1 native thermodynamic hot loop` for the first production slice:

- native AD substrate under `src/epcsaft/native/autodiff/**`
- native thermodynamic row evaluator for `ReactiveSpeciation` and `ReactiveElectrolyteBubble`
- implicit sensitivities for those converged nested solves
- Ceres-owned parameter iteration in C++
- Python limited to validation, marshalling, and serialization
- finite difference allowed only behind `EPCSAFT_ALLOW_FINITE_DIFFERENCE_DEBUG=1`

Scope control: do not attempt whole-package CppAD conversion in one pass. Templatize only the call graph required by the supported production slice. Unsupported derivative paths must fail closed as `backend_unavailable` or an equally honest status.

## Ordered Tasks

1. `T015`: audit finite-difference, numerical Jacobian, Python objective, and Python residual-packing paths.
2. `T016`: add slice-scoped native AD scaffolding and focused derivative tests.
3. `T017`: implement native regression problem, parameter map, row workspace, and thermodynamic evaluator for `ReactiveSpeciation` and `ReactiveElectrolyteBubble`.
4. `T018`: add implicit sensitivities for converged nested speciation and fixed-liquid bubble solves.
5. `T019`: wrap the native evaluator in a Ceres `DynamicCostFunction` production solve loop.
6. `T020`: enforce debug-only finite-difference policy across public production derivative entry points.
7. `T021`: update capabilities, docs, and benchmarks proving `optimizer_backend=ceres`, `native_hot_loop=true`, `python_objective_used=false`, and non-finite-difference derivatives for the supported slice.

## Stop Conditions

- Stop if a worker keeps Python packing or evaluating production thermodynamic rows.
- Stop if a production path uses finite differences, SciPy/Python optimizer loops, or Ceres numeric differentiation as fallback.
- Stop if implementation expands into a broad whole-EOS templating rewrite beyond the first production-slice call graph.
- Stop if evidence still only proves residual-record solving rather than production thermodynamic regression.
