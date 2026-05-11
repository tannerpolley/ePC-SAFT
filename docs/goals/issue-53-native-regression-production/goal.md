# Issue #53 Native Regression Production Goal

## Objective

Fully implement GitHub issue #53: replace package-owned production regression fitting with a native C++ stack centered on Ceres plus CppAD/analytic/implicit derivatives, remove Python/SciPy/NumPy finite-difference optimizer paths from production fitting, and provide native mixed pressure/speciation reactive-electrolyte regression that downstream projects can call directly.

Issue URL: https://github.com/tannerpolley/ePC-SAFT/issues/53

## Completion Proof

The goal is complete only when all of these are true:

- Production regression fitting is native C++.
- Ceres is the default package-owned bounded nonlinear least-squares backend when available.
- CppAD, analytic derivatives, or implicit sensitivity derivatives own production derivative paths.
- Python remains an API/data/serialization layer and does not optimize production regression.
- SciPy is not used for package-owned production fitting.
- Finite differences are rejected for production regression unless an explicit debug gate is enabled.
- Native mixed pressure/speciation reactive electrolyte regression runs with fixed-shape residuals, bounded parameters, status diagnostics, row diagnostics, and no `bounded_incomplete`.
- At least one SSM+DS/Born-related parameter and at least one `k_ij` parameter are represented in native regression tests or benchmark fixtures.
- Package-owned native benchmark cases include tiny neutral/binary/reactive cases and an MEA-style 35-row public surrogate.
- Docs explain native C++ regression architecture, Ceres, CppAD/autodiff/implicit derivative policy, finite-difference debug-only policy, status contract, benchmark commands, and downstream usage.
- Required tests, benchmarks, docs build, lint, and formatting pass.
- A PR is created with the issue-required evidence, checks pass, and the issue is closed or clearly linked to the merge.

## Non-Goals And Guardrails

- Do not rewrite the entire EOS from scratch.
- Do not hard-code MEA-specific chemistry into generic optimizer code.
- Do not remove Python APIs for inputs, examples, docs, or result serialization.
- Do not remove legacy parameter dictionaries.
- Do not make IPOPT the default before a true constrained NLP backend exists.
- Do not weaken electrolyte LLE, reactive speciation, or reactive bubble acceptance gates.
- Do not silently relabel incomplete work as production.
- Do not modify generated or analysis-only artifacts unless a task explicitly requires it.

## Required Implementation Phases

1. Map the current regression and native build state against issue #53 without editing files.
2. Decide the minimum viable native backend slice that can be completed without destabilizing the whole EOS.
3. Add native regression data contracts and pybind/Python conversion boundaries.
4. Add CMake dependency handling for Ceres and CppAD with robust fallback/capability reporting.
5. Implement native residual evaluation with fixed-shape outputs and penalty residuals for recoverable row failures.
6. Implement Ceres bounded least-squares solving, bounds, robust losses, statuses, diagnostics, and result serialization.
7. Add the production derivative policy: no finite differences in production; CppAD/analytic/implicit where implemented; explicit debug-only finite-difference gate.
8. Migrate Python public wrappers so production fit calls invoke native C++ once and no longer optimize in Python.
9. Add package-owned native benchmark cases and `scripts/benchmark_native_regression.py`.
10. Update docs, capabilities, tests, and downstream guidance.
11. Run validation, benchmark evidence, final audit, PR, merge, and branch cleanup.

## Required Validation Commands

Run at minimum before completion:

```powershell
uv sync --no-install-project
uv run python scripts/build_epcsaft.py
uv run python run_pytest.py tests/native/test_native_regression_types.py tests/native/test_native_ceres_regression.py tests/native/test_native_regression_autodiff.py tests/native/test_native_reactive_regression.py -q
uv run python run_pytest.py tests/api/test_reactive_regression.py tests/api/test_runtime.py tests/api/test_parameter_schema.py -q
uv run python scripts/benchmark_native_regression.py --warmup 1 --repeat 3
uv run python scripts/validate_project.py quick
uv run python scripts/validate_project.py docs
uv run ruff check src tests docs
uv run black --check src tests docs
```

If native IPOPT scaffolding is enabled beyond docs/capabilities:

```powershell
uv run python run_pytest.py tests/native/test_native_ipopt_regression.py -q
```

## PR Reporting Requirements

The final PR must include:

```markdown
## Summary

## Native Backend Architecture

## Optimizer And Derivative Policy

## Regression Status Contract

## Benchmark Evidence

## Validation Commands

## Known Limitations

## Downstream Impact
```

## Starter Command

```text
/goal Follow docs/goals/issue-53-native-regression-production/goal.md.
```
