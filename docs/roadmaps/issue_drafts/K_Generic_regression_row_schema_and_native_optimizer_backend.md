# Generic regression row schema and native optimizer backend

## Purpose

Make regression generic around target rows, parameter maps, and native optimizer loops with real Jacobians.

This issue is part of the general reactive/electrolyte equilibrium readiness roadmap. It must keep `epcsaft` general-purpose and must not introduce application-specific public APIs.

## Branch

```text
codex/generic-regression-backend
```

## Phase

```text
Regression
```

## Dependencies

```text
B, C, E
```

If any dependency is not merged into `origin/main`, do not start implementation.

## Scope

- generic target-row compilation
- generic parameter maps
- Ceres preferred when it owns native loop
- other native optimizers allowed only with analytic/CppAD/implicit derivatives
- row diagnostics
- source summaries
- active bounds
- parameter movement

## Out of scope

- do not add fit_lithium_extraction_parameters or fit_mea_absorption APIs
- do not use finite difference or Python-owned production objective loops

## Required policy

```text
No finite difference.
CppAD for explicit algebraic derivatives.
Analytic formulas where exact and validated.
Implicit sensitivities for solved states.
No backend_unavailable for required workflows.
backend_unavailable only for explicitly out-of-scope workflows.
No application-specific public APIs.
```

## Required pre-read

```text
docs/roadmaps/general_reactive_electrolyte_equilibrium_readiness.md
docs/roadmaps/agent_dependency_plan.md
docs/roadmaps/agent_prompts/index.yaml
```

## Validation

- `uv run python scripts/build_epcsaft.py --enable-ceres --enable-cppad`
- `uv run python run_pytest.py tests/native/test_ceres_pure_regression.py tests/native/test_ceres_binary_regression.py tests/api/test_regression_api.py -q`
- `uv run python scripts/validate_project.py quick`

## Acceptance criteria

```text
focused PR
dependency gate documented
no finite difference
generic public API only
tests added or updated
capabilities/coverage updated where relevant
diagnostics honest
limitations documented
```
