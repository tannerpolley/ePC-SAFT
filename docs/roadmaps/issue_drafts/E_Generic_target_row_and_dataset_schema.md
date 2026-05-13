# Generic target-row and dataset schema

## Purpose

Create application-neutral regression and validation target-row schemas for future generic regression/equilibrium workflows.

This issue is part of the general reactive/electrolyte equilibrium readiness roadmap. It must keep `epcsaft` general-purpose and must not introduce application-specific public APIs.

## Branch

```text
codex/generic-target-row-schema
```

## Phase

```text
Regression
```

## Dependencies

```text
A
```

If any dependency is not merged into `origin/main`, do not start implementation.

## Scope

- pure density
- pure vapor pressure
- P-rho-T
- binary VLE
- binary LLE
- osmotic coefficient
- MIAC
- relative permittivity
- activity/fugacity
- speciation
- VLE partial pressure
- LLE phase composition
- regularization

## Out of scope

- do not add lithium-specific or MEA-specific public APIs
- do not implement optimizer internals here

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

- `uv run python run_pytest.py tests/api/test_regression_problem_schema.py tests/api/test_regression_api.py -q`
- `uv run python scripts/validate_project.py docs`

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
