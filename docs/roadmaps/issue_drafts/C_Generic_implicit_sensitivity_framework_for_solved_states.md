# Generic implicit sensitivity framework for solved states

## Purpose

Create reusable implicit-sensitivity machinery for solved internal states across EOS, equilibrium, and regression.

This issue is part of the general reactive/electrolyte equilibrium readiness roadmap. It must keep `epcsaft` general-purpose and must not introduce application-specific public APIs.

## Branch

```text
codex/generic-implicit-sensitivity-framework
```

## Phase

```text
Backend
```

## Dependencies

```text
A, B
```

If any dependency is not merged into `origin/main`, do not start implementation.

## Scope

- association site fractions
- density root
- speciation solve
- VLE root
- LLE phase split
- reactive LLE solve
- ImplicitSolveResult with state/residual/jacobians/sensitivity/backend/status/diagnostics

## Out of scope

- do not tape iterative solver loops as production derivatives
- do not merge broad equilibrium rewrites into this issue

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

- `uv run python run_pytest.py tests/native/test_association_implicit_derivative_contract.py tests/native/test_reactive_speciation_implicit_sensitivity.py -q`
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
