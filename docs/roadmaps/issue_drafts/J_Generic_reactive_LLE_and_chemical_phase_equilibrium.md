# Generic reactive LLE and chemical phase equilibrium

## Purpose

Combine reaction/speciation and LLE in a generic chemical-phase-equilibrium framework.

This issue is part of the general reactive/electrolyte equilibrium readiness roadmap. It must keep `epcsaft` general-purpose and must not introduce application-specific public APIs.

## Branch

```text
codex/generic-reactive-lle
```

## Phase

```text
Reactive
```

## Dependencies

```text
I, D, F
```

If any dependency is not merged into `origin/main`, do not start implementation.

## Scope

- reaction coordinates
- element balances
- phase split
- fugacity equality
- reaction equilibrium
- nonnegativity
- Ascani 2023 benchmark attempt

## Out of scope

- do not create extraction-specific or MEA-specific public APIs
- do not start with a full constrained NLP unless the issue explicitly decides that path

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

- `uv run python run_pytest.py tests/equilibrium/test_reactive_lle.py tests/api/test_reactive_staged_workflow_contract.py -q`
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
