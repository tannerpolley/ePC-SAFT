# General reaction and equilibrium-constant convention layer

## Purpose

Define reaction constants, bases, standard states, and optional fitting/correction terms without guessing conventions.

This issue is part of the general reactive/electrolyte equilibrium readiness roadmap. It must keep `epcsaft` general-purpose and must not introduce application-specific public APIs.

## Branch

```text
codex/reaction-constant-conventions
```

## Phase

```text
Equilibrium
```

## Dependencies

```text
A
```

If any dependency is not merged into `origin/main`, do not start implementation.

## Scope

- ideal mole fraction constants
- molality constants
- activity-based thermodynamic constants
- apparent constants
- regularized correction terms
- fitted K coefficients
- warnings/errors for incompatible conventions

## Out of scope

- do not make MEA-specific APIs
- do not force reaction-constant fitting as the default

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

- `uv run python run_pytest.py tests/api/test_reactive_staged_workflow_contract.py tests/api/test_reaction_constant_conventions.py -q`
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
