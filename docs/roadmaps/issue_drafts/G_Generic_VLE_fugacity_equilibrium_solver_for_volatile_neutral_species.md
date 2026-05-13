# Generic VLE/fugacity-equilibrium solver for volatile neutral species

## Purpose

Implement or harden generic fugacity-equilibrium routes for volatile neutral species.

This issue is part of the general reactive/electrolyte equilibrium readiness roadmap. It must keep `epcsaft` general-purpose and must not introduce application-specific public APIs.

## Branch

```text
codex/generic-vle-fugacity-equilibrium
```

## Phase

```text
Equilibrium
```

## Dependencies

```text
C
```

If any dependency is not merged into `origin/main`, do not start implementation.

## Scope

- direct volatile partial pressure from liquid fugacity where valid
- bubble pressure
- bubble temperature
- dew pressure
- dew temperature
- TP flash
- diagnostics report route used

## Out of scope

- do not create MEA bubble-pressure-specific public APIs
- do not assume ions distribute to vapor unless problem explicitly models it

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

- `uv run python run_pytest.py tests/equilibrium/test_derivative_policy.py tests/native/test_cppad_bubble_derivatives.py -q`
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
