# Generic electrolyte LLE with distributed ions

## Purpose

Build electrolyte LLE on the generic LLE foundation using ion-based notation and distributed ions.

This issue is part of the general reactive/electrolyte equilibrium readiness roadmap. It must keep `epcsaft` general-purpose and must not introduce application-specific public APIs.

## Branch

```text
codex/generic-electrolyte-lle
```

## Phase

```text
LLE
```

## Dependencies

```text
H, C
```

If any dependency is not merged into `origin/main`, do not start implementation.

## Scope

- ion-based species
- phase electroneutrality
- distributed ions
- mixed solvents
- mixed electrolytes
- charge-balance diagnostics
- Ascani 2022 benchmark attempt

## Out of scope

- do not create lithium-extraction-specific public APIs
- if Ascani 2022 Case Study 2 is inconsistent, document it rather than forcing a fake pass

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

- `uv run python run_pytest.py tests/equilibrium/test_electrolyte_lle.py tests/native/test_cppad_lle_derivatives.py -q`
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
