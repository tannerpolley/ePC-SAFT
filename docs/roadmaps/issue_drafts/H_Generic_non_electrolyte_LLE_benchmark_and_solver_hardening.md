# Generic non-electrolyte LLE benchmark and solver hardening

## Purpose

Prove ordinary two-liquid-phase splitting before layering electrolyte constraints.

This issue is part of the general reactive/electrolyte equilibrium readiness roadmap. It must keep `epcsaft` general-purpose and must not introduce application-specific public APIs.

## Branch

```text
codex/generic-non-electrolyte-lle
```

## Phase

```text
LLE
```

## Dependencies

```text
C
```

If any dependency is not merged into `origin/main`, do not start implementation.

## Scope

- phase split
- fugacity equality
- stability checks or anti-trivial-solution strategy
- clear phase diagnostics
- simple literature or repo benchmark

## Out of scope

- do not include electrolyte accounting in this issue
- do not force a poor benchmark to pass

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

- `uv run python run_pytest.py tests/equilibrium/test_lle.py tests/native/test_cppad_lle_derivatives.py -q`
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
