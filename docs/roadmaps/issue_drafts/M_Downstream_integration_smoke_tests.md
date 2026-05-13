# Downstream integration smoke tests

## Purpose

Prove downstream projects can use generic package APIs without private workaround code.

This issue is part of the general reactive/electrolyte equilibrium readiness roadmap. It must keep `epcsaft` general-purpose and must not introduce application-specific public APIs.

## Branch

```text
codex/downstream-integration-smokes
```

## Phase

```text
Downstream
```

## Dependencies

```text
F, G, I, J, K
```

If any dependency is not merged into `origin/main`, do not start implementation.

## Scope

- MEA-Thermodynamics smoke
- Lithium_Extraction smoke
- MEA-Absorption-Column smoke
- generic problem construction
- generic outputs consumed downstream
- no copied EOS implementation

## Out of scope

- do not add downstream-application-specific public APIs to epcsaft
- do not compute downstream metrics inside epcsaft package APIs

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

- `uv run python scripts/validate_project.py quick`
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
