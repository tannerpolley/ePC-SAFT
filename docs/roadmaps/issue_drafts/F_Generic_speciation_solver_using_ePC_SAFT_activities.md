# Generic speciation solver using ePC-SAFT activities

## Purpose

Harden generic liquid reactive/speciation solving with ideal/apparent and activity-based modes.

This issue is part of the general reactive/electrolyte equilibrium readiness roadmap. It must keep `epcsaft` general-purpose and must not introduce application-specific public APIs.

## Branch

```text
codex/generic-activity-speciation
```

## Phase

```text
Equilibrium
```

## Dependencies

```text
C, D
```

If any dependency is not merged into `origin/main`, do not start implementation.

## Scope

- ideal/apparent mode
- activity-based mode
- fixed K mode
- fitted K mode support where schema exists
- diagnostics for residual blocks and derivative status
- ePC-SAFT activities/fugacities where requested

## Out of scope

- do not implement MEA-specific speciation APIs
- do not fit reaction constants by default

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

- `uv run python run_pytest.py tests/api/test_reactive_staged_workflow_contract.py tests/native/test_reactive_speciation_implicit_sensitivity.py -q`
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
