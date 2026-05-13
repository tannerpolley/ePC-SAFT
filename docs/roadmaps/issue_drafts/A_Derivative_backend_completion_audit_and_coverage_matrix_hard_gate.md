# Derivative backend completion audit and coverage matrix hard gate

## Purpose

Audit and harden derivative coverage and capability reporting before any new equilibrium/regression implementation starts.

This issue is part of the general reactive/electrolyte equilibrium readiness roadmap. It must keep `epcsaft` general-purpose and must not introduce application-specific public APIs.

## Branch

```text
codex/backend-coverage-hard-gate
```

## Phase

```text
Backend
```

## Dependencies

```text
none
```

If any dependency is not merged into `origin/main`, do not start implementation.

## Scope

- expand derivative_coverage_matrix and runtime capability row-family matrix
- classify derivative paths as production_supported, blocker, or out_of_scope
- ensure capabilities match coverage
- add tests for coverage semantics
- produce follow-up blocker list if gaps remain

## Out of scope

- do not implement large solver features
- do not add application-specific public APIs

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

- `uv run python run_pytest.py tests/native/test_derivative_coverage_matrix.py tests/native/test_property_derivative_backend_contract.py tests/native/test_association_implicit_derivative_contract.py -q`
- `uv run python run_pytest.py tests/api/test_runtime_capabilities_dependency_gates.py -q`
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
