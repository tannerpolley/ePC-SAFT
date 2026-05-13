# Explicit CppAD parameter derivatives for EOS/property APIs

## Purpose

Implement or verify CppAD/analytic derivatives for explicit algebraic parameter effects used by properties and regression.

This issue is part of the general reactive/electrolyte equilibrium readiness roadmap. It must keep `epcsaft` general-purpose and must not introduce application-specific public APIs.

## Branch

```text
codex/cppad-explicit-parameter-derivatives
```

## Phase

```text
Backend
```

## Dependencies

```text
A
```

If any dependency is not merged into `origin/main`, do not start implementation.

## Scope

- m, sigma, epsilon explicit effects
- k_ij, l_ij, k_hb_ij direct effects where algebraic
- d_born, f_solv, dielectric parameters
- pressure/fugacity/activity/chemical-potential parameter derivatives where explicit
- regression APIs request derivatives by parameter name

## Out of scope

- do not implement implicit solved-state sensitivities in this issue
- do not tape solver loops

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

- `uv run python scripts/build_epcsaft.py --enable-cppad`
- `uv run python run_pytest.py tests/native/test_cppad_pressure_derivatives.py tests/native/test_cppad_fugacity_derivatives.py tests/native/test_cppad_activity_derivatives.py tests/native/test_cppad_relative_permittivity_derivatives.py -q`
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
