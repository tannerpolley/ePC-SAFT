# Literature benchmark suite

## Purpose

Inventory and then build generic literature benchmarks for package-level confidence.

This issue is part of the general reactive/electrolyte equilibrium readiness roadmap. It must keep `epcsaft` general-purpose and must not introduce application-specific public APIs.

## Branch

```text
codex/literature-benchmark-suite
```

## Phase

```text
Benchmarks
```

## Dependencies

```text
none
```

If any dependency is not merged into `origin/main`, do not start implementation.

## Scope

- fixture inventory can start early
- implementation waits on relevant solver/regression issues
- MEA simple workflow benchmark
- MDEA ePC-SAFT benchmark
- Figiel 2025 SSM+DS Born benchmark
- Held 2014 revised ePC-SAFT benchmark
- non-electrolyte LLE benchmark
- Ascani 2022 electrolyte LLE benchmark
- Ascani 2023 reactive LLE benchmark
- Khudaida salting-out LLE benchmark
- Hubach/Yu lithium-related equilibrium benchmark

## Out of scope

- do not implement benchmark tests that depend on missing APIs
- do not require downstream repo access for package-level benchmarks

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

- `uv run python scripts/validate_project.py docs`
- `uv run python run_pytest.py tests/regression/test_literature_pure_parameter_regression.py tests/regression/test_literature_binary_kij_regression.py tests/regression/test_figiel_2025_born_parameter_parity.py -q`

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
