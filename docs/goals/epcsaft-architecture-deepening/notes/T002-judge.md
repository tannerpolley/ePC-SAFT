# T002 Judge Receipt

Result: done

Decision: approved

Full outcome complete: false

## Rationale

Start with **Production Solver Path** because it is the dependency gate for route acceptance, rejection, and capability honesty. This is executable behavior, not a pass-through or documentation slice.

Keep **Capability Contract** last so claims follow proven route, regression, and parameter behavior.

ADR needed: no. The reordered start is reversible and grounded in Scout evidence, not a surprising public contract change.

## Ordered Candidates

1. Production Solver Path
2. Equilibrium Problem
3. Target Dataset and Regression Problem
4. Parameter Family
5. Capability Contract

## First Worker Package

Objective: deepen the Production Solver Path result gate so accepted/rejected native route translation has one clear interface, focused tests, and no duplicated diagnostics/capability glue that can overclaim unsupported routes.

Allowed files:

- `src/epcsaft/equilibrium.py`
- `src/epcsaft/equilibrium_core/native_results.py`
- `src/epcsaft/equilibrium_core/native_requests.py`
- `src/epcsaft/equilibrium_core/`
- `src/epcsaft/native/equilibrium_nlp/result_builder.cpp`
- `src/epcsaft/native/equilibrium_nlp/result_builder.h`
- `src/epcsaft/native/equilibrium_nlp/route_builders.cpp`
- `src/epcsaft/native/equilibrium_nlp/route_builders.h`
- `tests/equilibrium/`
- `tests/native/equilibrium/`
- `tests/api/equilibrium/`
- `tests/workflows/paper_validation/`

Verify:

- `uv run python run_pytest.py tests/native/equilibrium tests/equilibrium/electrolyte tests/api/equilibrium -q`
- `uv run python scripts/dev/check_text_gates.py`

Stop if:

- A clean `_core` rebuild is required while another task may import `_core`.
- Acceptance/rejection semantics would change without a focused regression test.
- The Worker would need regression, parameter schema, downstream-specific, or broad capability-contract files.
- The slice devolves into pass-through modules, inventory, diagnostics-only proof, staged-workflow-only proof, or synthetic-only route assertions.
