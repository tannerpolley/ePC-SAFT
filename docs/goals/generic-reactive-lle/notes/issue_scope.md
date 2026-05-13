# Issue Scope

GitHub issue: https://github.com/tannerpolley/ePC-SAFT/issues/93

Issue state: `OPEN`

Issue title: Generic reactive LLE and chemical phase equilibrium

## Current Scope

Combine reaction/speciation and LLE in a generic chemical-phase-equilibrium framework.

In-scope items from the current issue body:

- reaction coordinates
- element balances
- phase split
- fugacity equality
- reaction equilibrium
- nonnegativity
- Ascani 2023 benchmark attempt

## Current Constraints

- Keep `epcsaft` general-purpose.
- Do not introduce application-specific public APIs.
- Do not create extraction-specific or MEA-specific public APIs.
- Do not start with a full constrained NLP unless the issue explicitly decides that path.
- No finite difference.
- Use CppAD for explicit algebraic derivatives.
- Use analytic formulas where exact and validated.
- Use implicit sensitivities for solved states.
- Do not use `backend_unavailable` for required workflows.
- Use `backend_unavailable` only for explicitly out-of-scope workflows.

## Dependency Gate From Issue

The issue names dependencies `I`, `D`, and `F` and says not to start implementation if any dependency is not merged into `origin/main`.

Current dependency observations from this `/goal` run:

- I / issue #92, "Generic electrolyte LLE with distributed ions": `CLOSED`; PR #110 is `MERGED` with merge commit `4ff975c1f4760c27f235fb92d3f6d153cba3925f`.
- D / issue #87, "General reaction and equilibrium-constant convention layer": `CLOSED`; PR #99 is `MERGED` with merge commit `a81f950069d01834746150fd673e1107a2d2dcf2`.
- F / issue #89, "Generic speciation solver using ePC-SAFT activities": `CLOSED`; PR #105 is `MERGED` with merge commit `7443a5fee8affbeda3df94b9c0268d0149bfa783`.

The branch `codex/generic-reactive-lle` was fast-forwarded cleanly to `origin/main` at `d7485537aa472febe51cbbd7ad50ca62ba2ab165`, so the dependency gate is `GATE_PASSED`.

## Required Validation After Gate Pass

- `uv run python run_pytest.py tests/equilibrium/test_reactive_lle.py tests/api/test_reactive_staged_workflow_contract.py -q`
- `uv run python scripts/validate_project.py quick`
