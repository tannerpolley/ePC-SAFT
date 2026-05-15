# Staged And Optional Module Inventory

Phase 5 source scope from issue #120:

- `src/epcsaft/reactive_staged.py`
- `src/epcsaft/electrolyte_bubble.py`
- `src/epcsaft/reactive_electrolyte.py`
- `src/epcsaft/ipopt_backend.py`

## Classification

| Original path | Classification | Phase 5 decision |
| --- | --- | --- |
| `src/epcsaft/electrolyte_bubble.py` | Stable native public contract for fixed-liquid electrolyte bubble pressure. | Kept in place; already exposed through `epcsaft.electrolyte` and the package root. |
| `src/epcsaft/reactive_electrolyte.py` | Public sequential workflow: native chemical speciation followed by native fixed-liquid electrolyte bubble pressure. | Kept in place to preserve public imports and monkeypatchable test seams; module wording now says sequential. |
| `src/epcsaft/reactive_staged.py` | Public staged workflow helper for chemical equilibrium followed by an explicit phase route. | Kept in place to preserve public imports; added to the organized `epcsaft.reactive` namespace and tightened module wording. |
| `src/epcsaft/ipopt_backend.py` | Optional dependency bridge for explicit cyipopt-backed routes. | Implementation moved to `src/epcsaft/_optional_backends/ipopt.py`; top-level path is a compatibility module alias. |

## Public Surface

`ElectrolyteBubbleOptions`, `ElectrolyteBubbleResult`,
`ReactiveElectrolyteBubbleOptions`, `ReactiveElectrolyteBubbleResult`,
`solve_reactive_electrolyte_bubble`, `solve_reactive_electrolyte_bubble_sweep`,
`ReactiveStagedEquilibriumResult`, and `solve_reactive_staged_equilibrium`
remain public because they are imported by `epcsaft.__init__`, used by current
API tests, and represented in capability metadata or problem-object dispatch.

`epcsaft.ipopt_backend` remains importable only as a compatibility alias for
existing tests and explicit user imports. New code should treat IPOPT support as
an optional backend bridge owned by `epcsaft._optional_backends.ipopt`.

## Acceptance Notes

No solver algorithms, equations, reactive speciation logic, bubble-pressure
logic, or regression targets were changed. Phase 5 only clarified ownership and
relocated the optional backend bridge implementation.
