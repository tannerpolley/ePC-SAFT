# Staged And Optional Module Inventory

Phase 5 source scope from issue #120:

- `src/epcsaft/reactive_staged.py`
- `src/epcsaft/electrolyte_bubble.py`
- `src/epcsaft/reactive_electrolyte.py`

## Classification

| Original path | Classification | Phase 5 decision |
| --- | --- | --- |
| `src/epcsaft/electrolyte_bubble.py` | Stable native public contract for fixed-liquid electrolyte bubble pressure. | Kept in place; already exposed through `epcsaft.electrolyte` and the package root. |
| `src/epcsaft/reactive_electrolyte.py` | Public sequential workflow: native chemical speciation followed by native fixed-liquid electrolyte bubble pressure. | Kept in place to preserve public imports and monkeypatchable test seams; module wording now says sequential. |
| `src/epcsaft/reactive_staged.py` | Public staged workflow helper for chemical equilibrium followed by an explicit phase route. | Kept in place to preserve public imports; added to the organized `epcsaft.reactive` namespace and tightened module wording. |

## Public Surface

`ElectrolyteBubbleOptions`, `ElectrolyteBubbleResult`,
`ReactiveElectrolyteBubbleOptions`, `ReactiveElectrolyteBubbleResult`,
`solve_reactive_electrolyte_bubble`, `solve_reactive_electrolyte_bubble_sweep`,
`ReactiveStagedEquilibriumResult`, and `solve_reactive_staged_equilibrium`
remain public because they are imported by `epcsaft.__init__`, used by current
API tests, and represented in capability metadata or problem-object dispatch.

Historical Python IPOPT bridge notes are obsolete. New IPOPT work belongs in the
native constrained-NLP adapter, not a Python optional backend module.

## Acceptance Notes

No solver algorithms, equations, reactive speciation logic, bubble-pressure
logic, or regression targets were changed. Phase 5 only clarified ownership and
clarified the optional backend bridge ownership at the time.
