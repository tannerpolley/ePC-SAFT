# T011 Judge Receipt

Result: done

Decision: rejected

Full outcome complete: false

## Rationale

The slice improves typed problem diagnostics, but it remains mostly a wrapper over mixture methods and shared kind strings. `mixture.equilibrium(...)` still owns substantial route dispatch and validation complexity, and the `T010` receipt did not create a concrete follow-up to remove that remaining dispatcher ownership.

## Evidence

- `T010` adds `equilibrium_route` and `route_reason` for typed problem solves.
- `src/epcsaft/equilibrium.py` classifies a kind string, calls the existing mixture method, then adds diagnostics.
- `src/epcsaft/epcsaft.py` still contains route classification, kind branching, validation, and diagnostic threading.
- Direct typed problem-object coverage is still TPFlash-heavy.

## Missing Evidence

- Focused tests for typed `BubblePoint`, `DewPoint`, `StabilityAnalysis`, `LLEProblem`, `ElectrolyteLLEProblem`, and `ElectrolyteBubblePoint` diagnostics/error paths.
- A follow-up implementation that removes duplicated route ownership from `mixture.equilibrium(...)`, or a board-approved contract that makes that remaining dispatcher complexity explicit.

## Required Follow-Up

Add `T012` before `T030`.
