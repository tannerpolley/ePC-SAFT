# Issue #117 Stage 0 Intake

Date: 2026-05-15

Issue: https://github.com/tannerpolley/ePC-SAFT/issues/117

## Dependency Status

Issue #117 remains dependency-gated behind issue #116 production electrolyte LLE completion. This intake records the current route map and fixture targets only; it does not open reactive source edits.

## Current Reactive Route Map

- `ReactivePhaseEquilibriumProblem.solve` in `src/epcsaft/equilibrium.py` delegates to `mixture.reactive_staged_equilibrium(...)`.
- `solve_reactive_staged_equilibrium` in `src/epcsaft/reactive_staged.py` solves reactive speciation first and then sends the resulting composition to a phase route.
- `_solve_phase_route` can route to LLE, VLE, electrolyte LLE, or flash-style phase routes, but the reaction and phase residuals are not solved in one coupled state.
- `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp` owns the native homogeneous chemical-equilibrium solve from issue #115. It is a subcomponent candidate, not a final coupled reactive phase-equilibrium solver.
- For electrolyte reactive LLE, issue #117 must wait for issue #116's production electrolyte LLE route rather than building on the current predictive electrolyte path.

## Existing Staged Diagnostics That Must Not Be Production Completion

- `workflow = chemical_equilibrium_then_phase_equilibrium`
- `reactive_workflow_class = staged`
- `reactive_phase_method = chemical_equilibrium_then_phase_equilibrium`
- `coupling_level = staged_not_full_simultaneous_nlp`
- `full_simultaneous_reactive_nlp = False`

These labels are valid for compatibility or diagnostic history only. They cannot be the accepted production route for issue #117.

## Required #117 Production Direction

Issue #117 requires one native coupled residual solve whose variables and residuals include:

- Phase fractions or phase amounts.
- Phase compositions or species amounts.
- Reaction extents.
- Ion-combination variables where charged species exist.
- Density, pressure, or closure variables where required by the selected model.
- Material or element balance residuals.
- Reaction equilibrium residuals from ePC-SAFT activities.
- Neutral interphase equilibrium residuals.
- Ionic interphase equilibrium residuals for electrolyte cases.
- Charge balance and normalization residuals or enforced transform diagnostics.
- A real Jacobian path with solved-state sensitivity diagnostics.
- Ceres solver diagnostics from the accepted coupled state.

## Candidate Fixture Targets

- Neutral reactive LLE: create a repo-contained esterification-style benchmark or smoke fixture following the issue #117 Ascani-style requirement. The current methanol/cyclohexane staged tests are useful route guards but are not enough because they monkeypatch or accept staged behavior.
- Reactive electrolyte LLE: select or create a repo-contained generic ion-exchange or reactive-transfer electrolyte LLE fixture after issue #116 is complete. Existing candidate data surfaces include Ascani-style water/butanol/salt fixtures and Khudaida 2026 salting-out data, but a source-backed reaction fixture still needs selection before implementation.

## Tests That Currently Would Not Prove Completion

- `tests/equilibrium/reactive/test_reactive_lle.py` currently proves staged route behavior.
- `tests/api/reactive/test_reactive_staged_equilibrium.py` currently proves the chemical-then-phase workflow.
- `tests/api/reactive/test_reactive_staged_workflow_contract.py` currently asserts staged diagnostics and unavailable coupled sensitivity status.

These tests should become compatibility or route-guard coverage. New production tests must fail if the accepted reactive LLE route remains staged.

## Source-Edit Entry Condition

Reactive source edits may start only after issue #116 is audited complete or after the board explicitly records a same-branch continuation where issue #116 was completed first.
