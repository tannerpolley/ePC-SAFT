# T009 Issue 117 Route Map

## Scope

This scout pass maps the current issue #117 reactive phase-equilibrium routes after the issue #116 audit unlocked #117 work.

## Current Public Reactive Route

- `ReactivePhaseEquilibriumProblem.solve(...)` is still staged: `src/epcsaft/equilibrium.py:263` defines the problem and `src/epcsaft/equilibrium.py:271` delegates to `mixture.reactive_staged_equilibrium(...)`.
- `ePCSAFTMixture.equilibrium(...)` maps `reactive_lle`, `reactive_lle_flash`, `reactive_electrolyte_lle`, and `reactive_electrolyte_lle_flash` onto `self.reactive_staged_equilibrium(...)` in `src/epcsaft/epcsaft.py:646` through `src/epcsaft/epcsaft.py:673`.
- `reactive_lle` forces `phase_kind = "lle_flash"`; `reactive_electrolyte_lle` forces `phase_kind = "electrolyte_lle"`. That means the route can now call the #116 Ceres electrolyte LLE solver as a phase subroute, but the reactive solve itself remains sequential.

## Current Staged Workflow

- `solve_reactive_staged_equilibrium(...)` first calls `solve_reactive_speciation(...)`, then builds a new mixture from the solved chemical composition, then calls `_solve_phase_route(...)`.
- Its diagnostics explicitly report `workflow = chemical_equilibrium_then_phase_equilibrium`, `reactive_workflow_class = staged`, `coupling_level = staged_not_full_simultaneous_nlp`, and `full_simultaneous_reactive_nlp = False`.
- This is useful as a compatibility/reference workflow, but it does not satisfy issue #117 production behavior because reaction and phase residuals are not evaluated in one solved state.

## Native Chemical-Equilibrium Path

- `src/epcsaft/reactive_speciation.py:538` calls `_core._solve_chemical_equilibrium_native(...)`.
- Native entrypoints are present and tested: `_solve_chemical_equilibrium_native` and `_evaluate_chemical_equilibrium_residual_native` are bound and asserted in `tests/native/equilibrium/test_chemical_equilibrium_native_api.py:93` through `tests/native/equilibrium/test_chemical_equilibrium_native_api.py:119`.
- Native chemical-equilibrium diagnostics identify `_solve_chemical_equilibrium_native` in `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp:1051`, `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp:1485`, and `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp:1612`.
- This native path is a valid #117 subcomponent, but the current API uses it as a pre-phase speciation stage rather than as a coupled residual block inside a reactive phase-equilibrium solve.

## Issue 116 Solver Availability

- T008 approved #116 as available for #117.
- The native electrolyte LLE route now reports `solver_backend = ceres`, `solver_method = ceres_trust_region_residual_solve`, `jacobian_backend = cppad_implicit`, and `derivative_backend = cppad_implicit` in `src/epcsaft/native/epcsaft_equilibrium.cpp:3372` through `src/epcsaft/native/epcsaft_equilibrium.cpp:3381`.
- Public fallback diagnostics in `src/epcsaft/equilibrium.py:2030` through `src/epcsaft/equilibrium.py:2036` report the same production route when native execution is not available in the current runtime.
- `tests/api/equilibrium/test_electrolyte_lle_problem_production_route.py:9` through `tests/api/equilibrium/test_electrolyte_lle_problem_production_route.py:28` prove `ElectrolyteLLEProblem` routes through the generic Ceres electrolyte LLE production diagnostics.

## Tests That Still Pass On Staged Behavior

- `tests/equilibrium/reactive/test_reactive_lle.py:80` through `tests/equilibrium/reactive/test_reactive_lle.py:116` assert that `reactive_lle` stages chemical equilibrium into a neutral LLE split and reports `full_simultaneous_reactive_nlp = False`.
- `tests/equilibrium/reactive/test_reactive_lle.py:122` through `tests/equilibrium/reactive/test_reactive_lle.py:144` assert `ReactivePhaseEquilibriumProblem` routes to the generic staged LLE path.
- `tests/api/reactive/test_reactive_staged_equilibrium.py:27` through `tests/api/reactive/test_reactive_staged_equilibrium.py:53` assert the exported staged helper returns chemical and phase results with `workflow = chemical_equilibrium_then_phase_equilibrium`.
- `tests/api/reactive/test_reactive_staged_workflow_contract.py:44` through `tests/api/reactive/test_reactive_staged_workflow_contract.py:64` assert staged workflow policy and `full_simultaneous_reactive_nlp = False`.
- These tests should be retained only for explicit staged compatibility behavior or replaced/augmented with tests that reject staged behavior as the production route for `ReactivePhaseEquilibriumProblem`, `reactive_lle`, and `reactive_electrolyte_lle`.

## Candidate #117 Fixtures

- Neutral reactive LLE: reuse the existing methanol/cyclohexane-style synthetic fixture only as an API and residual-coupling proof. It is not a literature benchmark.
- Reactive electrolyte LLE: use the #116 distributed-ion electrolyte LLE fixture as the phase-equilibrium base, with a simple repo-contained reaction that shifts neutral/ionic composition and can be checked by element/material balance, reaction residuals, phase residuals, charge balance, and phase distance.
- Do not add downstream extraction, selectivity, absorber, or solvent-screening metrics to the package API.

## Boundary Recommendation For T010

Approve #117 source edits only if the next Worker package replaces production reactive phase-equilibrium behavior with one native coupled residual solve. The staged workflow can remain as an explicitly named compatibility/helper route, but it must not be accepted production behavior for `ReactivePhaseEquilibriumProblem.solve(...)`, `kind="reactive_lle"`, or `kind="reactive_electrolyte_lle"`.
