# T010 Issue 117 Boundary Decision

## Decision

Approved for issue #117 source edits, with the existing T011 through T014 sequence as the implementation boundary.

The current branch satisfies the dependency gate because T008 approved issue #116 as available for #117. Source edits may proceed only within the queued #117 Worker slices and must preserve the issue requirement that production reactive phase equilibrium is one native coupled residual solve.

## Basis

- GitHub issue #117 remains open and was refreshed during this task. Its current title is "Replace staged reactive phase equilibrium with a production coupled native reactive LLE solver".
- T009 shows `ReactivePhaseEquilibriumProblem.solve(...)`, `kind="reactive_lle"`, and `kind="reactive_electrolyte_lle"` still route through staged chemical-equilibrium-then-phase-equilibrium behavior.
- T009 shows the native chemical-equilibrium solver exists as a subcomponent but is currently used as a pre-phase speciation stage.
- T008 shows the issue #116 electrolyte LLE Ceres route is available as a phase subcomponent with production residual/Jacobian diagnostics.

## Approved Worker Boundaries

### T011: Coupled Residual Foundation

T011 may implement the native problem representation, transformed variables, and residual-vector surface for neutral reactive LLE and reactive electrolyte LLE.

Required proof:

- One unknown vector represents reaction and phase state together.
- Residual blocks include material or element balance, reaction equilibrium, neutral phase equilibrium, ionic phase equilibrium where charged species exist, charge balance or an explicit transform guarantee, and normalization or an explicit transform guarantee.
- Result diagnostics cannot claim accepted production Ceres completion until T012.

Stop if reaction and phase residuals cannot be evaluated from one physical state.

### T012: Jacobian And Ceres Production Solve

T012 may add the coupled residual Jacobian, solved-state sensitivity diagnostics, and Ceres trust-region accepted solve.

Required proof:

- Accepted diagnostics report `solver_backend = ceres`.
- Accepted diagnostics report `solver_method = ceres_trust_region_coupled_reactive_phase_equilibrium`.
- Accepted diagnostics report `jacobian_available = true` and a real analytic, CppAD, or implicit Jacobian provenance.
- Staged/speciation-only and phase-only solves may initialize or compare, but cannot produce the accepted result.

Stop if the coupled residual Jacobian cannot be produced for the required variables.

### T013: Public API And Benchmarks

T013 may route the generic Python production APIs to the native coupled solve and add the required benchmark tests.

Required proof:

- `ReactivePhaseEquilibriumProblem.solve(...)` no longer calls the staged workflow for production reactive LLE or reactive electrolyte LLE.
- Any staged route remains explicitly named as staged/compatibility behavior.
- Neutral reactive LLE and reactive electrolyte LLE tests assert reaction residuals, phase residuals, material/element balance, charge balance where charged species exist, phase distance, solver backend, solver method, Jacobian backend, and derivative backend.
- Public APIs remain generic and do not expose downstream application metrics.

Stop if repo-contained benchmark data cannot be represented without fabrication.

### T014: Validation And Guard Repair

T014 must run the issue #117 validation ladder and route guards, then repair only narrow failures within #117 scope.

Required proof:

- Clean Ceres+CppAD build passes.
- Native equilibrium, reactive equilibrium, API reactive, quick validation, docs validation, and diff checks pass.
- Remaining route-guard matches are compatibility-only, initialization-only, or tests proving staged behavior is not production.

Stop if validation failure indicates the production coupled solver or production Jacobian is still absent.

## Rejected Completion Modes

Do not approve completion based on route inventory, schema changes, diagnostic labels, staged-only behavior, mocked payloads, or benchmark fixtures that do not exercise the native coupled solver. Do not approve production completion if the accepted result is assembled from separate chemical and phase solved states.
