# Non-Negotiable Solver Correction

This issue update supersedes the earlier broad wording for the equilibrium implementation work.

The current package already has native-looking equilibrium paths. That is not enough. This issue is complete only when the accepted production solve is performed by a mature native solver backend with production Jacobians and strict thermodynamic residual checks.

Do not close this issue with any of the following:

```text
inventory
manifest
schema-only support
diagnostic-only route
staged-only route
synthetic-only fixture
mocked result payload
documented limitation
capability label
hand-coded simplex minimizer
hand-coded Powell/Nelder-Mead production solver
derivative-approximation Jacobian
```

The old native residual minimizers may be used only as deleted code, test-only historical reference, or explicitly named seed generators if the production result does not depend on them. They must not be the accepted production solver.

The accepted production solver must use one of:

```text
Ceres trust-region residual solve with explicit analytic / CppAD / implicit Jacobian
IPOPT constrained thermodynamic NLP with analytic / CppAD derivative data
another mature external nonlinear solver with the same derivative and acceptance requirements
```

For this issue, the expected implementation is **Ceres trust-region residual solving on transformed variables**, not direct Gibbs-energy minimization. Constraints must be enforced by variable transformations and residual blocks. Gibbs / g-hat and TPD calculations are allowed as stability, seed-selection, and acceptance checks. They are not the production solve by themselves.

Allowed derivative provenance:

```text
analytic
cppad
analytic_implicit
cppad_implicit
ceres_solver_with_cppad_jacobian
```

Forbidden derivative provenance for production completion:

```text
derivative approximation
not_required for an accepted nonlinear equilibrium solve
missing required solve-path implementation
```

Package APIs must remain generic. Do not add public APIs named after MEA, lithium extraction, extraction efficiency, selectivity, distribution coefficient, absorber columns, or solvent screening.


# Issue #117 Replacement Body — Native Coupled Reactive LLE and Reactive Electrolyte Phase Equilibrium

## Suggested title

Replace staged reactive phase equilibrium with a production coupled native reactive LLE solver

## Purpose

Implement a production coupled reactive phase-equilibrium solver.

The required chain is:

```text
reaction residuals
→ phase-equilibrium residuals
→ electrolyte constraints when ions exist
→ one native coupled Ceres trust-region residual solve
→ implicit solved-state sensitivity diagnostics
→ generic Python reactive phase-equilibrium API
→ neutral reactive LLE and reactive electrolyte LLE benchmarks
```

This issue is complete only when chemical equilibrium and phase equilibrium are solved together in one native production residual problem.

## Current code that must be corrected

The implementation agent must inspect and explicitly address:

```text
src/epcsaft/equilibrium.py
src/epcsaft/reactive.py
src/epcsaft/reactive_speciation.py
src/epcsaft/reactive_staged.py
src/epcsaft/native/epcsaft_equilibrium.cpp
src/epcsaft/native/epcsaft_equilibrium.h
src/epcsaft/native/epcsaft_chemical_equilibrium.cpp
src/epcsaft/native/epcsaft_chemical_equilibrium.h
tests/equilibrium/reactive/
tests/api/reactive/
tests/native/equilibrium/
```

Known current problem:

```text
ReactivePhaseEquilibriumProblem is currently documented as reactive speciation followed by a generic phase-equilibrium route.
ReactivePhaseEquilibriumProblem.solve currently calls mixture.reactive_staged_equilibrium(...).
Existing staged reactive tests can pass without one coupled native reactive phase-equilibrium solve.
```

This issue must replace that production route. The staged route may remain only as an explicitly named initialization/comparison helper. It must not be the production implementation of `ReactivePhaseEquilibriumProblem`.

## Dependency gate

This issue depends on:

```text
#115 native activity speciation / reactive VLE merged through PR #124
#116 production electrolyte LLE solver completed in the same branch or already merged
#118 liquid-electrolyte derivative/regression support merged through PR #125
```

If #116 is not merged, this issue may be implemented in the same branch only if the #116 production solver is implemented first and tests pass before #117 work begins.

Before source edits:

```powershell
git fetch origin --prune
git merge --ff-only origin/main
git status --short
```

Record dependency status in:

```text
docs/goals/native-coupled-reactive-lle/notes/dependency_gate.md
```

## Stage 0 — Intake

Create:

```text
docs/goals/native-coupled-reactive-lle/notes/intake.md
```

It must include:

```text
origin/main SHA
current branch
confirmation that #115 and #118 are merged
confirmation that #116 production solver exists or is being implemented first in this branch
current ReactivePhaseEquilibriumProblem.solve route
current reactive_staged_equilibrium route
current native chemical-equilibrium route from #115
current neutral/electrolyte LLE production solver route from #116
current reactive LLE tests that accept staged behavior
chosen neutral reactive LLE benchmark fixture
chosen reactive electrolyte LLE benchmark fixture
```

No source edit may happen before this file is complete.

## Stage 1 — Native solver ownership

Create a native coupled reactive phase-equilibrium solver.

Preferred files:

```text
src/epcsaft/native/equilibrium/reactive_phase_equilibrium_problem.h
src/epcsaft/native/equilibrium/reactive_phase_equilibrium_problem.cpp
src/epcsaft/native/equilibrium/reactive_phase_equilibrium_residual.h
src/epcsaft/native/equilibrium/reactive_phase_equilibrium_residual.cpp
src/epcsaft/native/equilibrium/reactive_phase_equilibrium_solver.h
src/epcsaft/native/equilibrium/reactive_phase_equilibrium_solver.cpp
```

If a smaller patch must keep the implementation inside `epcsaft_equilibrium.cpp`, the code still must separate these responsibilities:

```text
problem construction
phase/species mask handling
variable pack/unpack
reaction residual evaluation
phase-equilibrium residual evaluation
charge/electroneutral residual evaluation
Jacobian evaluation
Ceres solve
diagnostic construction
```

## Stage 2 — Define coupled variables

The solver must represent all active phases and reactions in one unknown vector.

Required concepts:

```text
phase amounts or phase fractions
species amounts or transformed phase compositions
reaction extents or transformed reaction coordinates
electroneutral ion-combination variables where charged species exist
density/pressure closure variables if required
```

A valid implementation pattern is:

```text
u_beta      = unconstrained phase-fraction variables
u_x[p]      = unconstrained phase-composition variables
u_xi[p,r]   = reaction-coordinate variables or direct species-amount variables
```

The physical state is unpacked as:

```text
phase fractions beta_p
phase species mole fractions x_i,p
phase species amounts n_i,p
reaction extents xi_r,p
```

All mole numbers must be nonnegative. Use log, softplus, or softmax-style transforms. Do not rely on clamping after an infeasible solve as the production method.

## Stage 3 — Coupled residual vector

Implement:

```cpp
ReactivePhaseResidualValue evaluate_reactive_phase_equilibrium_residual(
    const ReactivePhaseVariables& u,
    const ReactivePhaseEquilibriumProblemNative& problem
);
```

Required residual blocks:

### Element balance

For element matrix `A` and total element vector `b`:

```math
R^{elem}=A\sum_p n_p-b
```

### Species material balance when species amounts are used directly

If direct species-balance form is used:

```math
R_i^{mat} =
n_i^{feed}
-
\sum_p n_{i,p}
+
\sum_r \nu_{ir}\xi_r
```

If element-basis variables enforce material balance, report this as a diagnostic instead.

### Reaction equilibrium

For reaction `r` in reacting phase `p`:

```math
R_{r,p}^{rxn}
=
\sum_i \nu_{ir}\ln a_{i,p}
-
\ln K_r
```

Activities must come from the ePC-SAFT state evaluated inside the coupled residual.

### Neutral interphase equilibrium

For neutral species `i` allowed in phases `alpha` and `beta`:

```math
R_i^{phase}
=
\ln f_i^\alpha
-
\ln f_i^\beta
```

### Ionic interphase equilibrium

For each independent electroneutral ion-combination vector `b_q`:

```math
\sum_i z_i b_{iq}=0
```

```math
R_q^{ion}
=
\sum_i b_{iq}
\left(
\mu_i^\alpha-\mu_i^\beta
\right)/RT
```

### Charge balance per liquid phase

For each liquid phase with charged species:

```math
R_p^{charge}
=
\sum_i z_i n_{i,p}
```

If the transform enforces charge by construction, report this as a diagnostic and assert it in tests.

### Phase normalization

If composition variables do not enforce normalization:

```math
R_p^{norm}
=
\sum_i x_{i,p}-1
```

### Closure residuals

Add density/pressure closure residuals only where the chosen state variables require them. If density roots are solved inside property evaluation, use implicit sensitivity diagnostics.

## Stage 4 — Jacobian and solved-state sensitivities

The production residual solve must have a real Jacobian.

Expected implementation:

```text
CppAD residual Jacobian with respect to transformed variables
analytic block derivatives where simpler
implicit sensitivity for density, association, speciation, and phase-split solved states
```

Do not differentiate through iteration histories.

Required solved-state relation:

```math
F(u,\theta)=0
```

```math
\frac{du}{d\theta}
=
-F_u^{-1}F_\theta
```

The PR must report:

```text
jacobian_available = true
jacobian_backend = cppad or analytic_implicit/cppad_implicit
solved_state_sensitivity_available = true
solved_state_sensitivity_backend
```

## Stage 5 — Ceres production solve

Use Ceres for the coupled residual equations.

Do not use Ceres as a direct Gibbs minimizer. Do not use it as a substitute for a constrained Gibbs NLP. Use it as the native trust-region nonlinear residual solver after constraints are represented through transforms and residual blocks.

Required cost-function form:

```cpp
class ReactivePhaseEquilibriumCostFunction final : public ceres::CostFunction {
public:
    explicit ReactivePhaseEquilibriumCostFunction(ReactivePhaseResidualEvaluator evaluator);

    bool Evaluate(
        double const* const* parameters,
        double* residuals,
        double** jacobians
    ) const override;
};
```

Required solver diagnostics for accepted production results:

```text
solver_backend = ceres
solver_method = ceres_trust_region_coupled_reactive_phase_equilibrium
ceres_trust_region_strategy
ceres_linear_solver
ceres_termination_type
ceres_initial_cost
ceres_final_cost
ceres_iteration_count
jacobian_backend
derivative_backend
```

Accepted production results must not report:

```text
reactive_staged_equilibrium as production route
staged chemical-then-phase as production result
not_required derivative status for the accepted nonlinear solve
hand-coded simplex/Powell/Nelder-Mead solve path
```

## Stage 6 — Use #115 and #116 only as subcomponents

Allowed:

```text
use #115 native activity speciation as an initial guess
use #116 electrolyte LLE solver as an initial phase-split guess
use TPD/g-hat stability analysis to choose initial phase count and seeds
use staged workflow as comparison in tests
```

Not allowed:

```text
speciation result accepted first and phase equilibrium solved after as final result
phase equilibrium accepted first and reaction equilibrium solved after as final result
staged route presented as coupled reactive phase equilibrium
```

The final result must come from one coupled residual solve and one consistent solved state.

## Stage 7 — Python API

Update the generic API so:

```python
ReactivePhaseEquilibriumProblem(...).solve(mixture)
```

routes to the native coupled production solve for reactive LLE / reactive electrolyte LLE cases.

If compatibility with the staged route is retained, expose it as a clearly separate helper. Example:

```python
StagedReactivePhaseEquilibriumProblem(...)
```

or an internal/private compatibility route.

The public production API must remain generic:

```python
ReactivePhaseEquilibriumProblem(...)
ReactiveLLEProblem(...)
ReactiveElectrolyteLLEProblem(...)
equilibrium(problem)
```

No application-specific APIs.

## Stage 8 — Benchmarks

Required benchmark A:

```text
neutral reactive LLE benchmark with real phase split and reaction equilibrium
```

Use an esterification-style fixture based on the Ascani 2023 pattern. The benchmark must have repo-contained data and numeric tolerances.

Required benchmark B:

```text
reactive electrolyte LLE benchmark with charged species and one reaction
```

Use a generic ion-exchange or reactive transfer fixture. The benchmark must check thermodynamic residuals. It must not compute or expose extraction efficiency, distribution coefficient, selectivity, or any downstream application metric in package public API.

Benchmark diagnostics must include:

```text
reaction_residual_norm
phase_equilibrium_residual_norm
ionic_equilibrium_residual_norm where ions exist
material_balance_norm
element_balance_norm
phase_charge_balance_norm where ions exist
phase_distance
phase_compositions
phase_amounts
reaction_extents
solver_backend
solver_method
jacobian_backend
derivative_backend
```

## Stage 9 — Tests

Required new or updated tests:

```text
tests/native/equilibrium/test_reactive_phase_equilibrium_ceres_solver.py
tests/native/equilibrium/test_reactive_phase_equilibrium_residual_jacobian.py
tests/equilibrium/reactive/test_reactive_lle_coupled_solver.py
tests/equilibrium/reactive/test_reactive_electrolyte_lle_coupled_solver.py
tests/api/reactive/test_reactive_phase_equilibrium_problem_routes_native.py
tests/api/reactive/test_staged_reactive_route_not_production.py
```

Tests must fail if:

```text
ReactivePhaseEquilibriumProblem.solve calls reactive_staged_equilibrium for the production route
accepted diagnostics report a staged route
reaction residual and phase residual are computed from different solved states
phase charge balance is missing for charged phases
Jacobian is absent for accepted production coupled solve
```

## Stage 10 — Validation commands

Run:

```powershell
uv run python scripts/dev/build_epcsaft.py --clean --enable-ceres --enable-cppad
uv run python run_pytest.py tests/native/equilibrium -q
uv run python run_pytest.py tests/equilibrium/reactive -q
uv run python run_pytest.py tests/api/reactive -q
uv run python scripts/dev/validate_project.py quick
uv run python scripts/dev/validate_project.py docs
git diff --check
```

Also run targeted route guards:

```powershell
rg "ReactivePhaseEquilibriumProblem|reactive_staged_equilibrium|staged chemical|native_derivative_free_nelder_mead|not_required.*phase" src tests docs
```

Any remaining match must be explicitly justified as compatibility-only, initialization-only, or a test that proves the old route is not used as production.

## Definition of done

This issue is complete only when every line is true:

```text
ReactivePhaseEquilibriumProblem production route is native coupled solve
reaction and phase residuals are evaluated in one coupled solved state
solver_backend is ceres for accepted production benchmarks
Jacobian is analytic / CppAD / implicit, not an approximate derivative route
staged reactive route is not the production result
#115 speciation and #116 LLE are used only as initialization/subcomponents
neutral reactive LLE benchmark passes
reactive electrolyte LLE benchmark passes
phase charge balance is enforced for charged phases
element/material balance is enforced
reaction residual norm is below tolerance
phase-equilibrium residual norm is below tolerance
Python API remains generic
no application-specific metric API is added
all validation commands pass
```
