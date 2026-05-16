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


# Issue #116 Replacement Body — Native Electrolyte LLE with Distributed Ions and Production Solver

## Suggested title

Replace hand-coded LLE minimizers with production Ceres-based native electrolyte LLE solver

## Purpose

Implement a production electrolyte LLE solver for distributed-ion mixed-solvent systems.

The required chain is:

```text
native transformed electrolyte variables
→ TPD/g-hat stability and seed generation
→ Ceres trust-region residual solve with CppAD/implicit Jacobian
→ distributed-ion phase equilibrium residuals
→ Python ElectrolyteLLEProblem API
→ Ascani/Khudaida-style benchmark proof
```

The issue is complete only when a distributed-ion electrolyte LLE benchmark solves through the production solver and reports residuals, charge balance, material balance, phase distance, and derivative provenance from the solved state.

## Current code that must be corrected

The implementation agent must inspect and explicitly address these current paths before writing source code:

```text
src/epcsaft/native/epcsaft_equilibrium.cpp
src/epcsaft/native/epcsaft_equilibrium.h
src/epcsaft/equilibrium.py
src/epcsaft/electrolyte_lle.py
tests/equilibrium/electrolyte/
tests/native/equilibrium/
```

Known current problems to replace:

```text
minimize_lle_residual_variables(...)
polish_formula_tpd_variables(...)
nelder_mead_variables(...)
newton_step(...) currently throwing because electrolyte LLE residual sensitivities are not implemented
accepted diagnostics that label a production solve as native_transformed_newton while the real Newton step is not implemented
accepted diagnostics that label neutral LLE as native_derivative_free_nelder_mead
Reactive/electrolyte tests that accept candidate diagnostics rather than a production solver
```

The PR must remove the possibility that an accepted electrolyte LLE result is produced by hand-coded simplex minimization. If any old helper remains, it must be renamed as seed-only/test-only and must not appear in accepted production diagnostics.

## Literature and thermodynamic target

This issue is Ascani-style in the following specific sense:

```text
ion species are retained explicitly
electroneutrality is enforced per liquid phase
mixed salts and common ions are supported
ePC-SAFT advanced fugacities/chemical potentials drive phase equilibrium
TPD/stability and g-hat/Gibbs checks guide phase count and candidate selection
phase compositions are determined from phase-equilibrium equations
```

This issue is not a direct Gibbs-energy NLP and not a RAND implementation. IPOPT/RAND may be future backends, but they are not required here.

The residual Helmholtz basis is:

```math
a^{res}=a^{HC}+a^{disp}+a^{assoc}+a^{DH}+a^{Born}
```

Electrolyte LLE must use phase chemical potentials or equivalent log fugacity forms derived from this ePC-SAFT state.

## Dependency gate

Start from current `origin/main` after these are merged:

```text
#114 via PR #123
#115 via PR #124
#118 via PR #125
```

Before source edits, run:

```powershell
git fetch origin --prune
git merge --ff-only origin/main
git status --short
```

If the branch cannot fast-forward to current `origin/main`, stop and rebase/merge cleanly before implementation.

## Stage 0 — Intake

Create:

```text
docs/goals/native-electrolyte-lle-distributed-ions/notes/intake.md
```

It must include:

```text
origin/main SHA
current branch
confirmation that PR #123, #124, and #125 merge commits are in ancestry
current neutral LLE production route and diagnostics
current electrolyte LLE production route and diagnostics
all current functions that do hand-coded simplex/Powell/Nelder-Mead-style minimization
current newton_step(...) behavior
current electrolyte basis implementation
current TPD/g-hat seed implementation
current Python ElectrolyteLLEProblem route
current tests that would still pass with the old non-production route
chosen Ascani-style fixture
chosen salting-out fixture
```

No source edit may happen before this file is complete.

## Stage 1 — Refactor native solver ownership

Create or update a clear solver module. Preferred layout:

```text
src/epcsaft/native/equilibrium/electrolyte_lle_basis.h
src/epcsaft/native/equilibrium/electrolyte_lle_basis.cpp
src/epcsaft/native/equilibrium/electrolyte_lle_residual.h
src/epcsaft/native/equilibrium/electrolyte_lle_residual.cpp
src/epcsaft/native/equilibrium/electrolyte_lle_solver.h
src/epcsaft/native/equilibrium/electrolyte_lle_solver.cpp
```

If the repo structure makes a smaller patch safer, implementation may stay inside `epcsaft_equilibrium.cpp`, but the PR must still separate these responsibilities in code comments and helper functions:

```text
basis construction
variable pack/unpack
residual evaluation
Jacobian evaluation
Ceres solve
diagnostic construction
benchmark fixture routing
```

## Stage 2 — Build a real distributed-ion basis

The public problem uses explicit species, including ions.

The internal solver may use transformed formula variables, but only as a coordinate system. The result and residual diagnostics must always report explicit species mole fractions.

Required basis behavior:

```text
supports at least two neutral species plus cations and anions
supports common-ion mixed-salt systems
does not collapse the public model into salt-only components
enforces phase electroneutrality by construction or by explicit residual
documents the exact independent ion-combination basis used
```

For charges `z_i`, every independent ion-combination vector `b_q` must satisfy:

```math
\sum_i z_i b_{iq}=0
```

The PR must add a basis report in diagnostics:

```text
variable_model
basis_rank
basis_vectors
charged_species_indices
neutral_species_indices
phase_charge_balance_aq
phase_charge_balance_org
```

## Stage 3 — Define transformed variables

For two liquid phases, use variables that keep the solve feasible:

```text
u_beta       = unconstrained phase-fraction coordinate
u_org        = unconstrained organic/formula composition coordinates
x_org(u)     = explicit organic-phase species mole fractions
x_aq(u)      = explicit aqueous-phase species mole fractions from material balance
beta_org(u)  = physical organic phase fraction
```

A valid pattern is:

```math
\beta = \sigma(u_\beta)
```

```math
q^{org} = \operatorname{softmax}(u_{org})
```

```math
q^{aq}_k =
\frac{q^{feed}_k-\beta_q q^{org}_k}{1-\beta_q}
```

where the implementation must reject infeasible candidates before they can be reported as accepted.

If a different transform is used, document it in:

```text
docs/goals/native-electrolyte-lle-distributed-ions/notes/variable_transform.md
```

## Stage 4 — Residual vector

Implement a residual evaluator:

```cpp
ElectrolyteLLEResidualValue evaluate_electrolyte_lle_residual(
    const ElectrolyteLLEVariables& u,
    const ElectrolyteLLEProblemNative& problem
);
```

Required residual blocks:

### Neutral species phase equilibrium

For neutral species `i` present in both liquid phases:

```math
R_i^{neutral} =
\left[\ln x_i^{org}+\ln \varphi_i^{org}\right]
-
\left[\ln x_i^{aq}+\ln \varphi_i^{aq}\right]
```

### Ionic electroneutral-combination equilibrium

For every independent electroneutral ion-combination vector `b_q`:

```math
R_q^{ion} =
\sum_i b_{iq}
\left[
\mu_i^{org}-\mu_i^{aq}
\right]/RT
```

Equivalent log-fugacity form is acceptable if documented and tested.

### Material balance diagnostic

If the transform enforces material balance by construction, material balance is a diagnostic. If it does not, it must be a residual block.

```math
R_i^{mat} =
z_i^{feed}
-
\sum_p \beta_p x_{i,p}
```

### Charge balance diagnostic or residual

If the transform enforces charge balance by construction, charge balance is a diagnostic. If it does not, it must be a residual block.

```math
R_p^{charge} =
\sum_i z_i x_{i,p}
```

### Scaling

All residuals must be scaled. Add explicit scaling constants and tests so no one component dominates only because of units.

Required diagnostic norms:

```text
neutral_fugacity_residual_norm
ionic_equilibrium_residual_norm
material_balance_norm
phase_charge_balance_norm
phase_distance
scaled_solver_residual_norm
unscaled_solver_residual_norm
```

## Stage 5 — Jacobian evaluation

Implement the production Jacobian. No derivative-approximation Jacobian.

Expected implementation:

```cpp
struct ResidualAndJacobian {
    std::vector<double> residual;
    std::vector<double> jacobian_row_major;
    int rows;
    int cols;
};
```

Use CppAD over the transformed variable vector `u` wherever the residual expression is explicit.

For solved internal states inside each residual evaluation:

```text
density roots
association site fractions
other internal implicit states
```

use existing analytic or CppAD implicit sensitivity paths. Do not differentiate through iteration loops.

If a CppAD path cannot be made for a required residual, the PR is incomplete. Do not mark the issue ready.

## Stage 6 — Ceres production solve

Use Ceres as the production nonlinear residual solver.

Do not use Ceres as a direct Gibbs minimizer. Use it to solve the transformed thermodynamic residual equations.

Required cost-function shape:

```cpp
class ElectrolyteLLECostFunction final : public ceres::CostFunction {
public:
    explicit ElectrolyteLLECostFunction(ElectrolyteLLEResidualEvaluator evaluator);

    bool Evaluate(
        double const* const* parameters,
        double* residuals,
        double** jacobians
    ) const override;
};
```

Inside `Evaluate`:

```text
read transformed variables u
evaluate residual vector R(u)
if jacobians requested:
    evaluate dR/du through CppAD / implicit sensitivity
    copy row-major Jacobian into Ceres output
return false only when the candidate is numerically invalid
```

Required Ceres options:

```cpp
ceres::Solver::Options options;
options.linear_solver_type = ceres::DENSE_QR;
options.max_num_iterations = equilibrium_options.max_iterations;
options.function_tolerance = equilibrium_options.tolerance;
options.gradient_tolerance = equilibrium_options.tolerance;
options.parameter_tolerance = equilibrium_options.tolerance;
```

Dogleg or sparse linear solvers may be used if tests prove better behavior.

Accepted production diagnostics must include:

```text
solver_backend = ceres
solver_method = ceres_trust_region_residual_solve
ceres_trust_region_strategy
ceres_linear_solver
ceres_termination_type
ceres_initial_cost
ceres_final_cost
ceres_iteration_count
jacobian_backend
derivative_backend
jacobian_available = true
```

Accepted production diagnostics must not include the old accepted solver labels for the production result:

```text
native_derivative_free_nelder_mead
native_transformed_newton when newton_step is not implemented
not_required derivative status for the accepted nonlinear solve
```

## Stage 7 — TPD and g-hat usage

Keep TPD and g-hat only in their correct roles.

Allowed:

```text
stability precheck
trial phase generation
phase-count decision support
seed selection
post-solve stability check
acceptance check that split is thermodynamically favored
```

Not allowed:

```text
TPD-only accepted LLE result
g-hat-only accepted LLE result
candidate with small g-hat but failed phase-equilibrium residuals
```

Use this diagnostic expression for the g-hat-style check:

```math
\hat{g}
=
RT
\sum_p \beta_p
\sum_i x_{i,p}\ln f_{i,p}
```

The PR must report:

```text
g_hat_feed
g_hat_split
g_hat_delta
tpd_min
tpd_trial_count
post_solve_stability_checked
```

## Stage 8 — Python API

Keep the public API generic:

```python
ElectrolyteLLEProblem(...)
equilibrium(problem)
mixture.electrolyte_lle_tp(...)
```

The Python layer may validate inputs and format results. It must not own the production solve loop.

Required Python result fields:

```text
phases
diagnostics
attempt_diagnostics
phase_labels
split_detected
stable
```

## Stage 9 — Benchmarks and tests

Required benchmark A:

```text
Ascani-style distributed-ion electrolyte LLE fixture with mixed solvent and at least one salt
```

Required benchmark B:

```text
water + ethanol + isobutanol + NaCl salting-out LLE fixture or equivalent repo-contained quaternary salting-out electrolyte LLE fixture
```

Benchmark checks must assert:

```text
two liquid phases exist when expected
each phase has finite density
each phase composition sums to one
each phase charge balance is below tolerance
material balance is below tolerance
neutral phase-equilibrium residual norm is below tolerance
ionic equilibrium residual norm is below tolerance
g_hat_delta supports accepted split
phase_distance is above anti-trivial threshold
solver_backend is ceres
jacobian_available is true
```

Required new or updated tests:

```text
tests/native/equilibrium/test_electrolyte_lle_ceres_solver.py
tests/native/equilibrium/test_electrolyte_lle_residual_jacobian.py
tests/equilibrium/electrolyte/test_distributed_ion_lle_production_solver.py
tests/equilibrium/electrolyte/test_salting_out_lle_benchmark.py
tests/api/equilibrium/test_electrolyte_lle_problem_production_route.py
```

Tests must fail if accepted production diagnostics report the old hand-coded simplex route.

## Stage 10 — Validation commands

Run:

```powershell
uv run python scripts/dev/build_epcsaft.py --clean --enable-ceres --enable-cppad
uv run python run_pytest.py tests/native/equilibrium -q
uv run python run_pytest.py tests/equilibrium/electrolyte -q
uv run python run_pytest.py tests/equilibrium/core -q
uv run python run_pytest.py tests/api/equilibrium -q
uv run python scripts/dev/validate_project.py quick
uv run python scripts/dev/validate_project.py docs
git diff --check
```

Also run targeted guard searches:

```powershell
rg "native_derivative_free_nelder_mead|not_required.*phase split solve|newton_step.*missing sensitivity|reactive_staged_equilibrium" src tests docs
```

Any remaining match must be explicitly justified as a rejected old path, compatibility-only route, or test asserting the old path is not used.

## Definition of done

This issue is complete only when every line is true:

```text
accepted electrolyte LLE result is solved by Ceres trust-region residual solve
production residuals use explicit ePC-SAFT fugacity/chemical-potential evaluations
production Jacobian is analytic / CppAD / implicit, not an approximate derivative route
old hand-coded simplex route cannot produce an accepted production result
newton_step missing-sensitivity path is removed or no longer reachable for accepted production
distributed-ion phase variables are explicit in public result and diagnostics
each liquid phase is electroneutral
mixed-solvent and common-ion mixed-electrolyte cases are supported
TPD/g-hat are stability and seed/acceptance tools, not sole production solve
Ascani-style benchmark passes
salting-out benchmark passes
Python API remains generic
all validation commands pass
```
