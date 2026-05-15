# PR #126 Repair Audit and Completion Contract

Repository: `tannerpolley/ePC-SAFT`  
Target branch for repair agent: create a fresh branch from current `origin/main`  
Primary files to inspect: native equilibrium C++, Python equilibrium API, reactive convention handling, and tests added or edited by PR #126

This document is an implementation audit and a completion contract. It is intended to be attached to the prompt for the next Codex agent. The agent must read this entire document before planning, coding, testing, or opening a PR.

---

## 0. Mission

PR #126 merged real work, but it did not finish the job. The repair agent must fix the incomplete and incorrect parts of PR #126 without reducing scope, hiding failures behind diagnostics, replacing missing logic with toy tests, or closing with a checklist-only audit.

The repair PR is complete only when the package has:

```text
native C++ production equilibrium solvers
honest Ceres acceptance gates
correct reaction standard-state residuals
real derivative/Jacobian provenance
generic Python API routing
no accepted hand-coded simplex LLE route
real source-backed benchmarks where the issues requested literature proof
strict tests that fail when old shortcuts reappear
```

The repair PR must not merely update comments, docs, labels, or diagnostics. It must change the implementation and prove the changes with tests.

---

## 1. Source facts that triggered this audit

### 1.1 PR #126 closed the issues but overclaimed completion

PR #126 is titled:

```text
Native LLE reactive production solvers
```

It merged into `main` and closed both issues #116 and #117. Its PR body claimed that it replaced accepted electrolyte LLE behavior with native Ceres trust-region residual solving, routed reactive LLE and reactive electrolyte LLE through a native coupled Ceres solve, and added validation/audit receipts.

GitHub source:

```text
https://github.com/tannerpolley/ePC-SAFT/pull/126
```

Merge commit reported by GitHub:

```text
869e3354ddc0b52075ddc9efe687b34d6aa98316
```

The next agent must not accept the merge status as evidence of correctness. It must audit the code.

### 1.2 Codex review found two P1 bugs immediately after merge

PR #126 received two post-merge P1 review comments.

#### P1-A: reaction standard states are passed but ignored

Review comment summary:

```text
The reactive residual path treats every reaction as sum(nu_i * ln a_i) - ln K = 0 using only ln_activity.
The new production route passes reaction_standard_states for conventions like ideal_mole_fraction and concentration.
Those codes are never applied, so non-activity conventions are silently solved on the wrong thermodynamic basis.
```

GitHub source:

```text
https://github.com/tannerpolley/ePC-SAFT/pull/126#discussion_r3249116078
```

#### P1-B: reactive Ceres solve is marked accepted even if Ceres fails

Review comment summary:

```text
After ceres::Solve, the function always returns a normal result and unconditionally sets ceres_accepted_solve=true.
It does not check summary.IsSolutionUsable() or termination type.
No-convergence or failure cases can therefore propagate accepted two-phase outputs.
```

GitHub source:

```text
https://github.com/tannerpolley/ePC-SAFT/pull/126#discussion_r3249116082
```

The repair PR must address both P1s in code and tests.

---

## 2. Verified current code problems

The repair agent must confirm these from the current checkout before editing. The files and function names below are the expected current state after PR #126.

### 2.1 `reaction_standard_states` are ignored in native residual equations

Current file:

```text
src/epcsaft/native/equilibrium/reactive_phase_equilibrium_problem.cpp
```

Current function of concern:

```cpp
std::vector<double> reaction_residuals(
    const std::vector<double>& stoichiometry_row_major,
    int reaction_rows,
    int ncomp,
    const std::vector<double>& log_equilibrium_constants,
    const std::vector<double>& ln_activity
)
```

Current behavior:

```text
R_r = sum_i nu_ri * ln_activity_i - lnK_r
```

The function does not receive or use `reaction_standard_states`. The native evaluator accepts `reaction_standard_states`, validates the vector length, stores it in diagnostics, and passes it through the Ceres cost function, but the actual reaction residual formula does not branch by convention.

This is not acceptable.

### 2.2 Ceres acceptance is unconditional in native reactive phase equilibrium

Current file:

```text
src/epcsaft/native/equilibrium/reactive_phase_equilibrium_problem.cpp
```

Current behavior after `ceres::Solve`:

```cpp
result.diagnostics_bool["ceres_accepted_solve"] = true;
result.diagnostics_bool["jacobian_available"] = true;
result.diagnostics_bool["derivative_available"] = true;
result.diagnostics_bool["solved_state_sensitivity_available"] = true;
return result;
```

The current route builds and returns a normal two-phase result even when the Ceres termination is not usable. This must be replaced with strict solver and thermodynamic acceptance gates.

### 2.3 Electrolyte tests allow no-convergence as accepted production evidence

Current test of concern:

```text
tests/native/equilibrium/test_electrolyte_lle_ceres_solver.py
```

Current bad assertion pattern:

```python
assert diagnostics["ceres_termination_type"] in {"convergence", "no_convergence"}
```

An accepted production equilibrium test must not allow `no_convergence`.

### 2.4 Old hand-coded neutral LLE minimizer still exists and is still production

Current file:

```text
src/epcsaft/native/epcsaft_equilibrium.cpp
```

Current functions of concern:

```cpp
minimize_lle_residual_variables(...)
solve_lle_attempt(...)
lle_two_phase_result(...)
```

Current accepted neutral LLE behavior still calls the hand-coded simplex minimizer and can report:

```text
nonlinear_solver = native_derivative_free_nelder_mead
derivative_backend = not_applicable
derivative_status = not_required
```

This conflicts with the original #116 instruction to replace hand-coded LLE minimizers as accepted production solvers. Do not classify this as “neutral-only therefore acceptable” unless a new issue explicitly narrows #116 after maintainer approval. For this repair, implement the replacement.

### 2.5 TPD seed polishing still uses a hand-coded simplex-like optimizer

Current file:

```text
src/epcsaft/native/epcsaft_equilibrium.cpp
```

Current function of concern:

```cpp
polish_formula_tpd_variables(...)
```

This may remain only if it is strictly seed-generation support and cannot return an accepted equilibrium. Preferred repair: remove hand-coded optimizer behavior from production source or replace it with deterministic non-optimizing seed generation or a packaged solver path. If it remains, tests must prove it is not used as a production equilibrium solve.

### 2.6 Failure diagnostics currently look like successful production derivatives

Current file:

```text
src/epcsaft/native/epcsaft_equilibrium.cpp
```

Current function of concern:

```cpp
electrolyte_lle_failure_result(...)
```

Problem behavior:

```text
acceptance_gate = predictive_solve_failed
solver_backend = ceres
solver_method = ceres_trust_region_residual_solve
jacobian_backend = cppad_implicit
derivative_backend = cppad_implicit
jacobian_available = true
derivative_available = true
```

A failed solve may report which solver was attempted. It must not report availability or accepted derivative status as though a valid accepted equilibrium state exists.

Required distinction:

```text
solver_attempted = ceres
solver_attempt_result = failed
accepted_solver_backend = none
accepted_derivative_backend = none
jacobian_available_for_accepted_state = false
solution_accepted = false
```

Use exact field names only after checking existing diagnostics style.

### 2.7 Benchmarks added by PR #126 are too weak

#### Electrolyte salting-out benchmark is not the requested Khudaida system

Current test:

```text
tests/equilibrium/electrolyte/test_salting_out_lle_benchmark.py
```

Current fixture:

```python
species = ["H2O", "Butanol", "Na+", "Cl-"]
mix = ePCSAFTMixture.from_dataset("2022_Ascani", species, feed, 298.15)
```

This is not the requested Khudaida 2026 quaternary system:

```text
water + ethanol + isobutyl alcohol + NaCl
293.15 K and 313.15 K
5 wt% and 10 wt% NaCl
```

#### Neutral reactive LLE benchmark is a model-consistent toy

Current test:

```text
tests/equilibrium/reactive/test_reactive_lle_coupled_solver.py
```

Current fixture is effectively:

```text
Methanol ⇌ Cyclohexane
```

with a log-equilibrium constant chosen from the same current model state.

This is a smoke test, not the requested Ascani 2023 esterification-style reactive LLE benchmark.

#### Reactive electrolyte benchmark is not actually ion-coupled

Current test:

```text
tests/equilibrium/reactive/test_reactive_electrolyte_lle_coupled_solver.py
```

Current fixture reaction:

```text
H2O ⇌ Butanol
```

Ions are spectators. That does not prove a reactive electrolyte LLE coupled through charged species.

### 2.8 No shared reusable residual-solver architecture was created

PR #126 added separate one-off Ceres patterns in large files. The repair agent must introduce or finish a shared native residual-solver abstraction so every production residual solve uses the same acceptance, diagnostics, and Jacobian contract.

Expected reusable shape:

```text
src/epcsaft/native/equilibrium/residual_solver.h
src/epcsaft/native/equilibrium/residual_solver.cpp
```

If the repo already added an equivalent after this document was created, use that equivalent. Otherwise, add it.

---

## 3. Scientific and package scope that must guide the repair

### 3.1 Package identity

`epcsaft` is a general-purpose thermodynamic package.

Allowed public concepts:

```python
ParameterSet
State
EquilibriumProblem
ReactiveEquilibriumProblem
ElectrolyteLLEProblem
ReactiveLLEProblem
RegressionProblem
TargetRow
TargetDataset
model.equilibrium(problem)
epcsaft.equilibrium(model, problem)
epcsaft.regress_parameters(problem)
```

Forbidden public API concepts:

```python
fit_lithium_extraction_parameters(...)
screen_lithium_extractants(...)
fit_mea_absorption(...)
calculate_extraction_efficiency(...)
calculate_distribution_coefficient(...)
calculate_selectivity(...)
fit_absorption_column(...)
```

Downstream projects compute application metrics from generic package outputs.

### 3.2 Thermodynamic basis

The ePC-SAFT advanced residual Helmholtz basis is:

```math
a^{res}=a^{HC}+a^{disp}+a^{assoc}+a^{DH}+a^{Born}
```

All equilibrium residuals must derive activities, fugacities, or chemical potentials from the ePC-SAFT state. Do not invent surrogate activity formulas.

### 3.3 Electrolyte LLE target

The electrolyte LLE solver must follow the Ascani-style distributed-ion purpose:

```text
explicit ions in the public species model
mixed solvents
mixed salts / common ions
per-phase electroneutrality
ePC-SAFT fugacity or electrochemical-potential equality
TPD and Gibbs-like checks as stability/initialization/acceptance support
production equilibrium from residual equations, not seed minimization alone
```

### 3.4 Reactive LLE target

The reactive LLE solver must follow the Ascani 2023-style purpose:

```text
chemical equilibrium and phase equilibrium in one coupled state
reaction residuals and phase residuals evaluated on the same current state
not staged speciation followed by phase equilibrium
not a two-species artificial conversion test as literature proof
```

Staged routes may initialize or compare. They must not be production completion proof.

---

## 4. Required implementation repairs

The repair agent must implement all sections below. Do not split them into “future work” unless the PR remains draft/open and does not claim completion.

---

### Repair A — Reaction standard-state conventions

#### A1. Inspect existing convention definitions

Before coding, inspect:

```text
src/epcsaft/reactive_speciation.py
src/epcsaft/reactive.py
src/epcsaft/equilibrium.py
src/epcsaft/bindings.cpp
src/epcsaft/native/epcsaft_chemical_equilibrium.cpp
src/epcsaft/native/equilibrium/reactive_phase_equilibrium_problem.cpp
```

Find the exact mapping for:

```text
mole_fraction_activity
ideal_mole_fraction
concentration
any other currently accepted convention
```

Do not guess convention codes or formulas. Use the existing `ReactionConstantConvention` / `native_standard_state_code` definitions as the source of truth. If the existing definitions are ambiguous, fix them in Python and native code together.

#### A2. Add a native standard-state residual helper

Add a helper that computes the per-species log reaction term for a given phase and reaction convention.

Expected conceptual shape:

```cpp
std::vector<double> reaction_log_terms_for_standard_state(
    const ReactivePhaseState& state,
    int standard_state_code,
    double temperature,
    double pressure,
    double concentration_reference,
    double floor
);
```

Then compute:

```math
R_{r,p} = \sum_i \nu_{ri} \, \Lambda_{i,p}^{(standard\ state)} - \ln K_r
```

Where `Lambda` is not always `ln(x_i) + ln(phi_i)`. The correct expression depends on the convention.

#### A3. Required convention behavior

The agent must implement exact behavior according to the existing package definitions. At minimum, tests must cover:

```text
mole_fraction_activity
ideal_mole_fraction
concentration
```

The expected values must be calculated independently in the tests, not copied from the native result.

#### A4. Native and Python diagnostics

Diagnostics must include:

```text
reaction_standard_state_codes
reaction_constant_conventions
reaction_residual_basis
reaction_residual_standard_state_applied = true
```

Do not report `reaction_residual_standard_state_applied=true` unless the native residual actually applies the codes.

#### A5. Required tests

Add tests that fail before the repair:

```text
tests/native/equilibrium/test_reactive_phase_equilibrium_standard_states.py
tests/equilibrium/reactive/test_reactive_lle_standard_state_conventions.py
```

Minimum test cases:

1. Same stoichiometry and same K, but different standard states, must produce different residuals when the composition/density makes the conventions different.
2. A native residual surface request with `reaction_standard_states=[ideal_mole_fraction_code]` must match a hand-calculated ideal mole-fraction residual.
3. A native residual surface request with `reaction_standard_states=[concentration_code]` must match a hand-calculated concentration-basis residual.
4. Existing mole-fraction activity behavior must remain correct.

---

### Repair B — Strict Ceres acceptance gates

#### B1. Shared acceptance policy

Add a shared acceptance helper used by both electrolyte LLE and reactive LLE.

Expected conceptual shape:

```cpp
struct SolverAcceptanceResult {
    bool solver_usable;
    bool termination_accepted;
    bool residuals_accepted;
    bool physical_gates_accepted;
    bool accepted;
    std::string rejection_reason;
};

SolverAcceptanceResult evaluate_ceres_acceptance(
    const ceres::Solver::Summary& summary,
    const ResidualEvaluation& final_eval,
    const PhysicalGateEvaluation& gates,
    const EquilibriumOptionsNative& options
);
```

Accepted production solve requires all:

```text
summary.IsSolutionUsable()
accepted termination type
all required residual norms <= tolerance
all material/element balance norms <= tolerance
all charge balance norms <= tolerance where charged species exist
phase distance gate satisfied for split result
phase fractions physical
density/pressure closure valid
all reported phase compositions finite, nonnegative, and normalized
```

#### B2. Do not accept these as production convergence

Forbidden for accepted production result:

```text
NO_CONVERGENCE
FAILURE
USER_FAILURE
result with residual norm above tolerance
result with only phase_distance proving split
result with missing final residual re-evaluation
```

`USER_SUCCESS` may only be accepted if the implementation explicitly uses a user callback to stop after all physical gates are satisfied. If no such callback exists, do not accept `USER_SUCCESS`.

#### B3. Reactive route must throw or return failed state honestly

In `reactive_phase_equilibrium_native(...)`, after `ceres::Solve`:

1. Re-evaluate final residuals.
2. Compute all physical gates.
3. If acceptance fails:
   - Do not return two accepted phases.
   - Do not set `ceres_accepted_solve=true`.
   - Return a failed result or raise a native solution error following existing package patterns.
   - Include diagnostics explaining the failure.

#### B4. Required tests

Add tests:

```text
tests/native/equilibrium/test_reactive_phase_equilibrium_ceres_acceptance.py
tests/equilibrium/reactive/test_reactive_phase_equilibrium_failure_gates.py
tests/native/equilibrium/test_electrolyte_lle_ceres_acceptance.py
```

Minimum test cases:

1. Force `max_iterations=0` or `1` so Ceres cannot solve.
2. Assert no accepted result is returned.
3. Assert `ceres_accepted_solve` is false or absent.
4. Assert `solver_usable` is false when Ceres says unusable.
5. Assert no phases are returned as a valid equilibrium unless gates pass.
6. Remove any test that allows `no_convergence` as an accepted state.

---

### Repair C — Replace accepted neutral LLE hand-coded minimizer

#### C1. Current route to replace

Current accepted neutral LLE route:

```text
solve_lle_attempt(...)
→ minimize_lle_residual_variables(...)
→ lle_two_phase_result(...)
→ nonlinear_solver = native_derivative_free_nelder_mead
```

This must no longer be an accepted production route.

#### C2. Required production replacement

Use the same shared native residual solver abstraction as electrolyte and reactive LLE.

Required residual blocks for neutral LLE:

```math
R_i^{phase} =
\left(\ln x_i^\alpha + \ln\phi_i^\alpha\right)
-
\left(\ln x_i^\beta + \ln\phi_i^\beta\right)
```

```math
R_i^{mat} =
(1-\beta)x_i^\alpha+\beta x_i^\beta-z_i
```

Use transformed variables for:

```text
phase fraction
phase compositions
density or density closure if required by chosen formulation
```

Ceres must receive a real Jacobian from analytic/CppAD/implicit derivatives.

#### C3. Delete or quarantine old helpers

The repair PR must do one of the following:

Preferred:

```text
delete minimize_lle_residual_variables(...)
```

Allowed only if clearly non-production:

```text
rename to legacy_seed_candidate_for_tests_only(...)
move out of production solve path
add tests proving it cannot produce accepted results
```

Do not leave the same function name on the accepted route.

#### C4. Required diagnostics

Accepted neutral LLE diagnostics must report:

```text
solver_backend = ceres
solver_method = ceres_trust_region_residual_solve
jacobian_backend = analytic or cppad or analytic_implicit or cppad_implicit
jacobian_available = true
phase_equilibrium_residual_norm
material_balance_norm
phase_distance
acceptance_gate = residual_and_physical_gates
```

Accepted neutral LLE diagnostics must not report:

```text
native_derivative_free_nelder_mead
not_applicable derivative backend
not_required derivative status
```

#### C5. Required tests

Add or update:

```text
tests/native/equilibrium/test_neutral_lle_ceres_solver.py
tests/native/equilibrium/test_neutral_lle_residual_jacobian.py
tests/equilibrium/core/test_neutral_lle_production_solver.py
```

Tests must fail if accepted neutral LLE still reports the old solver label.

---

### Repair D — Honest diagnostics for failed solves

#### D1. Separate attempted solver from accepted solver

For every equilibrium result, distinguish:

```text
solver_attempted
solver_attempt_result
accepted_solver_backend
accepted_solver_method
accepted_derivative_backend
solution_accepted
```

Use existing result schema names if better names already exist, but preserve the distinction.

#### D2. Failure diagnostics rules

If an equilibrium solve fails:

```text
solution_accepted = false
split_detected = false unless explicitly documented as unstable/feed-only diagnostic
ceres_accepted_solve = false
jacobian_available_for_accepted_state = false
derivative_available_for_accepted_state = false
accepted_solver_backend = none
accepted_derivative_backend = none
```

It may report:

```text
solver_attempted = ceres
attempted_jacobian_backend = cppad_implicit
last_residual_norm = ...
failure_reason = ...
```

Do not make failed results look accepted.

#### D3. Required tests

Add tests that intentionally fail solves and assert diagnostics are honest:

```text
tests/native/equilibrium/test_failed_solve_diagnostics_are_not_accepted.py
tests/equilibrium/electrolyte/test_failed_electrolyte_lle_diagnostics.py
tests/equilibrium/reactive/test_failed_reactive_lle_diagnostics.py
```

---

### Repair E — Real benchmark proof

#### E1. Replace weak salting-out benchmark

The #116 benchmark must use a sourced salting-out electrolyte LLE fixture. Preferred system:

```text
water + ethanol + isobutyl alcohol + NaCl
293.15 K and/or 313.15 K
5 wt% and/or 10 wt% NaCl
```

Data source:

```text
Khudaida et al. 2026, Journal of Chemical & Engineering Data 71, 708–717
```

Required fixture file:

```text
tests/fixtures/literature/khudaida_2026_salting_out_lle.json
```

or equivalent repo-standard fixture path.

Required fixture metadata:

```json
{
  "source": "Khudaida et al. 2026",
  "doi_or_citation": "...",
  "temperature_K": ...,
  "pressure_Pa": ...,
  "salt_wt_percent": ...,
  "species": ["H2O", "Ethanol", "Isobutanol", "Na+", "Cl-"],
  "phase_data_basis": "...",
  "extraction_notes": "..."
}
```

Do not invent values. If the table values are not yet extracted, extract them from the paper/SI or leave the benchmark task open.

#### E2. Replace weak neutral reactive LLE benchmark

The #117 neutral reactive LLE benchmark must use an Ascani 2023-style esterification system.

Preferred systems:

```text
acetic acid + 1-pentanol ⇌ ester + water
acetic acid + 1-hexanol ⇌ ester + water
```

The benchmark must include:

```text
at least four species
one reversible esterification reaction
two liquid phases
reaction residual checks
phase-equilibrium residual checks
element/material balance checks
source-backed compositions or tolerances
```

Do not use two-component artificial transformations like:

```text
Methanol ⇌ Cyclohexane
```

as literature proof.

#### E3. Reactive electrolyte benchmark

A reactive electrolyte LLE benchmark must include a reaction involving charged species or changing charged/neutral distribution. A neutral solvent conversion with spectator ions is not enough.

Acceptable generic form for a real coupled test:

```text
Li+_aq + HR_org ⇌ RLi_org + H+_aq
Na+_aq + HR_org ⇌ RNa_org + H+_aq
```

This may be a source-backed lithium-related fixture, but the package API must remain generic. Do not expose extraction efficiency, distribution coefficient, or selectivity in the package API.

If a source-backed reactive electrolyte fixture cannot be completed in this PR, then:

```text
do not claim reactive electrolyte literature benchmark complete
keep a dedicated follow-up issue open
state clearly what was implemented and what remains unvalidated
```

Do not close with a synthetic-only proof.

---

### Repair F — Shared native residual solver architecture

#### F1. Add reusable solver abstraction

Add or finish:

```text
src/epcsaft/native/equilibrium/residual_solver.h
src/epcsaft/native/equilibrium/residual_solver.cpp
```

Expected interface:

```cpp
struct NativeResidualSolveOptions {
    int max_iterations;
    double residual_tolerance;
    double step_tolerance;
    double function_tolerance;
};

struct NativeResidualEvaluation {
    bool success;
    std::vector<double> residual;
    std::vector<double> jacobian_row_major;
    int rows;
    int cols;
    std::vector<std::string> residual_names;
    std::vector<std::string> variable_names;
    std::map<std::string, double> diagnostics_double;
    std::map<std::string, int> diagnostics_int;
    std::map<std::string, bool> diagnostics_bool;
    std::map<std::string, std::string> diagnostics_string;
};

class NativeResidualProblem {
public:
    virtual ~NativeResidualProblem() = default;
    virtual int variable_count() const = 0;
    virtual int residual_count() const = 0;
    virtual NativeResidualEvaluation evaluate(
        const std::vector<double>& variables,
        bool need_jacobian
    ) const = 0;
};

struct NativeResidualSolveResult {
    bool solver_usable;
    bool accepted;
    std::string solver_backend;
    std::string solver_method;
    std::string termination_type;
    std::string rejection_reason;
    int iterations;
    double initial_cost;
    double final_cost;
    std::vector<double> variables;
    NativeResidualEvaluation final_evaluation;
};

NativeResidualSolveResult solve_native_residual_problem(
    const NativeResidualProblem& problem,
    const std::vector<double>& initial_variables,
    const NativeResidualSolveOptions& options
);
```

Exact names may differ, but the responsibilities must be separated.

#### F2. All three production routes must use it

Routes that must use the shared abstraction:

```text
neutral LLE
electrolyte LLE
reactive LLE / reactive electrolyte LLE
```

Do not maintain three separate one-off Ceres acceptance implementations.

#### F3. Centralize diagnostics

The shared solver must populate consistent solver diagnostics, including:

```text
solver_backend
solver_method
termination_type
solver_usable
accepted
rejection_reason
initial_cost
final_cost
iteration_count
residual_norm_linf
residual_norm_l2
jacobian_available
jacobian_backend
```

Problem-specific routes add physical residual norms.

---

## 5. Mandatory audit commands

The repair agent must run these before coding and write results to:

```text
docs/goals/pr126-repair-audit/notes/current_state_audit.md
```

Commands:

```powershell
git fetch origin --prune
git merge --ff-only origin/main
git status --short
```

Search old routes:

```powershell
rg -n "minimize_lle_residual_variables|native_derivative_free_nelder_mead|polish_formula_tpd_variables|nelder_mead|simplex" src tests docs
```

Search reactive standard-state handling:

```powershell
rg -n "reaction_standard_states|native_standard_state_code|standard_state|ReactionConstantConvention|reaction_residuals" src tests docs
```

Search Ceres acceptance:

```powershell
rg -n "ceres_accepted_solve|IsSolutionUsable|NO_CONVERGENCE|ceres_termination_type|summary.termination_type" src tests docs
```

Search benchmark shortcuts:

```powershell
rg -n "model-consistent|Methanol|Cyclohexane|water_to_butanol|H2O.*Butanol|repo-contained model-consistent" tests data docs
```

The audit file must include:

```text
where each old route still exists
whether it is production, seed-only, test-only, or docs-only
where each repair will be made
which tests currently permit the bad behavior
```

Do not edit source before this audit file exists.

---

## 6. Required tests and their purpose

The repair PR must add or update tests so the previous mistakes are impossible to repeat.

### 6.1 Tests for reaction standard states

Required files:

```text
tests/native/equilibrium/test_reactive_phase_equilibrium_standard_states.py
tests/equilibrium/reactive/test_reactive_lle_standard_state_conventions.py
```

Required assertions:

```python
assert residual_for_mole_fraction_activity != residual_for_ideal_mole_fraction
assert residual_for_concentration != residual_for_mole_fraction_activity
assert native_residual == independently_calculated_expected_residual
assert diagnostics["reaction_residual_standard_state_applied"] is True
```

### 6.2 Tests for Ceres acceptance

Required files:

```text
tests/native/equilibrium/test_reactive_phase_equilibrium_ceres_acceptance.py
tests/native/equilibrium/test_electrolyte_lle_ceres_acceptance.py
tests/equilibrium/reactive/test_failed_reactive_lle_diagnostics.py
```

Required assertions:

```python
assert diagnostics["ceres_accepted_solve"] is True only when solver usable and physical gates pass
assert diagnostics["ceres_termination_type"] == "convergence" for accepted benchmark solves
assert diagnostics["ceres_termination_type"] != "no_convergence" for accepted benchmark solves
assert failed_result_or_exception.diagnostics["solution_accepted"] is False
```

### 6.3 Tests for neutral LLE solver replacement

Required files:

```text
tests/native/equilibrium/test_neutral_lle_ceres_solver.py
tests/native/equilibrium/test_neutral_lle_residual_jacobian.py
tests/equilibrium/core/test_neutral_lle_production_solver.py
```

Required assertions:

```python
assert diagnostics["solver_backend"] == "ceres"
assert diagnostics["jacobian_available"] is True
assert "native_derivative_free_nelder_mead" not in json.dumps(result.to_dict())
assert "not_required" not in json.dumps(result.to_dict())
```

### 6.4 Tests for benchmark correctness

Required files or equivalent:

```text
tests/equilibrium/electrolyte/test_khudaida_2026_salting_out_lle.py
tests/equilibrium/reactive/test_ascani_2023_esterification_reactive_lle.py
```

Do not reuse the current toy fixtures as benchmark proof. They may remain as smoke tests if clearly named as smoke tests.

### 6.5 Tests for diagnostic honesty

Required assertions:

```python
assert failed_diagnostics["accepted_solver_backend"] in {"", None, "none"}
assert failed_diagnostics["solution_accepted"] is False
assert failed_diagnostics["jacobian_available_for_accepted_state"] is False
assert "accepted" not in failed_diagnostics["acceptance_gate"].lower()
```

Use actual field names from the implemented diagnostic schema.

---

## 7. Required validation commands

The PR must run and report these exact commands unless the repo has changed commands after this document. If commands changed, update the audit note with the new repo-standard command and explain why.

```powershell
uv run python scripts/dev/build_epcsaft.py --clean --enable-cppad
uv run python scripts/dev/build_epcsaft.py --clean --enable-ceres --enable-cppad
uv run python run_pytest.py tests/native/equilibrium -q
uv run python run_pytest.py tests/equilibrium/electrolyte -q
uv run python run_pytest.py tests/equilibrium/reactive -q
uv run python run_pytest.py tests/equilibrium/core -q
uv run python run_pytest.py tests/api/equilibrium -q
uv run python run_pytest.py tests/api/reactive -q
uv run python scripts/dev/validate_project.py quick
uv run python scripts/dev/validate_project.py docs
git diff --check
```

Also run final route guards:

```powershell
rg -n "native_derivative_free_nelder_mead|not_required|model-consistent reactive|Methanol.*Cyclohexane|water_to_butanol|no_convergence" src tests docs
```

Any remaining match must be one of:

```text
historical audit text
test proving bad route is absent
explicitly named smoke test not used as benchmark proof
compatibility helper not on accepted production route
```

No remaining match may appear in accepted production diagnostics or benchmark completion claims.

---

## 8. Required PR body

The repair PR body must use this structure.

```markdown
## Summary

- Fixes PR #126 incomplete completion.
- Fixes reactive standard-state residual handling.
- Fixes Ceres acceptance gates.
- Replaces accepted neutral LLE hand-coded minimizer route.
- Makes failed-solve diagnostics honest.
- Replaces weak benchmark smokes with source-backed benchmark proof or leaves explicit open follow-up.

## Fixed P1 review items

- [ ] P1: reaction standard-state codes applied in native residual equations.
- [ ] P1: reactive Ceres solve accepts only usable Ceres solutions and passed physical gates.

## Old route removal

- [ ] `minimize_lle_residual_variables` no longer produces accepted production LLE.
- [ ] `native_derivative_free_nelder_mead` absent from accepted production diagnostics.
- [ ] `polish_formula_tpd_variables` removed or seed-only/test-only and cannot accept equilibrium.
- [ ] failed solves no longer report accepted derivative status.

## Benchmarks

- [ ] Khudaida-style salting-out fixture source and tolerances listed.
- [ ] Ascani 2023-style esterification reactive LLE fixture source and tolerances listed.
- [ ] Reactive electrolyte fixture status stated honestly.

## Validation

Paste exact command output summary here.

## Remaining work

If any item remains, do not write "Closes". Leave the relevant issue open or create a follow-up issue before requesting merge.
```

Do not include `Closes #116` or `Closes #117` in this repair PR unless every item in this audit is complete.

---

## 9. Completion checklist

The repair PR is complete only if every item below is true.

### Reactive residual correctness

- [ ] Native reaction residuals use `reaction_standard_states`.
- [ ] Tests cover at least mole-fraction activity, ideal mole fraction, and concentration conventions.
- [ ] Non-activity conventions no longer silently behave like mole-fraction activity.

### Ceres acceptance correctness

- [ ] Reactive route checks `summary.IsSolutionUsable()`.
- [ ] Accepted production routes reject no-convergence.
- [ ] Final residuals are re-evaluated after Ceres.
- [ ] Physical gates must pass before phases are returned as accepted.
- [ ] Failure cases do not return accepted two-phase outputs.

### LLE solver replacement

- [ ] Accepted neutral LLE no longer uses the hand-coded simplex minimizer.
- [ ] Accepted electrolyte LLE uses shared native residual solver.
- [ ] Accepted reactive LLE uses shared native residual solver.
- [ ] No accepted production result reports old solver labels.

### Derivative and diagnostics honesty

- [ ] Accepted result has real Jacobian provenance.
- [ ] Failed result does not claim accepted derivative availability.
- [ ] Attempted solver fields are separate from accepted solver fields.
- [ ] No diagnostic label is used as proof without a residual/Jacobian implementation.

### Benchmarks

- [ ] Salting-out benchmark uses source-backed water + ethanol + isobutanol + NaCl data or an explicitly approved equivalent.
- [ ] Reactive LLE benchmark uses source-backed esterification-style data or remains open.
- [ ] Reactive electrolyte test involves charged species in the reaction or is explicitly not claimed as literature benchmark completion.

### API and scope

- [ ] Public API remains generic.
- [ ] No MEA-specific, lithium-specific, extraction-specific, selectivity, distribution-coefficient, or absorber-column public APIs are added.

### Validation

- [ ] All required build/test/doc commands pass.
- [ ] Final route guards pass or remaining hits are documented and harmless.
- [ ] PR body lists exact validation results.

---

## 10. Stop rules

Stop and do not request merge if any of these happen:

```text
reaction standard-state convention remains ignored
Ceres no-convergence can still be accepted
failed result still reports accepted derivative status
neutral LLE still accepts native_derivative_free_nelder_mead
benchmark proof is still only Methanol ⇌ Cyclohexane or H2O ⇌ Butanol
manual numeric perturbation is used for production derivatives
source-backed benchmark data cannot be curated
```

In a stopped state, create:

```text
docs/goals/pr126-repair-audit/notes/stopped_state.md
```

and include:

```text
exact blocker
files inspected
partial changes made
tests added
what remains
why the PR must stay draft or open
```

Do not write a completion audit if a stop rule is triggered.

---

## 11. Suggested one-shot prompt for the repair agent

Copy and paste this prompt into the Codex agent after attaching this Markdown file.

```text
Thread name: ePC-SAFT PR126 Repair Audit and Completion

You are the repair agent for tannerpolley/ePC-SAFT after PR #126 merged incomplete work for issues #116 and #117.

Read the attached file `pr126_repair_audit_and_completion_contract.md` completely before planning or editing.

Your job is not to defend PR #126. Your job is to audit current main, fix the incomplete implementation, and prove the fixes with strict tests.

Start on a fresh branch from current origin/main.

Before editing source, create:
docs/goals/pr126-repair-audit/notes/current_state_audit.md

Run and record the audit commands from the attached document.

Then implement all required repairs:
1. Apply reaction standard-state conventions in native reactive residuals.
2. Gate reactive Ceres acceptance on solver usability and physical residual gates.
3. Replace accepted neutral LLE hand-coded minimizer route with the shared derivative-backed residual solver.
4. Make failed-solve diagnostics honest.
5. Replace weak benchmark smokes with source-backed benchmark proof or leave explicit open follow-up without claiming completion.
6. Add a shared native residual-solver abstraction and route neutral LLE, electrolyte LLE, and reactive LLE through it.
7. Add the tests listed in the audit document.

Do not ask me whether to reduce scope. Do not close with docs-only, diagnostics-only, smoke-only, or "future work" completion.

If any stop rule in the attached document is triggered, stop, write the stopped-state note, and keep the PR draft/open.

The PR may be ready only when every item in the completion checklist passes.
```
