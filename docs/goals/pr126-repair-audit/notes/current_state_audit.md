# PR126 Repair Current-State Audit

Date: 2026-05-15

Branch: `pr126-review-remediation`

Scope: Audit required by `docs/plans/pr126_repair_audit_and_completion_contract.md` before treating the repair branch as admissible for completion. The contract document was supplied after initial P1 repair edits had already started; this note records the current branch state and keeps completion claims blocked until every stop rule is cleared.

## Commands Run

```powershell
git fetch origin --prune
git merge --ff-only origin/main
git status --short --branch
```

Result: `origin/main` was already up to date. The working tree was dirty with active repair edits in native reactive equilibrium code, native reactive tests, and the reopened GoalBuddy state.

```text
## pr126-review-remediation...origin/main
 M docs/goals/native-lle-reactive-production-solvers/state.yaml
 M src/epcsaft/native/equilibrium/reactive_phase_equilibrium_problem.cpp
 M tests/native/equilibrium/test_reactive_phase_equilibrium_ceres_solver.py
 M tests/native/equilibrium/test_reactive_phase_equilibrium_residual_jacobian.py
 M tests/native/equilibrium/test_reactive_phase_equilibrium_residual_surface.py
?? docs/goals/native-lle-reactive-production-solvers/notes/T028-pr126-post-merge-review.md
?? docs/handoffs/
?? docs/plans/pr126_repair_audit_and_completion_contract.md
```

## Old Solver Route Search

Command:

```powershell
rg -n "minimize_lle_residual_variables|native_derivative_free_nelder_mead|polish_formula_tpd_variables|nelder_mead|simplex" src tests docs -S
```

Findings:

- Production: `src/epcsaft/native/epcsaft_equilibrium.cpp` still defines `minimize_lle_residual_variables(...)`, calls it from the neutral LLE attempt path, and reports `nonlinear_solver = native_derivative_free_nelder_mead` on accepted neutral LLE results.
- Seed support: `src/epcsaft/native/epcsaft_equilibrium.cpp` still defines and calls `polish_formula_tpd_variables(...)` during formula-TPD seed polishing. This must either become clearly seed-only or be replaced.
- Tests permitting old behavior: `tests/equilibrium/core/test_lle.py` and `tests/native/cppad/test_cppad_lle_derivatives.py` still assert or allow `native_derivative_free_nelder_mead`.
- Docs and goal notes contain historical mentions. Those may remain only if clearly historical and not completion evidence.

Repair locations:

- `src/epcsaft/native/epcsaft_equilibrium.cpp`
- new shared residual solver files under `src/epcsaft/native/equilibrium/`
- neutral LLE tests under `tests/native/equilibrium/` and `tests/equilibrium/core/`

## Reactive Standard-State Search

Command:

```powershell
rg -n "reaction_standard_states|native_standard_state_code|standard_state|ReactionConstantConvention|reaction_residuals" src tests docs -S
```

Findings:

- Production: the original PR126 route accepted `reaction_standard_states`, but the native reactive reaction residual builder used only activity terms. The branch now has in-progress edits in `src/epcsaft/native/equilibrium/reactive_phase_equilibrium_problem.cpp`; they must be verified against source-backed, non-synthetic tests before acceptance.
- Source-of-truth mapping must be checked in `src/epcsaft/reactive.py`, `src/epcsaft/reactive_speciation.py`, `src/epcsaft/equilibrium.py`, `src/epcsaft/bindings.cpp`, and `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp`.
- Tests currently touched by this branch include synthetic fixtures and cannot satisfy the user hard requirement. They must be removed or replaced with source-backed tests.

Repair locations:

- `src/epcsaft/native/equilibrium/reactive_phase_equilibrium_problem.cpp`
- `src/epcsaft/equilibrium.py` diagnostics mapping if native diagnostics gain new fields
- source-backed standard-state tests only

## Ceres Acceptance Search

Command:

```powershell
rg -n "ceres_accepted_solve|IsSolutionUsable|NO_CONVERGENCE|ceres_termination_type|summary\.termination_type" src tests docs -S
```

Findings:

- Production reactive route in `src/epcsaft/native/equilibrium/reactive_phase_equilibrium_problem.cpp` had unconditional `ceres_accepted_solve=true`; branch edits added an initial acceptance gate, but it still requires review against the contract because `USER_SUCCESS` must not be accepted without a callback and physical gates must be explicit.
- Production electrolyte route in `src/epcsaft/native/epcsaft_equilibrium.cpp` uses Ceres and `summary.IsSolutionUsable()`, but failed-solve diagnostics need to separate attempted and accepted solver state.
- Test weakness: `tests/native/equilibrium/test_electrolyte_lle_ceres_solver.py` allows `ceres_termination_type` to be either `convergence` or `no_convergence`, which is invalid for an accepted production benchmark.

Repair locations:

- shared acceptance helper under `src/epcsaft/native/equilibrium/`
- reactive route acceptance in `src/epcsaft/native/equilibrium/reactive_phase_equilibrium_problem.cpp`
- electrolyte failure diagnostics in `src/epcsaft/native/epcsaft_equilibrium.cpp`
- Ceres acceptance tests under `tests/native/equilibrium/`

## Benchmark Shortcut Search

Command:

```powershell
rg -n "model-consistent|Methanol|Cyclohexane|water_to_butanol|H2O.*Butanol|repo-contained model-consistent" tests data docs -S
```

Findings:

- `tests/equilibrium/reactive/test_reactive_lle_coupled_solver.py` uses a Methanol/Cyclohexane model-consistent fixture. It is smoke coverage only and cannot be benchmark proof.
- `tests/equilibrium/reactive/test_reactive_lle.py` uses Methanol/Cyclohexane staged reactive fixtures. They cannot be cited as literature proof.
- `tests/equilibrium/reactive/test_reactive_electrolyte_lle_coupled_solver.py` uses a neutral `H2O -> Butanol` reaction with spectator ions. This does not satisfy the reactive electrolyte coupling requirement.
- In-progress branch tests under `tests/native/equilibrium/test_reactive_phase_equilibrium_residual_surface.py`, `tests/native/equilibrium/test_reactive_phase_equilibrium_residual_jacobian.py`, and `tests/native/equilibrium/test_reactive_phase_equilibrium_ceres_solver.py` also use synthetic fixtures and must be removed or replaced before completion.

Repair locations:

- source-backed Khudaida 2026 electrolyte salting-out fixture/tests
- source-backed Ascani 2023 esterification reactive LLE fixture/tests, if enough data can be curated
- explicit stopped-state note if source-backed benchmark data cannot be curated in this PR

## Current Stop-Rule Status

Completion is blocked while any of these remain true:

- neutral LLE accepted production route still reports `native_derivative_free_nelder_mead`
- accepted electrolyte tests still permit `no_convergence`
- failed electrolyte diagnostics still conflate attempted solver state with accepted derivative availability
- benchmark proof is still based only on Methanol/Cyclohexane or H2O/Butanol shortcuts
- branch-local synthetic tests remain as proof of the new behavior
- source-backed reactive benchmark data is not curated or explicitly left open without a completion claim

This audit is not a completion receipt.
