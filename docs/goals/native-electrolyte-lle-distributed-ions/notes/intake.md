# Issue #116 Stage 0 Intake

Date: 2026-05-15

Issue: https://github.com/tannerpolley/ePC-SAFT/issues/116

## Dependency Gate

- Branch: `#116-117`
- HEAD before source edits: `4f1efb046b4614e9a941ded5371349d9246924bd`
- `origin/main` after `git fetch origin --prune`: `f3ce3cc332bc06e8504ece122e478c8396bc27b8`
- Fast-forward gate: `git merge --ff-only origin/main` reported `Already up to date.`
- Worktree before source edits: clean by `git status --short`.
- Required dependency PRs are merged and ancestors of this branch:
  - PR #123, merge commit `5391b4011f46f6817081a81faec0f3c5f407580d`
  - PR #124, merge commit `b8cc04b095c0cd9d53d7de08c9315d4040e691ff`
  - PR #125, merge commit `acf6fc74c90734f2cbcab4f79dfa73a801297de6`

Source-edit gate: open for issue #116 only after this intake is recorded. Issue #117 remains dependency-gated behind #116 production electrolyte LLE completion unless the same branch completes #116 first.

## Current Route Map

- `src/epcsaft/electrolyte_lle.py` does not exist in the current tree. The public route is in `src/epcsaft/equilibrium.py`.
- `ElectrolyteLLEProblem.solve` delegates to `mixture.electrolyte_lle_tp(...)`.
- `electrolyte_lle_flash_native` builds an electrolyte formula basis and calls `_call_native_equilibrium(kind="electrolyte_lle", ...)`.
- The native implementation lives in `src/epcsaft/native/epcsaft_equilibrium.cpp`, with basis construction, TPD seeding, a predictive solve loop, and result diagnostics in the same large file.
- `build_electrolyte_basis_native` currently constructs neutral indices, cation indices, anion indices, salt pairs, formula feed, basis rank, and variable model `ascani_transformed_salt_pairs`.
- Current accepted electrolyte diagnostics can report `solver_method = native_transformed_newton` and `acceptance_gate = predictive_nonlinear_solve`.

## Old Solver Paths That Must Not Remain Production

- Neutral LLE still uses `minimize_lle_residual_variables`, which is a hand-coded derivative-free minimizer path.
- Electrolyte TPD polishing uses `polish_formula_tpd_variables`, another hand-coded derivative-free polishing route.
- `solve_lle_attempt` still calls the hand-coded neutral minimizer.
- `solve_predictive_electrolyte_attempt` still calls `nelder_mead_variables` and then `newton_step`.
- `newton_step` currently throws an unsupported-derivative error for electrolyte LLE residual sensitivities.
- `_evaluate_electrolyte_lle_residual_native` currently throws an unsupported-derivative error for electrolyte LLE residual sensitivities.
- Neutral LLE diagnostics can still report `nonlinear_solver = native_derivative_free_nelder_mead`, `jacobian_available = false`, and derivative status values that are not acceptable for the issue #116 production route.

## Required #116 Production Direction

Issue #116 requires replacing the accepted production electrolyte LLE solve with a native Ceres trust-region residual solve on transformed feasible variables. The production route must own:

- Explicit distributed-ion phase compositions in public species space.
- A charge-neutral independent coordinate basis only as coordinates, not as hidden public composition state.
- Residual blocks for phase equilibrium, ionic combination equilibrium, material accounting, charge accounting where not enforced by transform, and closure.
- A real Jacobian path with diagnostics proving availability and provenance.
- Ceres solver diagnostics from the accepted solve state.
- TPD and g-hat limited to stability, seed, phase-count, and acceptance support.

## Chosen Fixture Targets

- Ascani-style distributed-ion fixture: existing `2022_Ascani` water/butanol/Na+/K+/Cl- and mixed Na+/K+/Cl- routes used by `tests/equilibrium/electrolyte/test_electrolyte_lle_solver_contracts.py`, including `data/multiphase/ascani_case2_model_comparison.csv`.
- Salting-out benchmark fixture: existing Khudaida 2026 electrolyte LLE benchmark data in `data/reference/equilibrium_benchmarks/electrolyte_lle/khudaida_2026` with `data/reference/epcsaft_parameters/2026_Khudaida`.
- Optional hard-case proof: Hubach 2024 LiCl/TBP ionic-liquid fixture under `data/reference/equilibrium_benchmarks/electrolyte_lle/hubach_2024`, currently opt-in and not a substitute for a routine production benchmark.

## Tests That Currently Would Not Prove Completion

- Tests asserting only `solver_method == native_transformed_newton` or a split result can still pass on the old predictive path.
- Tests expecting unavailable residual sensitivities currently document the gap rather than proving production behavior.
- Opt-in hard-case tests are useful stress evidence but cannot be the only production definition-of-done proof.

## Source-Edit Entry Condition

The next implementation package must remove accepted production dependence on the hand-coded minimizer path, add transformed variable and basis diagnostics, then introduce Ceres residual/Jacobian proof before any issue #116 completion claim.
