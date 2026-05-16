# T002 Scout Route Map: Issue #116

Date: 2026-05-15

## Current Public Route

- `ElectrolyteLLEProblem.solve` delegates to `mixture.electrolyte_lle_tp(...)`.
- `electrolyte_lle_flash_native` normalizes the feed, builds the electrolyte basis, and calls `_call_native_equilibrium(kind="electrolyte_lle", ...)`.
- `src/epcsaft/electrolyte_lle.py` is absent, so any issue text naming that file maps to `src/epcsaft/equilibrium.py` plus native C++.

## Current Native Route

- `build_electrolyte_basis_native` constructs neutral, cation, and anion indices plus salt pairs and formula feed.
- `electrolyte_tpd_global_search` provides deterministic stability and seed candidates.
- `solve_predictive_electrolyte_attempt` currently uses `nelder_mead_variables` before `newton_step`.
- Successful electrolyte diagnostics can include `backend = electrolyte_lle`, `problem_kind = electrolyte_lle_flash`, `stability_analysis = electrolyte_tpd`, `variable_model = ascani_transformed_salt_pairs`, `solver_method = native_transformed_newton`, `solver_language = c++`, and `native_entrypoint = _solve_equilibrium_native`.

## Hand-Coded Optimizer Evidence

- `minimize_lle_residual_variables` implements the old neutral LLE derivative-free minimizer path.
- `polish_formula_tpd_variables` implements derivative-free TPD polishing.
- `nelder_mead_variables` remains part of the predictive electrolyte attempt.
- `solve_lle_attempt` still calls the old neutral LLE minimizer.

## Residual And Jacobian Gap

- `newton_step` throws an unsupported-derivative error for electrolyte LLE residual sensitivities.
- `_evaluate_electrolyte_lle_residual_native` throws the same unavailable sensitivity message.
- Contract tests currently assert unavailable electrolyte residual sensitivities, so they document a gap rather than proving a production Jacobian.

## Tests That Can Pass On The Old Route

- `tests/equilibrium/electrolyte/test_electrolyte_lle_solver_contracts.py` accepts current `native_transformed_newton` diagnostics and split detection.
- `tests/equilibrium/electrolyte/test_electrolyte_lle_smokes.py` accepts current predictive diagnostics for explicit initial phases.
- `tests/native/contracts/test_equilibrium_native_contracts.py` includes tests that currently expect unavailable residual sensitivities.

## Fixture Candidates

- Ascani-style distributed ions: existing `2022_Ascani` water/butanol/salt fixtures and `data/multiphase/ascani_case2_model_comparison.csv`.
- Salting-out benchmark: existing Khudaida 2026 data and parameters under `data/reference/equilibrium_benchmarks/electrolyte_lle/khudaida_2026` and `data/reference/epcsaft_parameters/2026_Khudaida`.
- Optional hard-case stress: Hubach 2024 LiCl/TBP ionic-liquid fixture, currently opt-in.

## Worker Boundary Recommendation

Start issue #116 source work with native basis and transformed-variable ownership first, then introduce residual/Jacobian/Ceres solve proof, then wire the public API and benchmark tests. Do not relabel existing predictive diagnostics as production completion.
