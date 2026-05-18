# T001 Scout Receipt

Result: done

Branch/status: `codex/rezaee-reactive-electrolyte-lle`, clean short status.

## Summary

The live worktree supports the five-candidate architecture plan. Evidence shows typed `EquilibriumProblem` objects, native Ipopt-gated production equilibrium routes, `TargetDataset` and regression APIs, parameter-family schema coverage, and runtime capability contracts.

Main stale-risk: binary regression exposes `l_ij` and `k_hb_ij` as parameter families, but the native Ceres binary production path currently supports only constant `k_ij`. Treat `l_ij` and `k_hb_ij` optimizer support as pending, not complete.

## Evidence Map

- Equilibrium Problem: `src/epcsaft/equilibrium.py` defines `EquilibriumProblem` and concrete problem objects; `src/epcsaft/epcsaft.py` exposes `mixture.solve_equilibrium(problem)`.
- Production Solver Path: `src/epcsaft/equilibrium.py` centralizes native Ipopt-required failures and calls native route bindings; `tests/native/equilibrium/test_route_builders.py` exercises native route builders.
- Target Dataset and Regression Problem: `src/epcsaft/regression.py` defines `TargetRow` and `TargetDataset` and exposes public fit helpers; `src/epcsaft/reactive_regression.py` defines reactive batch/context objects.
- Parameter Family: `src/epcsaft/parameter_schema.py` covers pure, association, Born, ion, and binary `k_ij`/`l_ij`/`k_hb_ij` records; `src/epcsaft/parameters.py` loads and populates those matrices.
- Capability Contract: `src/epcsaft/runtime.py` is the capability surface; tests under `tests/api/runtime/` check native backend claims and reactive diagnostic wording.

## Recommended Order

1. Production Solver Path
2. Equilibrium Problem
3. Target Dataset and Regression Problem
4. Parameter Family
5. Capability Contract

Capability updates should trail implementation and tests so they do not overclaim.

## Recommended First Worker Package

Name: Production Solver Path capability-tightening slice

Why: it is the dependency base for the other four candidates and minimizes overclaiming risk before broader API/schema work.

Allowed files:

- `src/epcsaft/equilibrium.py`
- `src/epcsaft/equilibrium_core/**`
- `src/epcsaft/runtime.py`
- `tests/equilibrium/**`
- `tests/native/equilibrium/**`
- `tests/api/runtime/test_runtime_exports_and_metadata.py`

Verify:

- `uv run python run_pytest.py tests/equilibrium/core tests/equilibrium/electrolyte tests/equilibrium/reactive tests/api/runtime/test_runtime_exports_and_metadata.py -q`
- `uv run python run_pytest.py tests/native/equilibrium/test_route_builders.py -q`

Stop if:

- A change would require rebuilding or deleting `_core` in this Worker package.
- Native `_core` lacks the route binding needed for the intended Production Solver Path contract.
- Any capability string would need to claim production support without a passing public-interface and native-route test.
- The work needs changes outside allowed files, especially regression or parameter schema files.
