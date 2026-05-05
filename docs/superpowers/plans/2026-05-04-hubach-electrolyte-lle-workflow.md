# General Electrolyte Phase-Equilibrium Workflow Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Improve the general electrolyte phase-equilibrium algorithm and workflow, using Hubach 2024 as the first hard fixed-species case study rather than as a solver special case.

**Architecture:** Preserve `ePCSAFTMixture.equilibrium(kind="electrolyte_lle")` as the public entrypoint and keep the native C++ backend as the main solver. Expand the charge-constrained formula basis beyond 1:1 salts, record seed-attempt diagnostics, try broader seed families, add continuation support, add Hubach/divalent-ion regression coverage, and preserve a diagnostic-only legacy candidate when strict native acceptance still fails.

**Tech Stack:** Python 3, uv, pytest through `run_pytest.py`, pybind11, native C++ equilibrium code, existing ePC-SAFT dataset and benchmark layouts.

---

## Summary

The current Hubach 2024 failure is a workflow and algorithm stress test: electrolyte TPD reports an unstable feed, but the predictive nonlinear solve satisfies residuals at a collapsed one-phase-like point. The fix should not loosen the distinct-phase acceptance gate. Instead, the implementation should improve the general electrolyte LLE path in this order:

1. Add richer diagnostics so failed cases explain what seeds were tried and why they were rejected.
2. Add robust seed families and continue after collapsed candidates when TPD says the feed is unstable.
3. Add continuation helpers for O/A and composition curves.
4. Generalize the charge-constrained electrolyte formula basis beyond 1:1 salts.
5. Add divalent-ion basis support and regression coverage for Mg-containing systems.
6. Add a diagnostic-only Python legacy-candidate fallback for native failures where strict acceptance remains invalid.

Hubach remains a fixed-species phase-equilibrium case study. Chemical/speciation equilibrium, Yu-style complexation, and extraction-efficiency matching are explicitly outside this implementation phase.

## Key Implementation Changes

- Update `src/epcsaft/equilibrium_core/electrolyte_basis.py` and the native basis in `src/epcsaft/native/epcsaft_equilibrium.cpp` so salt pairs carry stoichiometric cation/anion coefficients. The basis must support 1:1, 2:1, 1:2, and mixed cation/anion systems while preserving charge neutrality by construction.
- Keep `initial_phases={"aq": ..., "org": ..., "phase_fraction": ...}` unchanged. Improve its validation diagnostics and make it the first seed attempted when supplied.
- Add seed-attempt diagnostics to native result payloads: `seed_attempts`, `seed_attempt_count`, `best_noncollapsed_candidate`, and `unstable_feed_collapsed_all_candidates`.
- Broaden native seed generation: TPD trial, Gibbs-refined TPD, formula component-rich seeds, water/organic endpoint seeds, salt aqueous-rich seeds, salt partially organic seeds, and beta sweeps.
- Add `initial_phases_from_result(result)` and `mixture.equilibrium_curve(points, kind="electrolyte_lle", ...)` for continuation-first curve solving.
- Add the Hubach 2024 dataset fixture and row 0 regression as a case study, not a hard-coded exception.
- Align the `2024_Hubach` option surface with Lithium's canonical `huback_2024` settings: combined relative-permittivity rule, analytical derivative modes, Born enabled, fitted `d_Born_mode`, no SSM/DS, Born bulk mode `mix`, and Bjerrum disabled.
- Allow public `equilibrium(..., options={...})` dicts for legacy PC-SAFT option names. Map safe keys (`max_nfev`, `solver_tol`, `split_tol`, `solver_accept_norm`) and report unsupported legacy keys through diagnostics as `ignored_legacy_options`.
- Add `src/epcsaft/equilibrium_core/legacy_electrolyte_candidate.py` as a diagnostic candidate generator modeled on the old `_solve_two_phase_lle` shape: softmax composition, dependent material-balanced second phase, beta bounds, charge-neutral salt basis, multiple starts, split preference, and anchored/heuristic retry.
- Do not return the legacy candidate as an accepted `EquilibriumResult`. On strict native failure, merge it into `SolutionError.diagnostics` with `legacy_candidate_*` keys and keep `acceptance_gate="predictive_solve_failed"`.

## Test Cases

- Existing Ascani tests still pass and keep the current 1:1 behavior stable.
- New generalized-basis tests cover:
  - 1:1 LiCl or NaCl basis equivalence with current behavior.
  - 2:1 MgCl2 charge-balanced formula reconstruction.
  - 1:2 Na2SO4-style charge-balanced formula reconstruction.
  - mixed Li+/Mg2+/Cl- basis construction without the old 1:1-only rejection.
  - non-neutral feeds and invalid explicit salt labels still raise `InputError`.
- New Hubach tests cover:
  - `2024_Hubach` fixture loading.
  - `2024_Hubach` option normalization matching Lithium's canonical surface.
  - row 0 explicit-seed fixed-species LLE convergence to a distinct split.
  - cold-start failure diagnostics include all seed attempts if no split is found.
  - cold-start strict failure preserves a JSON-safe distinct legacy candidate under `EPCSAFT_RUN_HUBACH_LLE=1`.
  - continuation uses an accepted split as the next point's seed.
- Legacy option tests cover:
  - dict options with old PC-SAFT names map to `EquilibriumOptions`.
  - unsupported old option names are reported in diagnostics instead of silently disappearing.
  - fallback diagnostics distinguish `legacy_candidate_found` from strict equilibrium acceptance.
- Native diagnostics tests cover JSON-safe `seed_attempts` on success and failure.

## Verification Commands

Run focused tests after each implementation layer:

```powershell
uv run python run_pytest.py tests/equilibrium/test_electrolyte_lle.py -q
uv run python run_pytest.py tests/equilibrium/test_hubach_electrolyte_lle.py -q
uv run python run_pytest.py tests/native/test_equilibrium_native_contracts.py -q
```

Final verification:

```powershell
uv run python run_pytest.py tests/equilibrium -q
uv run python run_pytest.py tests/native/test_equilibrium_native_contracts.py -q
uv run python run_pytest.py --confidence -q
```

## Assumptions And Boundaries

- Do not weaken the strict phase-distance acceptance gate. A collapsed phase split remains invalid.
- Divalent-ion support means charge-constrained basis and solver workflow support, not full reactive/speciation equilibrium.
- Hubach convergence proves the general fixed-species electrolyte LLE workflow is more robust; it does not prove exact Table S11 extraction-efficiency reproduction.
- Generated plot/gallery scripts remain out of the default verification path.
