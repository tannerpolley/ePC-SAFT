# T001: Scout Inventory

Task: `T001`
Kind: `scout`
Status: `current`

## Summary

The dependency gate is open in substance for this branch: issue `#95` still lists `Dependencies: none`, `git rev-list --left-right --count origin/main...HEAD` returned `0 0`, and `origin/main` is already an ancestor of `HEAD`. The literal `git rebase origin/main` command was blocked by shell policy, but there is no commit divergence to rebase. The package already has several literature-tied tests and fixtures, but they are fragmented across regression, equilibrium, workflow, and analysis-adjacent paths instead of being presented as one package-owned literature suite.

## Details

- Existing package-owned or package-exercising literature coverage found:
  - `tests/regression/test_literature_pure_parameter_regression.py`
  - `tests/regression/test_literature_binary_kij_regression.py`
  - `tests/regression/test_figiel_2025_born_parameter_parity.py`
  - `tests/regression/test_miac_liquid_electrolyte_regression.py`
  - `tests/equilibrium/test_electrolyte_lle_confidence.py`
  - `tests/equilibrium/test_hubach_electrolyte_lle.py`
- Current benchmark harnesses exist for neutral equilibrium and reactive regression under `src/epcsaft/benchmarks/`, but there is no parallel package-owned literature-suite manifest or runner.
- Scope items already covered at least partially by tests:
  - MEA simple workflow benchmark
  - Figiel 2025 SSM+DS Born benchmark
  - Khudaida salting-out electrolyte LLE benchmark
  - Hubach/Yu lithium-related equilibrium benchmark
- Scope items that still lack a clear package-owned literature-suite entry or are blocked on broader generic routes:
  - MDEA ePC-SAFT benchmark
  - Held 2014 revised ePC-SAFT benchmark
  - non-electrolyte LLE benchmark as an explicit literature anchor
  - Ascani 2022 electrolyte LLE benchmark
  - Ascani 2023 reactive LLE benchmark
- One scope risk is visible immediately: current literature regression coverage still relies on the public helper `fit_mea_co2_h2o_electrolyte`, which is application-specific and conflicts with the issue rule that the package must stay general-purpose.

## Board Receipt Snippet

```yaml
receipt:
  result: done
  note: notes/T001-scout-inventory.md
  summary: "Inventory complete. Existing literature coverage is fragmented; several anchors are already tested, but the repo lacks a package-owned literature-suite manifest/runner and still exposes MEA-specific regression coverage through an application-specific public helper."
```

