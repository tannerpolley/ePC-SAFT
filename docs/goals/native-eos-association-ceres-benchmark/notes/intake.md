# Intake Notes - Native EOS, Association Sensitivities, And Ceres Benchmark

This note is the Stage 0 gate required by GitHub issue #114 before source edits. It exists now as a board setup artifact, but Stage 0 is not complete until every required evidence field below is populated by the active GoalBuddy task.

## Source Authority

- GitHub issue: https://github.com/tannerpolley/ePC-SAFT/issues/114
- Issue title: Complete native EOS/AD, association implicit sensitivities, and associating-binary Ceres regression
- Issue state at board creation: open
- Issue updated at board creation: 2026-05-15T04:15:06Z

## Board Prep Snapshot

- Prep branch: `issue#114`
- Prep HEAD: `255ca07fbce7a2842aaea353fd626ab297af6067`
- origin/main SHA at board creation: `255ca07fbce7a2842aaea353fd626ab297af6067`
- Local main SHA at board creation: `255ca07fbce7a2842aaea353fd626ab297af6067`
- Worktree status before board files: clean
- JetBrains MCP: initially unavailable, then reachable after the IDE was opened. Board-prep `ide_index_status` reported `isDumbMode=false` and `isIndexing=false`.
- JetBrains MCP board-prep file-index check:
  - Found `src/epcsaft/native/epcsaft_ares.cpp`.
  - Did not find `src/epcsaft/native/epcsaft_properties.cpp`.
  - Found `src/epcsaft/native/epcsaft_regression.cpp`.
  - Found `src/epcsaft/regression.py`.
  - Found `docs/pages/parameter_regression.rst`.
  - For `implicit_sensitivity`, found Python/API/test paths but no native file path in the first IDE-index result page.

## Stage 0 Evidence To Fill Before Source Edits

### Current branch

- Required value: not yet filled by Stage 0 executor.
- Board-prep observation: `issue#114`

### origin/main SHA

- Required value: not yet confirmed by Stage 0 executor after fresh fetch.
- Board-prep observation: `255ca07fbce7a2842aaea353fd626ab297af6067`

### C++ EOS contribution files inspected

- Required value: not yet filled.
- Minimum expectation: list the concrete native contribution and EOS/property files inspected, with short notes on scalar-templated and CppAD readiness for issue #114.
- Board-prep IDE note: `epcsaft_ares.cpp` exists; `epcsaft_properties.cpp` was not found by `ide_find_file` and must be reconciled during Stage 0.

### Association solver files inspected

- Required value: not yet filled.
- Minimum expectation: list the concrete association site-fraction solver files inspected, the solved-state variables found, and where implicit sensitivity must enter.
- Board-prep IDE note: native implicit-sensitivity paths were not confirmed by `ide_find_file`; Stage 0 must use IDE semantic lookup or `rg` to locate the actual native owner.

### Native regression files inspected

- Required value: not yet filled.
- Minimum expectation: list the concrete native regression and Python API files inspected, including ownership of residual packing, optimizer loop, result fields, and Python/native boundary.
- Board-prep IDE note: `src/epcsaft/native/epcsaft_regression.cpp`, `src/epcsaft/regression.py`, and regression tests were visible in the IDE file index.

### Current failing or skipped associating-binary regression path

- Required value: not yet filled.
- Minimum expectation: identify the exact current test, API path, skipped case, missing path, or failing reproducer that blocks associating-binary `k_ij` regression.

### Chosen benchmark fixture

- Required value: not yet chosen.
- Preferred fixture from issue #114: MEA + water binary VLE or excess-property data.
- Acceptable substitute from issue #114 if MEA-water data are not normalized: water + associating alcohol binary VLE or LLE data.
- Minimum expectation: record the repo-contained data source, target family, expected fit target, and why the fixture proves association sensitivity participation.

## Commands Stage 0 Should Run Or Justify

```powershell
git fetch origin --prune
git status --short
git rev-parse --abbrev-ref HEAD
git rev-parse origin/main
rg --files src/epcsaft/native
rg -n "CppAD|implicit_sensitivity|association|k_ij|k_hb_ij|l_ij|Ceres" src tests docs
uv run python scripts/dev/doctor.py
```

Stage 0 must add any narrower commands it actually used, including focused tests or build checks that reveal the current failing or skipped associating-binary regression path.

## Hard Boundaries

- Source edits are blocked until Stage 0 fills every required issue #114 intake field above.
- Main thread owns clean or repair builds for `_core`.
- Do not broaden package capabilities before executable proof exists.
- Do not add downstream-application public APIs.
- If files outside the issue-owned paths are needed, Judge must explicitly approve the scope expansion before edits.
