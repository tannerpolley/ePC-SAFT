# T058 Figiel 2025 parity handoff

Date: 2026-05-11  
Branch: `codex/issue-53-native-regression-production`  
Goal state: active  
Board truth: `active_task: T058`

## Why this handoff exists

The issue-53 branch already has substantial native Ceres/CppAD work landed, including the recent bubble-derivative tranche. The latest user redirect changes the priority:

1. stop treating bubble as the primary Born+SSM+DS benchmark,
2. use **Figiel 2025 liquid-electrolyte parity** as the main truth surface for `d_born` and `f_solv`,
3. use that parity lane to expose the remaining derivative gap.

This note is the exact restart point for the next thread.

## Current branch state

Last committed checkpoints:

- `604bc8e` `Advance bubble audit`
- `a2b57c2` `Add bubble derivative slice`
- `39f5cdc` `Advance bubble tranche`

Current uncommitted tracked changes:

- `src/epcsaft/benchmarks/native_ceres_thermo_regression.py`
- `src/epcsaft/native/epcsaft_regression.cpp`
- `src/epcsaft/native/regression/thermo_regression.cpp`
- `src/epcsaft/regression.py`
- `src/epcsaft/runtime.py`
- `tests/api/test_runtime.py`
- `tests/native/test_native_ceres_thermodynamic_regression.py`

Current untracked test:

- `tests/regression/test_figiel_2025_born_parameter_parity.py`

Do not assume the old top-level goal text is current. The old sentence about a missing native pressure-composition derivative substrate is stale and false. That substrate is already implemented.

## What is already proven

### 1. Bubble tranche is real

The multi-vapor reactive-electrolyte bubble slice was added and passed focused validation before the Figiel pivot.

Relevant changed files:

- `src/epcsaft/native/regression/thermo_regression.cpp`
- `src/epcsaft/benchmarks/native_ceres_thermo_regression.py`
- `src/epcsaft/runtime.py`
- `tests/native/test_native_ceres_thermodynamic_regression.py`
- `tests/api/test_runtime.py`

Focused validation that passed:

```powershell
uv run python run_pytest.py tests/native/test_native_ceres_thermodynamic_regression.py tests/api/test_runtime.py tests/workflows/test_benchmark_native_regression.py -q
```

Recorded result:

- `56 passed in 64.36s`

That slice should not be discarded, but it is no longer the primary Born+SSM+DS benchmark.

### 2. Figiel 2025 is the correct Born+SSM+DS benchmark surface

Verified from:

- `docs/papers/md/Figiel, Yu, Held - 2025 - Predicting Thermodynamic Properties of Ions in Single Solvents and in Mixe.md`

Key facts:

- SSM+DS is the chosen Born treatment for the paper's main calculations.
- `d_born` is ion-specific.
- `f_solv` is solvent-specific through `f_mix`.
- The direct regression target is **liquid-electrolyte NaBr MIAC behavior** in water, methanol, and ethanol.

### 3. Repo already contains the published Figiel values

Canonical catalog:

- `data/reference/epcsaft_parameters/catalog/2025.csv`

Examples already present:

- `Na+ d_born = 3.445`
- `Cl- d_born = 4.1`
- `Br- d_born = 4.48`
- `I- d_born = 4.985`
- `Li+ d_born = 2.784`
- `K+ d_born = 4.15`
- `H2O f_solv = 1.5`
- `Methanol f_solv = 1.4`
- `Ethanol f_solv = 1.6`

These defaults are also wired in:

- `src/epcsaft/parameters.py`

### 4. Public `fit_pure_ion(...)` already reproduces Figiel `Na+ d_born`

Synthetic MIAC records were generated from the current `2025_Figiel` setup for NaBr in water using:

- `analyses/miac_fits/scripts/validate_miac_fits.py`

Using molalities:

- `0.05, 0.1, 0.2, 0.4, 0.8`

Observed result:

- fitted `d_born ~= 3.445`
- fit success = `True`
- backend = `least_squares_native`
- `jacobian_backend = Backend_unavailable`

Interpretation:

- parity works,
- but the derivative path is still falling back to Backend unavailables.

### 5. Generic liquid-electrolyte native regression can now optimize `f_solv`

Before the patch, native generic regression rejected:

- `f_solv`
- `solvation_factor`

Root cause:

- `src/epcsaft/regression.py` did not include `solvation_factor` in `NATIVE_TARGET_KINDS`
- `src/epcsaft/native/epcsaft_regression.cpp` had no generic target enum/case for `f_solv`

Patched behavior:

- `src/epcsaft/regression.py`
  - added `solvation_factor: 6`
  - alias `f_solv: 6`
  - shifted `k_ij`, `l_ij`, `k_hb_ij` to `7`, `8`, `9`
- `src/epcsaft/native/epcsaft_regression.cpp`
  - added `kGenericTargetFSolv = 6`
  - added the `apply_generic_targets_cpp(...)` case for `f_solv`

After that patch, direct internal liquid-electrolyte regression for water returned:

- `x == [1.5]`
- success = `True`
- backend = `least_squares_native`
- `jacobian_backend = Backend_unavailable`

Interpretation:

- published Figiel `f_solv` parity works,
- but the derivative path is still backend-unavailable-backed.

## New test added

File:

- `tests/regression/test_figiel_2025_born_parameter_parity.py`

Tests inside:

1. `test_fit_pure_ion_recovers_figiel_2025_na_born_radius_from_synthetic_miac`
2. `test_native_generic_liquid_electrolyte_fit_recovers_figiel_2025_water_f_solv`

Purpose:

- freeze the Figiel 2025 parity lane as the liquid-electrolyte Born+SSM+DS benchmark
- expose that parity currently succeeds through a backend-unavailable Jacobian path

## Exact current blocker

The next thread should not chase chemistry first. The immediate problem is just two test failures in the new parity slice.

Command already run:

```powershell
uv run python run_pytest.py tests/regression/test_figiel_2025_born_parameter_parity.py tests/api/test_regression_api.py -q
```

Recorded result:

- `2 failed, 25 passed`

### Failure 1

File:

- `tests/regression/test_figiel_2025_born_parameter_parity.py`

Problem:

- the first test asserts `result.initial_cost > result.cost`
- `FitResult` does **not** expose `initial_cost`

Fix:

- remove that assertion
- keep the meaningful checks:
  - success
  - fitted value close to `3.445`
  - backend/jacobian reporting if desired

### Failure 2

File:

- `tests/api/test_regression_api.py`

Problem:

- expected target kinds are stale after inserting `solvation_factor`
- old expectation:

```python
[6, 7]
```

- new actual values:

```python
[7, 8]
```

Fix:

- update that assertion to `[7, 8]`

## Immediate next commands

1. Fix the two test failures above.

2. Re-run the focused slice:

```powershell
uv run python run_pytest.py tests/regression/test_figiel_2025_born_parameter_parity.py tests/api/test_regression_api.py -q
```

3. Then run the slightly broader confirmation slice:

```powershell
uv run python run_pytest.py tests/regression/test_figiel_2025_born_parameter_parity.py tests/api/test_regression_api.py tests/native/test_native_ceres_thermodynamic_regression.py tests/api/test_runtime.py tests/workflows/test_benchmark_native_regression.py -q
```

4. After that, update the issue-53 notes/board truth to say:

- Figiel parity benchmark is now tracked.
- Native liquid-electrolyte `solvation_factor` target-kind support exists.
- The remaining gap is still real:
  - generic liquid-electrolyte native regression Jacobians still use `Backend_unavailable`
  - the broader CppAD production objective is therefore still incomplete

## What not to do next

- do **not** mark the goal complete
- do **not** reopen the old pressure-composition-substrate complaint
- do **not** treat bubble as the primary Born+SSM+DS benchmark
- do **not** claim CppAD covers the full liquid-electrolyte regression path yet

## Recommended next engineering question

Once the two test failures are fixed and the parity tests are green, the next thread should answer:

> Why does the Figiel liquid-electrolyte native regression path still report `jacobian_backend = Backend_unavailable`, and what exact state/residual call in that path is still missing a production CppAD derivative surface?

That is the real next issue-53 blocker exposed by the Figiel benchmark.

