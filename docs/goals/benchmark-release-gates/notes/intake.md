# Issue 119 Intake

Source authority: GitHub issue #119, "Convert benchmark inventory and downstream smokes into executable release gates"

Goal root: `docs/goals/benchmark-release-gates`

## Phase 0 Baseline

### Branch and issue state

- `verified`: GitHub issue `#119` is open and had no comments when this intake was prepared.
- `verified`: current branch is `issue#119`.
- `verified`: current branch HEAD is `869e3354`.
- `verified`: `git status --short` shows only the new untracked `docs/goals/benchmark-release-gates/` goal root created for this board.
- `verified`: JetBrains MCP was probed once for this worktree and was unavailable because IntelliJ has the main checkout open, not `C:\Users\Tanner\.codex\worktrees\e3ab\ePC-SAFT`.

### Current issue-owned benchmark registry and command surfaces

- `verified`: `src/epcsaft/benchmarks/` currently contains:
  - `src/epcsaft/benchmarks/__init__.py`
  - `src/epcsaft/benchmarks/literature.py`
  - `src/epcsaft/benchmarks/neutral_equilibrium.py`
  - `src/epcsaft/benchmarks/reactive_regression.py`
- `verified`: `src/epcsaft/benchmarks/literature.py` is a compatibility shim that re-exports `scripts.benchmarks.helpers.literature`.
- `verified`: `scripts/benchmarks/` currently contains:
  - `scripts/benchmarks/benchmark_literature_suite.py`
  - `scripts/benchmarks/benchmark_neutral_equilibrium.py`
  - `scripts/benchmarks/benchmark_reactive_regression.py`
  - `scripts/benchmarks/profile_regression_runtime.py`
  - helper modules under `scripts/benchmarks/helpers/`
- `verified`: `scripts/validation/` currently contains:
  - `scripts/validation/validate_electrolyte_lle_confidence.py`
  - `scripts/validation/validate_hydrocarbon_regression.py`
  - `scripts/validation/equilibrium_core/confidence.py`
  - `scripts/validation/equilibrium_core/thermo_diagnostics.py`
- `verified`: `tests/workflows/benchmarks/` currently contains:
  - `tests/workflows/benchmarks/test_benchmark_literature_suite.py`
  - `tests/workflows/benchmarks/test_benchmark_neutral_equilibrium.py`
  - `tests/workflows/benchmarks/test_benchmark_reactive_regression.py`
- `verified`: `tests/regression/literature/` currently contains focused literature tests for ethanol/water binary VLE, Figiel/Held electrolyte checks, literature binary `k_ij`, literature pure-parameter regression, the MEA-CO2-H2O pure benchmark fixture, and Baygi Table 2 pure MEA association.
- `verified`: `tests/fixtures/literature/` currently contains only:
  - `tests/fixtures/literature/binary_kij/ethanol_water_jced2021_100kpa.json`
  - `tests/fixtures/literature/figiel_2025/miac_liquid_electrolyte.json`
  - `tests/fixtures/literature/pure_neutral/mea_co2_h2o_benchmark.json`

### Current benchmark schema reality

- `verified`: `scripts/validation/equilibrium_core/confidence.py` already defines a `BenchmarkCase` dataclass, but that schema is specific to the Khudaida confidence workflow rather than the issue `#119` literature-release gate.
- `verified`: `scripts/benchmarks/helpers/literature.py` currently uses a `LiteratureBenchmarkEntry` dataclass with fields:
  - `case`
  - `title`
  - `classification`
  - `coverage_kind`
  - `package_surface`
  - `validation_paths`
  - `notes`
- `verified`: the current literature helper hard-codes `payload["issue"] == 95` and models the suite as an inventory/classification table, not as the executable `BenchmarkCase(id, source, model_setup, input_records, expected, tolerances, command)` contract required by issue `#119`.

### Live baseline commands

- `verified`: `uv run python scripts/benchmarks/benchmark_literature_suite.py --case figiel_2025_ssm_ds_born`
  - current output is a one-row inventory/classification table:
    - case `figiel_2025_ssm_ds_born`
    - classification `already_supported_with_tests`
    - coverage kind `smoke_regression`
    - validation path count `1`
- `verified`: `uv run python run_pytest.py tests/workflows/benchmarks -q`
  - current result: `14 passed in 68.50s`
- `verified`: `uv run python run_pytest.py tests/api/package -q`
  - current result: `6 passed in 1.02s`

### Current numeric tolerance and pass/fail surfaces

- `verified`: `tests/regression/literature/test_ethanol_water_binary_vle_regression.py` asserts:
  - `result.objective_final < result.objective_initial`
  - fitted `k_ij` moves from the initial value
  - `binary_vle_fugacity_balance < 0.04`
  - fitted `k_ij` remains within `0.01` of the paper reference
- `verified`: `tests/regression/literature/test_figiel_held_electrolyte_benchmarks.py` asserts:
  - aqueous NaCl osmotic coefficient agreement within absolute tolerance `0.03`
  - aqueous and mixed-solvent MIAC agreement within absolute tolerance `0.08`
  - `result.objective_final <= result.objective_initial` for the Figiel regression probe
  - nonzero parameter movement for `d_born` and `f_solv`
- `verified`: `tests/regression/literature/test_literature_pure_parameter_regression.py` asserts:
  - local fixture provenance is preserved
  - backend and Jacobian route stay native/non-forbidden
  - fit targets stay on the intended pure-parameter families
- `verified`: `scripts/validation/equilibrium_core/thermo_diagnostics.py` and `scripts/validation/equilibrium_core/confidence.py` already contain explicit tolerance and threshold logic for Khudaida-style electrolyte LLE validation, but that logic is not yet routed through `tests/workflows/benchmarks/` or the issue `#119` literature suite.
- `inference`: the repo has multiple real tolerance-bearing checks already, but they are fragmented across regression tests, workflow tests, and validation helpers instead of being assembled into one release-gate command surface.

## Benchmark Family Coverage Matrix

### 1. Gross/Sadowski pure PC-SAFT nonassociating parameters

- `verified`: current anchor exists in `tests/regression/core/test_hydrocarbon.py`.
- `verified`: `scripts/validation/validate_hydrocarbon_regression.py` prints a Gross/Sadowski Table 2 target summary.
- `inference`: this family already has executable coverage, but it is not yet represented in the issue `#119` literature suite contract.

### 2. Gross/Sadowski associating systems

- `verified`: current executable associating-binary regression anchor exists in `tests/regression/literature/test_ethanol_water_binary_vle_regression.py`.
- `verified`: repo-local issue `#114` closeout recorded that the selected associating-binary benchmark is the repo-contained ethanol/water 100 kPa VLE fixture because no normalized MEA/water binary VLE fixture was found.
- `inference`: the issue `#119` registry should treat ethanol/water as the current executable associating-system surrogate unless a normalized repo-contained Gross/Sadowski-style associating benchmark is added.

### 3. Baygi MEA association and MEA-water binary baseline

- `verified`: pure MEA association coverage exists in `tests/regression/literature/test_mea_table2_associating_pure.py`.
- `verified`: `tests/fixtures/literature/pure_neutral/mea_co2_h2o_benchmark.json` points at the Baygi 2015 local asset and feeds the MEA/CO2/H2O pure benchmark path.
- `verified`: no normalized repo-contained MEA/water binary VLE fixture was found in the current checkout; the prior issue `#114` report records that same conclusion.
- `inference`: the Baygi pure-association piece is executable now, but the MEA/water binary baseline remains a fixture gap rather than a proven executable issue `#119` gate.

### 4. Cameretti/Held aqueous electrolyte density and MIAC

- `verified`: executable checks exist in `tests/regression/literature/test_figiel_held_electrolyte_benchmarks.py`.
- `verified`: current reference files include:
  - `data/reference/osmotic/water/NaCl.csv`
  - `data/reference/MIAC/water/water-NaCl.csv`
- `inference`: this family already has a real numeric benchmark path and should be promoted into the issue `#119` registry instead of being left behind the issue `#95` inventory layer.

### 5. Held alcohol/salt mixed-solvent density, osmotic coefficient, and MIAC

- `verified`: executable mixed-solvent Held coverage exists in `tests/regression/literature/test_figiel_held_electrolyte_benchmarks.py`.
- `verified`: current reference file includes `data/reference/MIAC/water-methanol/water-methanol-NaCl.csv`.
- `verified`: the current test checks MIAC numerically and constrains osmotic coefficient to a finite physical band, but does not yet report the full issue `#119` benchmark fields in a release report.

### 6. Bulow/Ascani concentration-dependent dielectric and Born behavior

- `verified`: dataset and paper-parameter extraction surfaces exist for `2020_Bulow`, `2021_Bulow`, and `2025_Figiel` in `scripts/data/extract_paper_parameter_csvs.py`.
- `verified`: figure-owned analysis assets exist under `analyses/paper_validation/native/2025_figiel/`.
- `inference`: current issue-owned surfaces expose data extraction and figure validation anchors, but no dedicated issue `#119` executable benchmark command or registry contract was found for the Bulow/Ascani dielectric/Born family.
- `blocked current gap`: executable issue `#119` benchmark wiring is missing in the current package release-gate suite even though parameter/data assets exist.

### 7. Figiel 2025 modified Born / SSM / DS

- `verified`: executable checks exist in:
  - `tests/regression/literature/test_figiel_2025_born_parameter_parity.py`
  - `tests/regression/literature/test_figiel_held_electrolyte_benchmarks.py`
  - `tests/regression/electrolyte/test_miac_liquid_electrolyte_regression.py`
  - `tests/regression/electrolyte/test_miac_liquid_electrolyte_parity.py`
- `verified`: current fixture file is `tests/fixtures/literature/figiel_2025/miac_liquid_electrolyte.json`.
- `verified`: the current literature-suite inventory marks this family as `already_supported_with_tests`, but only as a classification row rather than a full issue `#119` executable contract.

### 8. Ascani 2022 distributed-ion electrolyte LLE

- `verified`: roadmap and issue surfaces point to this family, and issue `#116` is the production distributed-ion LLE implementation issue.
- `verified`: current repo data includes `data/reference/multiphase/ascani_case2_model_comparison.md`.
- `verified`: no executable issue-owned benchmark test or fixture was found under `tests/workflows/benchmarks/`, `tests/regression/literature/`, or `tests/fixtures/literature/` for an Ascani 2022 distributed-ion LLE release gate.
- `inference`: production capability may exist after closed issue `#116`, but issue `#119` still lacks the fixture-to-command-to-report wiring in its owned validation surfaces.

### 9. Ascani 2023 reactive phase equilibrium

- `verified`: roadmap and issue surfaces point to this family, and issue `#117` is the production coupled reactive phase-equilibrium implementation issue.
- `verified`: the current issue `#95` literature inventory still classifies `ascani_2023_reactive_lle` as `blocker_requires_followup`.
- `verified`: no executable issue-owned benchmark test or fixture was found under `tests/workflows/benchmarks/`, `tests/regression/literature/`, or `tests/fixtures/literature/` for an Ascani 2023 reactive phase-equilibrium release gate.
- `inference`: issue `#119` needs new benchmark wiring and likely a repo-contained fixture/report surface even if core package capability exists after closed issue `#117`.

### 10. Khudaida 2026 salting-out LLE

- `verified`: executable Khudaida validation logic exists in:
  - `scripts/validation/equilibrium_core/confidence.py`
  - `scripts/validation/equilibrium_core/thermo_diagnostics.py`
  - `tests/workflows/validation/equilibrium_core/test_electrolyte_lle_confidence.py`
  - `tests/workflows/validation/equilibrium_core/test_electrolyte_thermo_diagnostics.py`
- `verified`: the current issue `#95` literature inventory marks `khudaida_2026_salting_out_lle` as already supported, but the path is still outside the issue `#119` `tests/workflows/benchmarks/` surface.
- `inference`: the main gap is release-gate consolidation, not obvious missing core package behavior.

### 11. Rezaee lithium extraction thermodynamic model inputs

- `verified`: no current issue-owned package benchmark fixture or executable literature test was found for Rezaee in this checkout.
- `verified`: the strongest current Rezaee command surfaces are downstream in `C:\Users\Tanner\Documents\git\Lithium_Extraction`.
- `inference`: this family is a real downstream integration requirement for issue `#119`, not a package-local literature benchmark already represented in the current issue-owned surfaces.

### 12. MEA true-species pressure/speciation workflow fixture

- `verified`: no current issue-owned package fixture was found under `tests/fixtures/literature/` for a true-species MEA pressure/speciation workflow release gate.
- `verified`: issue `#115` is the upstream production issue for native coupled activity speciation and reactive VLE pressure benchmarks, and it is currently closed.
- `inference`: the package capability owner is upstream, but issue `#119` still lacks either a package-local executable workflow fixture or a wired downstream integration command/report proving the capability through a real repo.

## Package-install and downstream-proof reality

- `verified`: `tests/api/package/test_downstream_integration_smokes.py` is package-local and synthetic.
- `verified`: that file constructs `ReactiveSpeciationProblem`, `ElectrolyteLLEProblem`, and `ReactiveElectrolyteBubbleProblem` objects, but it checks synthetic payloads returned by `_generic_result_payload(...)` and `_generic_chemical_payload(...)` rather than running real downstream repositories.
- `verified`: the test is useful as a generic contract guard, but it does not satisfy issue `#119` downstream proof by itself.
- `verified`: `docs/pages/downstream_local_installs.rst` exists and documents editable/path installs plus local build-directory/Ceres reuse behavior.
- `verified`: `docs/roadmaps/release_benchmark_report.md` does not exist.
- `verified`: `docs/roadmaps/downstream_integration_report.md` does not exist.

## Candidate downstream workflow commands

### MEA-Thermodynamics

- `verified` fast install/integration check:
  - `uv run python scripts/check_epcsaft_integration.py --mode dev`
- `verified` heavier repo workflow candidate:
  - `uv run python analyses/epcsaft_ionic_regression/scripts/generate_data.py`
- `inference`: the integration script is the fastest lane to prove local install resolution and generic public API presence; the ionic-regression analysis script is the better candidate for the later “real downstream workflow” proof.

### Lithium_Extraction

- `verified` fast install/integration check:
  - `uv run python scripts/check_epcsaft_integration.py --mode dev`
- `verified` heavier repo workflow candidates:
  - `uv run python analyses/rezaee_2026_pcsaft_epcsaft/scripts/rezaee_reactive_equilibrium_replay.py`
  - `uv run python analyses/rezaee_2026_pcsaft_epcsaft/scripts/rezaee_tds_li_oa_calibrated_surrogate.py`
- `inference`: `rezaee_reactive_equilibrium_replay.py` is the cleaner candidate when issue `#119` needs a thermodynamics-centered downstream proof rather than a broader surrogate/costing artifact.

### MEA-Absorption-Column

- `verified` fast install/integration check:
  - `uv run python scripts/check_epcsaft_integration.py --mode dev`
- `verified` heavier repo workflow candidate:
  - `.\.venv\Scripts\python.exe analyses\nccc_validation\scripts\run_epcsaft_electrolyte_config_matrix.py`
- `inference`: the config-matrix workflow is the stronger “real downstream command” because `analysis.yaml` marks it as ePC-SAFT-dependent and manuscript-facing.

## Current blockers and proceed decision

- `verified`: the main issue `#119` gap is not simply “no benchmarks exist.” The repo already has benchmark fragments, regression tests, workflow tests, validation helpers, and downstream contract checks.
- `verified`: the release-gate gap is that these surfaces are still split across:
  - issue `#95` inventory/classification tables
  - standalone regression tests
  - workflow tests outside `tests/workflows/benchmarks/`
  - package-local synthetic downstream smokes
  - downstream repos that are not yet wired into package-side reports
- `verified`: no behavior-changing edits to the forbidden core implementation paths are needed to begin Phase 1 registry work.
- `inference`: Phase 1 can proceed inside the issue-owned validation surfaces by replacing the issue `#95` inventory-only literature suite with an honest issue `#119` benchmark registry contract that:
  - promotes the already-executable families into explicit case records
  - records blocked families explicitly instead of classifying them as “supported”
  - leaves missing package behavior routed to the owning upstream issue when that case arises

## Scope Guard

Do not use issue `#119` to change:

- `src/epcsaft/native/**`
- `src/epcsaft/equilibrium.py`
- `src/epcsaft/regression.py`
- `src/epcsaft/reactive_speciation.py`

except for import or path fixes that do not change package behavior.

## Required Validation Later

```powershell
uv run python scripts/benchmarks/benchmark_literature_suite.py
uv run python run_pytest.py tests/workflows/benchmarks -q
uv run python run_pytest.py tests/regression/literature -q
uv run python run_pytest.py tests/api/package -q
uv run python scripts/dev/validate_project.py quick
uv run python scripts/dev/validate_project.py docs
git diff --check
```
