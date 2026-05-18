# Algorithm Registry Handoff

## Purpose

Create an algorithm traceability registry for `epcsaft` analogous to the current equation registry, focused on package-owned algorithms such as phase equilibrium, stability, chemical/speciation equilibrium, reactive phase equilibrium, and Ceres-backed regression.

This is documentation and traceability work first. Do not change solver behavior, numerical methods, public API semantics, or optimizer coupling while implementing the registry.

## Active Branch And Workspace

Use the `ipopt` branch as the source branch for this work.

Recommended new-thread startup:

```powershell
cd C:\Users\Tanner\Documents\git\ePC-SAFT
git status --short --branch
git log --oneline --decorate -5
git worktree list
```

The manual worktree that was briefly created at:

```text
C:\Users\Tanner\.codex\worktrees\ePC-SAFT-algorithm-registry
```

has been removed. A local branch named `codex/algorithm-registry` may still exist at the same commit as `ipopt` with no unique work. If a new Codex app-created worktree wants that branch name and it collides, verify it has no unique commits before deciding whether to reuse or delete it.

Do not target `main`. If a PR is ever used, target a non-main remote branch such as `origin/ipopt`. For local-only integration, no PR is needed: commit in the isolated worktree branch, then merge locally into `ipopt` or another non-main integration branch.

## Required Context

Read before implementing:

- `AGENTS.md`
- `docs/.codex-journal/user_preferences.md`
- `docs/roadmaps/FULL_ROADMAP.md`
- `C:\Users\Tanner\.codex\PROJECT_ARCHITECTURE.md`
- `C:\Users\Tanner\.codex\instructions\GIT_AND_SANDBOX.md`
- `C:\Users\Tanner\.codex\instructions\COMMAND_CONTEXT_VALIDATION.md`
- `C:\Users\Tanner\.codex\instructions\PROCESS_AND_CLEANUP.md`

Suggested skills:

- `chemical-engineer`
- `jetbrains` for semantic route tracing
- `superpowers:using-git-worktrees` if creating a fresh worktree
- `superpowers:verification-before-completion` before final handoff

## Existing Equation Registry Pattern

The algorithm registry should follow the existing equation-registry shape instead of inventing a separate maintenance style.

Current equation registry anchors:

- Source of truth: `docs/latex/equations.tex`
- Generator: `scripts/docs/sync_equation_registry.py`
- Generated views: `docs/equations.md`, `docs/equations_registry.yaml`
- Native owner comments: `// EqID: ...` under `src/epcsaft/native/**`
- Contract tests: `tests/native/contracts/test_equation_registry.py`

Do not manually edit generated equation files. Do not mix algorithm entries into the equation registry.

## Desired New Artifacts

Add these files unless a better repo-local naming convention emerges during implementation:

- `docs/latex/algorithms.tex`
- `docs/algorithms.md`
- `docs/algorithms_registry.yaml`
- `scripts/docs/sync_algorithm_registry.py`
- `tests/native/contracts/test_algorithm_registry.py`

The generated Markdown should be human-readable. The YAML should be machine-readable and stable enough for future contract tests, capability audits, and agent navigation.

## Registry Source Model

Use `docs/latex/algorithms.tex` as the source of truth. Each algorithm entry should have a stable `AlgID` and metadata comments similar to `EqID` metadata.

Recommended metadata fields:

- `AlgID`
- `Family`
- `Status`
- `Public API`
- `Backend`
- `Dependency`
- `Derivative backend`
- `Solver role`
- `Implementation owner`
- `Validation`
- `Capability key`
- `Description`
- `Change note`

Keep the registry descriptive, not aspirational. If an algorithm is only a diagnostic residual evaluator or orchestration layer, say that. Do not imply production optimizer support unless current code and tests prove it.

## Owner Comment Model

The equation registry currently parses C++ `// EqID:` comments. Algorithm ownership crosses Python, pybind, and native C++, so the algorithm parser should support at least:

```cpp
// AlgID: neutral_tp_flash_ipopt
```

```python
# AlgID: public_equilibrium_problem_dispatch
```

Recommended scan roots:

- `src/epcsaft`
- `src/epcsaft/native`
- `tests`

Do not scan `docs/papers/**` or generated build outputs.

The owner reference should attach to the next non-empty, non-comment code line, mirroring the equation registry behavior.

## Solver Boundary Rule

Ceres and Ipopt must remain separate algorithm families.

- Ceres owns supported native regression least-squares routes.
- Ipopt owns supported native equilibrium/speciation/stability constrained-NLP routes.
- They may both call shared ePC-SAFT state/property/derivative code.
- They must not call each other.
- Do not create registry wording that suggests a nested Ceres-inside-Ipopt or Ipopt-inside-Ceres workflow.

This boundary is intentional.

## Initial Algorithm IDs To Register

Start with a compact set that covers the major package-owned pathways. Keep the first pass broad enough to be useful, but avoid documenting every helper as a top-level algorithm.

### Core Optimizer/Adapter Algorithms

- `native_nlp_problem_contract`
  - Owner: `src/epcsaft/native/equilibrium_nlp/nlp_problem.h`
  - Role: shared native NLP interface for Ipopt routes.

- `ipopt_tnlp_adapter`
  - Owner: `src/epcsaft/native/equilibrium_nlp/ipopt_adapter.cpp`
  - Role: maps `NlpProblem` into `Ipopt::TNLP` and runs `solve_ipopt_nlp`.

- `native_ceres_regression_adapter`
  - Owner: `src/epcsaft/native/epcsaft_regression.cpp`
  - Role: maps native regression residual/Jacobian evaluators into `ceres::CostFunction` and `ceres::Problem`.

### Equilibrium Algorithms

- `neutral_tp_flash_ipopt`
  - Owners: `src/epcsaft/equilibrium.py`, `src/epcsaft/native/equilibrium_nlp/route_builders.cpp`, `src/epcsaft/bindings.cpp`
  - Public API: `mix.equilibrium(kind="tp_flash", ...)`, `mix.solve_equilibrium(TPFlash(...))`

- `neutral_lle_ipopt`
  - Owners: `src/epcsaft/equilibrium.py`, `src/epcsaft/native/equilibrium_nlp/route_builders.cpp`, `src/epcsaft/bindings.cpp`
  - Public API: `mix.equilibrium(kind="lle_flash", ...)`, `mix.solve_equilibrium(LLEProblem(...))`

- `electrolyte_lle_ipopt`
  - Owners: `src/epcsaft/equilibrium.py`, `src/epcsaft/native/equilibrium_nlp/route_builders.cpp`, `src/epcsaft/bindings.cpp`
  - Public API: `mix.equilibrium(kind="electrolyte_lle", ...)`, `mix.solve_equilibrium(ElectrolyteLLEProblem(...))`

- `bubble_dew_ipopt`
  - Owners: `src/epcsaft/equilibrium.py`, `src/epcsaft/native/equilibrium_nlp/bubble_dew_route_builders.cpp`, `src/epcsaft/bindings.cpp`
  - Public API: bubble/dew `kind` routes in `mix.equilibrium(...)`

- `stability_tpd_ipopt`
  - Owners: `src/epcsaft/equilibrium.py`, `src/epcsaft/native/equilibrium_nlp/stability_route_builders.cpp`, `src/epcsaft/bindings.cpp`
  - Public API: `mix.equilibrium(kind="stability", ...)`, `mix.equilibrium(kind="electrolyte_stability", ...)`

### Reactive And Speciation Algorithms

- `ideal_speciation_ipopt`
  - Owners: `src/epcsaft/reactive_speciation.py`, `src/epcsaft/native/equilibrium_nlp/ideal_speciation_problem.cpp`, `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp`, `src/epcsaft/bindings.cpp`
  - Public API: `mix.chemical_equilibrium(...)`, `ReactiveSpeciationProblem`

- `nonideal_speciation_ipopt`
  - Owners: `src/epcsaft/reactive_speciation.py`, `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp`, `src/epcsaft/bindings.cpp`
  - Public API: `mix.chemical_equilibrium(...)`

- `reactive_lle_liquid_root_ipopt`
  - Owners: `src/epcsaft/equilibrium.py`, `src/epcsaft/native/equilibrium/reactive_phase_equilibrium_problem.cpp`, `src/epcsaft/bindings.cpp`
  - Public API: `mix.equilibrium(kind="reactive_lle", ...)`, `ReactivePhaseEquilibriumProblem`

- `reactive_electrolyte_lle_liquid_root_ipopt`
  - Owners: `src/epcsaft/equilibrium.py`, `src/epcsaft/native/equilibrium/reactive_phase_equilibrium_problem.cpp`, `src/epcsaft/bindings.cpp`
  - Public API: `mix.equilibrium(kind="reactive_electrolyte_lle", ...)`

### Regression Algorithms

- `pure_neutral_ceres_regression`
  - Owners: `src/epcsaft/regression.py`, `src/epcsaft/native/epcsaft_regression.cpp`, `src/epcsaft/bindings.cpp`
  - Public API: `fit_pure_parameters(...)`, `fit_pure_neutral(...)`

- `pure_ion_ceres_regression`
  - Owners: `src/epcsaft/regression.py`, `src/epcsaft/native/epcsaft_regression.cpp`, `src/epcsaft/bindings.cpp`
  - Public API: `fit_liquid_electrolyte_parameters(...)` where the supported target set is pure-ion/native Ceres.

- `binary_kij_ceres_regression`
  - Owners: `src/epcsaft/regression.py`, `src/epcsaft/native/epcsaft_regression.cpp`, `src/epcsaft/bindings.cpp`
  - Public API: `fit_binary_parameters(...)` for current native constant-`k_ij` support.
  - Caveat: do not claim native optimizer support for every binary parameter family unless implementation proves it. Existing architecture notes indicate `l_ij` and `k_hb_ij` support must be checked carefully.

- `reactive_electrolyte_batch_residual_context`
  - Owners: `src/epcsaft/reactive_regression.py`, `src/epcsaft/runtime.py`
  - Role: structured residual/context evaluator, not a registered native Ceres production optimizer.

## Generator Requirements

The new `sync_algorithm_registry.py` should:

1. Parse `docs/latex/algorithms.tex`.
2. Parse `# AlgID:` and `// AlgID:` owner comments from selected source roots.
3. Fail on duplicate `AlgID` entries in the LaTeX source.
4. Fail when code owner comments reference unknown `AlgID`s.
5. Attach code references into generated YAML/Markdown.
6. Support `--check`.
7. Support `--strict-traceability`.
8. Keep generated output stable and deterministic.
9. Avoid broad output or full-source dumps in diagnostics.

Reuse structure from `scripts/docs/sync_equation_registry.py` where practical, but do not modify equation-registry semantics.

## Contract Tests

Add `tests/native/contracts/test_algorithm_registry.py` with focused coverage:

- Generated outputs are synced.
- Strict traceability passes.
- Unknown source-code `AlgID` comments fail.
- Documentation-only or planned entries can be exempted only with explicit status.
- Python and C++ owner comments are both parsed.
- Ceres and Ipopt algorithm IDs are present and distinct.
- Generated Markdown names public API entrypoints and backend/dependency roles.

Avoid running broad native builds for this registry-only work unless behavior changes accidentally require it.

## Validation Commands

Minimum validation for a docs/registry-only implementation:

```powershell
uv run python scripts/docs/sync_algorithm_registry.py --check --strict-traceability
uv run python run_pytest.py tests/native/contracts/test_algorithm_registry.py -q
```

If the implementation touches shared docs tooling or existing equation registry code, also run:

```powershell
uv run python scripts/docs/sync_equation_registry.py --check --strict-traceability
uv run python run_pytest.py tests/native/contracts/test_equation_registry.py tests/native/contracts/test_algorithm_registry.py -q
```

Before final completion, run the repo cleanup hook:

```powershell
pwsh.exe -NoProfile -ExecutionPolicy Bypass -File "$env:USERPROFILE\.codex\hooks\codex-cleanup.ps1" -RepoRoot .
```

## Non-Goals

- Do not change Ipopt options, scaling, Hessians, warm starts, diagnostics, or solver math.
- Do not change Ceres optimization behavior.
- Do not add Ceres/Ipopt cross-calls.
- Do not update generated equation docs as part of this work unless a test reveals existing drift and the user approves.
- Do not claim new package capability from the registry alone.
- Do not push or open a PR unless the user explicitly asks.

## Completion Criteria

The work is complete when:

- `docs/latex/algorithms.tex` is the source of truth for algorithm cards.
- `docs/algorithms.md` and `docs/algorithms_registry.yaml` are generated and synced.
- Relevant Python/C++ owners carry `AlgID` comments without changing behavior.
- Strict traceability passes.
- Focused tests pass.
- The final report states exactly which validation ran and whether any native behavior was intentionally untouched.
