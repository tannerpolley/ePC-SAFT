# Graphify-Guided Codebase Improvement Handoff

## Purpose

Use the local Graphify knowledge graph to guide a cleanup-oriented architecture pass in a separate worktree. The goal is to improve navigation, reduce legacy/shared-helper confusion, tighten docs and tests around the current project structure, and avoid changing scientific behavior unless an existing contract requires it.

## Current State

- This repository is a package repo: package code lives in `src/epcsaft`, package/API/native/workflow tests live in `tests`, reusable checkout data lives in `data/reference`, and scientific workflows live in `analyses`.
- The local Graphify graph was built from commit `378edd7f` and reports 1965 nodes, 4809 edges, and 106 communities.
- Local Graphify artifacts live under `graphify-out/`, including `GRAPH_REPORT.md`, `graph.html`, `GRAPH_TREE.html`, and `graph.json`.
- Treat `graphify-out/` as local generated architecture context. Do not commit it unless the user explicitly asks for versioned graph artifacts.

## Graphify Findings To Preserve

- The public API/native bridge is central. Changes touching `ePCSAFTMixture`, `ePCSAFTState`, density closure, pybind bindings, or native C++ need focused API and native contract tests.
- Analysis folders are still active users of shared helpers, not purely archival folders. Cleanup must preserve analysis imports, output paths, and self-contained analysis workflows.
- `scripts/plot_outputs.py` still acts as a cross-analysis bridge. Consider a neutral helper name such as `scripts/analysis_outputs.py`, but keep compatibility or update all call sites in one pass.
- Parameter and reference-data loading is more central than it looks. Keep `data/reference` behavior documented and covered by tests before moving or renaming anything around datasets.
- Generic graph god nodes such as `main()` and `ValueError()` are noisy. Use Graphify for navigation and relationship discovery, then verify with `rg`, file reads, and focused tests.

## Recommended Implementation Order

1. Confirm the main checkout is clean except expected local/generated artifacts, then create an isolated worktree branch such as `codex/graphify-codebase-cleanup`.
2. Read `graphify-out/GRAPH_REPORT.md` and run two or three targeted Graphify queries before editing, for example:
   - `uv run graphify query "Which helpers connect analyses to package outputs?" --graph graphify-out/graph.json --budget 2500`
   - `uv run graphify query "Which modules connect public Python APIs to native C++ density and equilibrium code?" --graph graphify-out/graph.json --budget 2500`
   - `uv run graphify query "Where do docs still reference legacy analysis output paths?" --graph graphify-out/graph.json --budget 2500`
3. Audit stale docs and scripts references with:
   - `rg "scripts/paper_validation|scripts/fits|docs/plots|plot gallery|plot_manifest|gallery server|/out/|\\\\out\\\\"`
   - `rg "plot_outputs|paper_validation_output_path|analysis_final_path"`
4. Make the smallest docs or architecture cleanup first. Prefer clarifying current structure before renaming helpers.
5. If replacing or wrapping `scripts/plot_outputs.py`, add `scripts/analysis_outputs.py` and either preserve a compatibility import layer or update every call site and test in the same branch.
6. Add or update focused tests only where imports, public helper names, command paths, or output-path behavior changes.
7. Run `graphify update .` after code changes so the local graph remains useful for the next pass.
8. Commit tracked source/docs/test changes. Leave `graphify-out/` untracked unless the user explicitly asks otherwise.

## Non-Goals

- Do not redesign the public thermodynamic API in this cleanup pass.
- Do not move native solver logic between Python and C++ unless a failing issue or test specifically requires it.
- Do not reintroduce plot gallery, plot server, plot index, plot manifest, or long default plot-generation workflows.
- Do not add broad package-owned batch/cache APIs such as `batch=True`, `FugacityCache`, or many-point fugacity helpers.
- Do not treat Graphify inferred edges as proof. They are navigation hints that need source verification.

## Acceptance Criteria

- `git status --short` has no accidental generated noise staged or committed.
- Docs describe the current package/analysis/data architecture without stale `out`, gallery, server, or manifest language.
- Any helper rename has a compatibility path or all call sites are updated.
- Focused tests pass for any changed behavior.
- `uv run python scripts/validate_project.py quick` passes before handoff.
- If Sphinx docs changed, `uv run python scripts/validate_project.py docs` passes.
- If native C++ or pybind files changed, `uv run python scripts/build_epcsaft.py` passes before tests.

## Suggested Validation Ladder

- For helper/docs cleanup:
  - `uv run python run_pytest.py tests/workflows -q`
- For package/API touchpoints:
  - `uv run python run_pytest.py tests/api/test_runtime.py tests/native/test_runtime_contracts.py -q`
- For final repo health:
  - `uv run python scripts/doctor.py`
  - `uv run python scripts/validate_project.py quick`
  - `uv run python scripts/validate_project.py docs` if docs changed
  - `uv run ruff check .`
  - `uv run black --check .`

## Handoff Notes For The Next Agent

- Work in a separate worktree. The main checkout should remain the coordination and merge surface.
- Keep the first PR small if possible: docs and navigation cleanup before behavior changes.
- Coordinate native rebuilds through the main thread if any C++ or pybind files are edited.
- Do not commit ignored local-only files such as `AGENTS.md`, `.codex/**`, or generated Graphify outputs.
- If stale references are intentionally retained for compatibility, document why and add a follow-up note rather than leaving them ambiguous.
