# JetBrains `.idea` Cleanup Tranche for `ePC-SAFT`

## Summary
Implement a deterministic JetBrains metadata cleanup workflow as `scripts/dev/configure_jetbrains_project.py`, not at the repo-root `scripts/` level. The script should normalize all discovered `.iml` files generically, use stable sorted output for relevant root entries, warn about stale module dependencies without failing, and update the plan note plus repo `AGENTS.md` so the documented path matches the implemented one.

## Key Changes
- Add `scripts/dev/configure_jetbrains_project.py` with:
  - `--dry-run` to report pending changes
  - `--apply` to rewrite files in place
  - idempotent behavior and explicit reporting of every proposed/applied change
- Parse `.idea/*.iml` and root-level `*.iml` as XML.
- Never touch `.idea/workspace.xml`.
- Apply cleanup rules across all discovered `.iml` files, not just the main Python module.
- If a module lacks editable content-root structure, create the missing structure even for non-Python module types in this repo, then apply the same cleanup model.
- Enforce a universal root model after cleanup:
  - ensure `src` is present as a non-test source root
  - ensure `tests` is present as a test source root
- Remove transient or generated `sourceFolder` entries under paths such as:
  - `build`, `dist`, `.venv`, `.worktrees`, `_codex`, `.pytest_cache`, `.ruff_cache`, `.mypy_cache`, `docs/_build`, `docs/latex/out`, `results/runs`
- Add `excludeFolder` entries for those generated/transient paths when they exist.
- Generalize redundant nested-exclude cleanup:
  - if a parent path is excluded, remove redundant nested excludes under that parent anywhere in the module, not only under `.worktrees`
- Sort rewritten `sourceFolder` and `excludeFolder` entries deterministically so reruns stay clean and diffs are stable.
- Inspect `.idea/modules.xml` and warn when `.iml` files reference missing module dependencies such as the current `feos` entry, but keep exit code `0` in both dry-run and apply modes.
- Update `docs/plans/PLAN_IDEA_CLEANUP.md` and any affected `AGENTS.md` command text to use the `scripts/dev/` path.

## Current State To Target
- `.idea/ePC-SAFT.iml` currently contains transient `build/uv-cache/.../numpy` source roots that should be removed.
- `src` exists and is already a source root in the main module.
- `tests` exists but is not currently marked as a test source root.
- `.idea/ePC-SAFT.iml` currently contains nested `.worktrees/.../.venv` excludes instead of excluding `.worktrees` broadly.
- `.idea/modules.xml` declares only the local project modules, while `.idea/ePC-SAFT.iml` still references a stale `feos` module dependency.
- `.idea/ePC-SAFT.CMake.iml` currently has minimal XML and should be brought into the same deterministic cleanup model rather than skipped.

## Test Plan
- Pre-check:
  - confirm the transient `build/uv-cache` source roots are present before cleanup
  - confirm `tests` exists and is not yet marked as a test root
  - confirm stale module dependency warnings are expected from the current `feos` entry
- Run:
  - `uv run python scripts/dev/configure_jetbrains_project.py --dry-run`
  - `uv run python scripts/dev/configure_jetbrains_project.py --apply`
  - `uv run python scripts/dev/configure_jetbrains_project.py --dry-run`
- Acceptance criteria:
  - first dry run reports the expected pending changes
  - apply mode reports each actual rewrite
  - second dry run reports no pending changes
  - stale module dependency warnings are visible but do not fail the run
  - `git diff -- .idea scripts/dev/configure_jetbrains_project.py docs/plans/PLAN_IDEA_CLEANUP.md AGENTS.md` shows only intended script/metadata/doc updates
- Optional IDE verification:
  - if IntelliJ is open, call `ide_sync_files` once on the touched `.idea` files after edits
  - use `ide_diagnostics` only if IDE-side metadata problems need confirmation

## Assumptions And Locked Decisions
- The documented path is intentionally changing from `scripts/configure_jetbrains_project.py` to `scripts/dev/configure_jetbrains_project.py`.
- Repo docs update is limited to `PLAN_IDEA_CLEANUP.md` plus affected `AGENTS.md` text; no broader docs sweep is part of this tranche.
- All `.iml` files should be normalized generically.
- Missing editable structure should be created for any discovered module type in this repo.
- `src` and `tests` should be enforced universally after cleanup.
- Redundant nested excludes should be collapsed under any excluded ancestor, not only `.worktrees`.
