# PLAN — JetBrains `.idea` Metadata Cleanup For `ePC-SAFT`

## Objective

Add a deterministic repo-owned cleanup workflow for JetBrains project metadata so agents do not manually edit `.idea` XML and IntelliJ indexing stays focused on real source roots instead of transient build or cache paths.

## Scope For The Follow-On Agent

Implement a repo-owned script at `scripts/dev/configure_jetbrains_project.py` and keep this work separate from the user-level Codex migration tranche.

The script should:

- run from repo root;
- support `--dry-run` and `--apply`;
- parse `.idea/*.iml` and root-level `*.iml` as XML;
- avoid touching `.idea/workspace.xml`;
- be idempotent and print every proposed or applied change;
- preserve unrelated XML as much as practical.

## Required Behavior

Ensure existing `src` is a source root and existing `tests` is a test source root.

Remove bad source-root entries under transient or generated paths such as:

- `build`
- `dist`
- `.venv`
- `.worktrees`
- `_codex`
- `.pytest_cache`
- `.ruff_cache`
- `.mypy_cache`
- `docs/_build`
- `docs/latex/out`
- `results/runs`

Add exclude-folder entries for those same generated paths when they exist.

If `.worktrees` is excluded, remove redundant nested exclusions beneath it.

Warn about stale module dependencies that no longer exist in `.idea/modules.xml`, but do not auto-remove them.

## Validation

Run:

```powershell
uv run python scripts/dev/configure_jetbrains_project.py --dry-run
uv run python scripts/dev/configure_jetbrains_project.py --apply
uv run python scripts/dev/configure_jetbrains_project.py --dry-run
```

The second dry run must report no pending changes.

Then inspect:

```powershell
git diff -- .idea scripts/dev/configure_jetbrains_project.py
```

## Repo Policy Updates To Include In That Tranche

- update `AGENTS.md` only if the script path or usage text changes materially;
- add any repo-local Codex rule or approval policy only if it is needed for the new script path;
- do not expand the repo into broad JetBrains automation beyond metadata cleanup.
