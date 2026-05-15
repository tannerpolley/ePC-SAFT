# Intake Notes - Package Clutter Pruning

This note is seeded during board preparation only. The active Phase 0 task must replace placeholders with live command evidence before any cleanup edit.

## Source Authority

- GitHub issue: https://github.com/tannerpolley/ePC-SAFT/issues/120
- Checklist comment: https://github.com/tannerpolley/ePC-SAFT/issues/120#issuecomment-4456462579
- Issue title: Prune generated artifacts, runtime clutter, and staged/optional API noise before vertical implementation work
- Issue state at board creation: open

## Board Prep Snapshot

- Prep branch: `codex/package-cleanup`
- Prep HEAD: `8f0ea57`
- Prep observation: local `main`, `origin/main`, and `origin/HEAD` pointed at `8f0ea57` when the board was created.
- Prep worktree status before board files: clean.

## Phase 0 Evidence To Fill

Record these before product edits:

- origin/main SHA:
- current branch:
- current tracked artifact scan output:
- current import/export surface from `epcsaft.__init__`:
- current list of benchmark modules under `src/epcsaft`:
- current list of top-level staged/optional modules:
- current list of broad exception sites in `src/epcsaft`:
- baseline `uv run python scripts/dev/validate_project.py quick` result:
- baseline `uv run python scripts/dev/validate_project.py docs` result:
- baseline `uv run python run_pytest.py tests/api/package -q` result:
- baseline `uv run python run_pytest.py tests/api/runtime -q` result:
- baseline `uv run python run_pytest.py tests/workflows/repo -q` result:

## Hard Boundaries

- No scientific implementation changes.
- No #114-#119 implementation work.
- No inventory-only closure.
- No staged/debug-route proof counted as production completion.
- No exact banned backend or derivative-status tokens in committed text.
