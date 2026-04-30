# Codex Handoff: ePC-SAFT Equilibrium V3

## Start Here

Open this project folder in the new Codex thread:

```text
C:\Users\Tanner\Documents\git\ePC-SAFT\.worktrees\ePC-SAFT-equilibrium-v3
```

Expected branch:

```powershell
git status --short --branch
# ## codex/equilibrium-v3
```

This is a Git worktree, not a separate clone. The main checkout can stay on `codex/equilibrium-v1-tp-flash` while this worktree stays on `codex/equilibrium-v3`.

## Current Baseline

The worktree was created from the current equilibrium branch after the `.worktrees/` ignore rule was added.

Recent validation already completed in this worktree before it was renamed:

```powershell
uv sync --no-install-project
uv run python scripts/build_epcsaft.py
uv run python run_pytest.py tests/equilibrium -q
uv run python scripts/codex_doctor.py
uv run python scripts/sync_equation_registry.py --check
uv run python run_pytest.py --confidence -q
```

Observed results:

- `tests/equilibrium`: `50 passed`
- `--confidence`: `115 passed`
- equation registry: up to date, all implementation equations had C++ owner comments

## Git/Sandbox Notes

This worktree is already trusted for Git safe-directory checks:

```text
C:/Users/Tanner/Documents/git/ePC-SAFT/.worktrees/ePC-SAFT-equilibrium-v3
```

Plain Git commands should work inside this folder. Under Codex `workspace-write`, Git metadata writes may still need sandbox escalation on the first attempt for commands such as `git add`, `git commit`, `git fetch`, and branch/worktree operations.

Future project-local worktrees should be created from the primary checkout with:

```powershell
.\scripts\create_codex_worktree.ps1 -Name <repo-specific-name> -Branch codex/<branch-name>
```

## Equilibrium V3 Scope

Start as a stabilization/API-contract pass, not a broad electrolyte rewrite.

Current implemented surfaces:

- `src/epcsaft/equilibrium.py`
  - `tp_flash`
  - `lle_flash`
  - `neutral_stability`
- `src/epcsaft/epcsaft.py`
  - `ePCSAFTMixture.equilibrium(kind=...)`
- `tests/equilibrium/*`
  - API validation
  - neutral VLE
  - neutral LLE
  - neutral TPD stability

Useful first tasks:

- Clarify that `EquilibriumPhase.fugacity_coefficient` currently stores natural-log fugacity coefficients, or rename/document/test the behavior.
- Improve diagnostics when `EquilibriumOptions(stability_precheck=False)` is used so skipped stability is not confused with proven physical stability.
- Add tests for `include_phase_diagnostics=True`.
- Strengthen `to_dict()` schema and JSON-compatibility tests.

## Avoid

- Do not touch `docs/plots/**` or gallery generation unless explicitly asked.
- Do not run clean native build repair actions while another thread may be importing `epcsaft._core`.
- Do not move/remove worktrees with active submodules unless coordinated.

## Suggested New-Thread First Commands

```powershell
git status --short --branch
uv run python scripts/codex_doctor.py
uv run python run_pytest.py tests/equilibrium -q
```

Before claiming a handoff-ready state:

```powershell
uv run python scripts/sync_equation_registry.py --check
uv run python run_pytest.py --confidence -q
```
