# Placeholder Branch Setup Commands

Run only after the setup PR containing roadmap docs/prompts has been merged into `origin/main`.

```powershell
git fetch origin --prune
git checkout main
git pull --ff-only origin main
```

```powershell
git branch codex/backend-coverage-hard-gate origin/main
git push -u origin codex/backend-coverage-hard-gate
```

```powershell
git branch codex/cppad-explicit-parameter-derivatives origin/main
git push -u origin codex/cppad-explicit-parameter-derivatives
```

```powershell
git branch codex/generic-implicit-sensitivity-framework origin/main
git push -u origin codex/generic-implicit-sensitivity-framework
```

```powershell
git branch codex/reaction-constant-conventions origin/main
git push -u origin codex/reaction-constant-conventions
```

```powershell
git branch codex/generic-target-row-schema origin/main
git push -u origin codex/generic-target-row-schema
```

```powershell
git branch codex/generic-activity-speciation origin/main
git push -u origin codex/generic-activity-speciation
```

```powershell
git branch codex/generic-vle-fugacity-equilibrium origin/main
git push -u origin codex/generic-vle-fugacity-equilibrium
```

```powershell
git branch codex/generic-non-electrolyte-lle origin/main
git push -u origin codex/generic-non-electrolyte-lle
```

```powershell
git branch codex/generic-electrolyte-lle origin/main
git push -u origin codex/generic-electrolyte-lle
```

```powershell
git branch codex/generic-reactive-lle origin/main
git push -u origin codex/generic-reactive-lle
```

```powershell
git branch codex/generic-regression-backend origin/main
git push -u origin codex/generic-regression-backend
```

```powershell
git branch codex/literature-benchmark-suite origin/main
git push -u origin codex/literature-benchmark-suite
```

```powershell
git branch codex/downstream-integration-smokes origin/main
git push -u origin codex/downstream-integration-smokes
```

## Important

Branches created before the setup PR merges will not contain the prompt registry. Do not create implementation branches from stale main.

Every agent must rebase before implementation:

```powershell
git fetch origin --prune
git rebase origin/main
```

If rebase conflicts occur, stop and report.


## Codex worktree caveat

Even if these remote branches exist, the Codex app may create the worktree from `main`.

That is acceptable only if the agent then runs the branch bootstrap in:

```text
docs/roadmaps/branch_bootstrap_instructions.md
```

before Goal Prep writes files.

Agents must not write Goal Prep files or implementation edits while still on `main`.
