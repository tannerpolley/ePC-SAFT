# Branch Bootstrap Instructions for Codex App Worktrees

## Purpose

The Codex app may create a new worktree from `main` even when the intended task branch already exists or should be created.

Therefore every setup prompt and issue-agent prompt must explicitly switch to the assigned branch before Goal Prep writes files or implementation begins.

## Required behavior

At the very beginning of every Codex thread, before Goal Prep creates files, run the branch bootstrap.

For a task branch:

```powershell
git fetch origin --prune
git branch --show-current
```

Then apply this logic:

```text
If current branch is already the assigned branch:
  continue.

If origin/<assigned-branch> exists:
  switch to that branch using the remote branch as source.

If origin/<assigned-branch> does not exist:
  create the assigned branch from origin/main.

If switching/creating the branch fails:
  stop and report BRANCH_BOOTSTRAP_FAILED.
```

## Powershell command template

Replace `<BRANCH>` with the assigned branch.

```powershell
git fetch origin --prune
$current = (git branch --show-current).Trim()
if ($current -ne "<BRANCH>") {
    git ls-remote --exit-code --heads origin <BRANCH>
    if ($LASTEXITCODE -eq 0) {
        git switch <BRANCH> 2>$null
        if ($LASTEXITCODE -ne 0) {
            git switch -c <BRANCH> --track origin/<BRANCH>
        }
    } else {
        git switch -c <BRANCH> origin/main
    }
}
$current = (git branch --show-current).Trim()
if ($current -ne "<BRANCH>") {
    Write-Error "BRANCH_BOOTSTRAP_FAILED expected=<BRANCH> actual=$current"
    exit 1
}
git status --short
```

## Rules

Do not create local `.worktrees/`.

Do not continue on `main` unless the assigned branch is actually `main`.

Do not write Goal Prep files to the wrong branch.

If branch bootstrap fails, stop immediately and report the exact command output.
