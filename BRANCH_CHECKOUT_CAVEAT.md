# Branch Checkout Caveat for All Codex Agents

Use this in any prompt where the Codex app creates the worktree from `main`.

```text
The Codex app may have created this worktree from `main`. That is okay.

Before Goal Prep writes files or implementation starts, you must switch/create your assigned branch.

Do not write Goal Prep files on `main`.

Run:

git fetch origin --prune

If origin/<assigned-branch> exists, switch to it.
If origin/<assigned-branch> does not exist, create it from origin/main.
If switching or creating the branch fails, stop with BRANCH_BOOTSTRAP_FAILED.

After switching, confirm:

git branch --show-current

The branch must equal the assigned branch.

Only then run Goal Prep.
```
