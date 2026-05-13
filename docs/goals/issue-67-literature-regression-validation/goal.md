# Issue #67 Literature-Backed Regression Validation

Prepare and execute GitHub issue #67 only: add literature-backed regression validation using the repository's existing local assets and publish the focused result through a PR that closes the issue.

## Execution Boundary

This goal starts from the Codex app worktree at `C:\Users\Tanner\.codex\worktrees\fb77\ePC-SAFT`, which the user states was created from an updated local `main`. Before implementation begins, the execution run must checkout or confirm branch `codex/issue-67-literature-regression-validation`.

Do not start implementation until issue #66 and its attached PR are confirmed merged to `origin/main`, local root `main` at `C:\Users\Tanner\Documents\git\ePC-SAFT` is confirmed fast-forwarded to include that merge, and this branch is rebased onto the updated `origin/main` or local `main`. The first execution task (`T001`) must set up or run the watcher/check script for this dependency gate before any issue #67 planning or editing begins.

## Hard Constraints

- Work on GitHub issue #67 only.
- Confirm issue #66 is merged into `main`; if not, stop and report that #67 is blocked.
- Read issue #67 completely, including the Detailed implementation addendum, before planning or editing.
- Do not use PR #56 as a base branch.
- Do not continue or merge PR #56.
- Use PR #56 only as reference if needed.
- Use only existing repo assets: `docs/papers`, `analyses`, `data/reference`, Figiel 2025, and MIAC outputs.
- Do not fetch data from the internet.
- Do not use finite-difference derivatives.
- Do not test vapor Born derivatives.
- Run the validation commands listed in issue #67.
- Open a focused PR linked to issue #67.
- The PR body must include `Closes #67`.
- Stop after this issue.

## Completion Proof

The goal is complete only when the focused issue #67 PR has passed validation and GitHub checks, is squash-merged if permitted, the branch is deleted, root checkout local `main` is fast-forwarded to match `origin/main`, issue #67 is closed, and the final report includes the PR URL, merge commit, validation summary, label status, and final issue state.

If merge or issue closure is not possible, leave the PR ready and report the PR URL plus the exact blocker.

## Required Post-Merge Main Refresh

Before continuing or reporting completion after a successful merge, run from the root checkout context:

```powershell
git -C C:\Users\Tanner\Documents\git\ePC-SAFT pull --ff-only origin main
git -C C:\Users\Tanner\Documents\git\ePC-SAFT rev-parse --short=12 main
git -C C:\Users\Tanner\Documents\git\ePC-SAFT rev-parse --short=12 origin/main
```

The `main` and `origin/main` hashes must match and include the PR merge commit. If this fails, report that future worktrees from main are blocked until local main is updated.

## Starter Command

`/goal Follow docs/goals/issue-67-literature-regression-validation/goal.md.`
