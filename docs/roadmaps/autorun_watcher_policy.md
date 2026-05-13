# Bounded Watcher Auto-Run Policy

## Purpose

Roadmap tasks A-M use bounded watcher auto-run by default so prepared agents can wait for dependencies and start implementation when the real gate passes.

## Defaults

```yaml
watcher_mode: bounded
auto_start_after_gate: true
poll_interval_seconds: 120
max_wait_minutes: 480
dependency_gate_path: docs/goals/<slug>/notes/dependency_gate.yaml
watch_script_path: docs/goals/<slug>/notes/watch_dependency.ps1
```

## Gate Source Of Truth

The GitHub Project board is a dashboard only. The dependency gate is:

- dependency issue closed
- dependency PR merged
- origin/main contains the dependency merge
- current branch is the assigned branch
- current branch updates cleanly from origin/main by rebase or fast-forward
- dependency_gate.yaml passes

## Stop States

Allowed stop states include:

- PREPARED_WAITING
- BLOCKED_DEPENDENCY
- BLOCKED_REBASE_CONFLICT
- BLOCKED_MISSING_INPUT
- BLOCKED_CHECKS_FAILED
- BLOCKED_GITHUB_POLICY
- BLOCKED_MERGE_CONFLICT
- BLOCKED_REMOTE_BRANCH_DELETE
- BLOCKED_LOCAL_BRANCH_DELETE
- WATCHER_TIMEOUT
- IMPLEMENTATION_STARTED
- PR_OPENED
- MERGED

Do not stop at PREPARED_READY only because `auto_start_after_gate` is false. Manual checkpoint mode is retired unless the user explicitly asks for it in a future issue.

## PR And Merge Automation

For roadmap tasks F-M, the task prompt pre-authorizes the agent to complete the issue lifecycle without asking for another yes when all gates pass.

Required completion behavior:

- Inspect the corresponding GitHub issue before implementation and record the issue scope in the GoalBuddy notes.
- Open a focused draft PR if one does not already exist for the assigned branch, and include `Closes #<corresponding issue number>` in the PR body.
- Review the PR against `origin/main` before marking ready.
- Confirm changed files match the task scope and do not include finite differences, application-specific public APIs, unrelated files, or silent scope narrowing.
- Mark the PR ready only after task validation, `git diff --check`, and self-review pass.
- Wait for GitHub checks to finish.
- Merge only when checks pass, the branch is current with `origin/main`, the PR is mergeable, and the final GoalBuddy audit says `full_outcome_complete: true`.
- After merge, confirm the corresponding GitHub issue is closed; if still open, close it with a comment naming the merged PR and merge commit.
- Delete both the remote and local task branch. If issue closure or branch cleanup cannot be completed safely, record `BLOCKED_ISSUE_CLOSE_FAILED`, `BLOCKED_REMOTE_BRANCH_DELETE`, or `BLOCKED_LOCAL_BRANCH_DELETE` with the exact next command.

## Process Rules

- Never ask repeated confirmation questions.
- Do not create local `.worktrees/`.
- Do not write files on main.
- Use bounded polling only; no infinite watchers.
- If dependencies are satisfied, update from origin/main and start implementation without asking.
- If dependencies are not satisfied, run the bounded watcher and stop with PREPARED_WAITING after timeout.
