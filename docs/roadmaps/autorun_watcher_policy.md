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

## Stop States

Allowed stop states include:

- PREPARED_WAITING
- BLOCKED_DEPENDENCY
- BLOCKED_REBASE_CONFLICT
- BLOCKED_MISSING_INPUT
- WATCHER_TIMEOUT
- IMPLEMENTATION_STARTED
- PR_OPENED
- MERGED

Do not stop at PREPARED_READY only because `auto_start_after_gate` is false. Manual checkpoint mode is retired unless the user explicitly asks for it in a future issue.

## Process Rules

- Never ask repeated confirmation questions.
- Do not create local `.worktrees/`.
- Do not write files on main.
- Use bounded polling only; no infinite watchers.
- If dependencies are satisfied, update from origin/main and start implementation without asking.
- If dependencies are not satisfied, run the bounded watcher and stop with PREPARED_WAITING after timeout.
