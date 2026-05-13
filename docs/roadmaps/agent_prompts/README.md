# Agent Prompt Registry

Use `index.yaml` to find the prompt for a task key A-M.

Each prompt follows this pattern:

```text
Goal Prep first.
No repeated confirmation questions.
Create GoalBuddy files.
Write `docs/goals/<slug>/notes/dependency_gate.yaml`.
Use bounded watcher auto-run by default. Start implementation immediately when the dependency gate passes. If the gate does not pass, create the per-goal watcher, poll every 120 seconds, and stop after 480 minutes with `PREPARED_WAITING`.

Do not stop at `PREPARED_READY` merely because `auto_start_after_gate` is false. Manual checkpoint mode is retired unless the user explicitly asks for it in a future issue.
```

Tasks F-M also include the full completion lifecycle: open a focused draft PR, self-review against `origin/main`, mark ready after validation, wait for checks, merge without another yes only when checks pass and the final GoalBuddy audit says `full_outcome_complete: true`, then delete both the remote and local task branches. If GitHub, mergeability, or branch cleanup blocks that lifecycle, record the explicit `BLOCKED_*` status in `state.yaml` instead of silently narrowing completion.

Use the short pointer prompt from `docs/COPY_PASTE_AGENT_PROMPTS.md` when starting a Codex thread.
