# Assumptions

- Goal Prep is allowed to stop in PREPARED_WAITING because dependency F is still blocked.
- The issue scope for GitHub issue #96 will be recorded during `/goal` execution after the issue is read.
- The local GoalBuddy board should use the installed user-level Scout, Worker, and Judge agent configs.
- `watcher_mode: bounded`, `auto_start_after_gate: true`, `poll_interval_seconds: 120`, and `max_wait_minutes: 480` are the required defaults.
