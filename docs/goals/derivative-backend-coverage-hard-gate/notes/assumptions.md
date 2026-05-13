# Assumptions

- The task has no issue URL because `issue: null` is recorded for task A in `docs/roadmaps/agent_prompts/index.yaml`; the safest default is to leave issue linkage unknown until implementation/PR work discovers or creates the correct artifact.
- The goal slug is `derivative-backend-coverage-hard-gate`, derived from the task title.
- The local live board is the selected visual board because the prompt explicitly requests it.
- Dedicated GoalBuddy agents are treated as installed because matching user-level agent configs were found under `C:\Users\Tanner\.codex\agents`.
