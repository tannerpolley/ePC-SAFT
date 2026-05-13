# Assumptions

- The dependency gate is conservatively treated as `PREPARED_WAITING` until `/goal` rechecks GitHub issue #93, the prerequisite PRs, and the branch rebase gate.
- The live GoalBuddy board URL is `http://goalbuddy.localhost:41737/generic-reactive-lle/`.
- The issue scope note will be created during `/goal` execution, not during Goal Prep.
