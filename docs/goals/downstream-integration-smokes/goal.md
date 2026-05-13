# Downstream Integration Smoke Tests

## Charter

- Original request: Run Goal Prep for Task M, "Downstream integration smoke tests".
- Interpreted outcome: Prepare a GoalBuddy board for downstream smoke-test work, record the dependency gate, and keep implementation blocked until the prerequisite generic roadmap tasks are ready.
- Input shape: existing_plan
- Audience: epcsaft maintainers and downstream integration consumers
- Authority: requested
- Proof type: test, artifact
- Completion proof: Downstream projects can exercise generic `epcsaft` APIs without private workaround code, copied EOS implementation, or application-specific public APIs.
- Likely misfire: Turning this into downstream-specific package methods or trying to implement smoke-test behavior before the prerequisite generic APIs are available.
- Blind spots:
  - The issue scope for GitHub issue #96 must be recorded during `/goal` execution after the issue is read.
  - The dependency gate is expected to stay blocked until the prerequisite roadmap tasks are merged.
  - The package must remain application-neutral.
  - No finite-difference derivative paths are allowed.
- Existing plan facts:
  - Task key: M
  - Title: Downstream integration smoke tests
  - Assigned branch: codex/downstream-integration-smokes
  - Dependencies: F, G, I, J, K
  - Prompt file: docs/roadmaps/agent_prompts/M_downstream_integration_smoke_tests.md
  - Local live GoalBuddy board is required
  - Do not create local .worktrees/
  - Do not add application-specific public APIs
  - No finite difference
  - Dependency gate must pass before implementation starts

## Tranche

Prepare the dependency gate, keep the board watcher bounded, and do not start implementation until the prerequisite generic tasks are satisfied.

## Dependency Gate

- Gate file: `docs/goals/downstream-integration-smokes/notes/dependency_gate.yaml`
- Watcher script: `docs/goals/downstream-integration-smokes/notes/watch_dependency.ps1`
- Watcher mode: bounded
- Auto-start after gate: true
- Poll interval: 120 seconds
- Max wait: 480 minutes
- Branch update mode: ff-only
