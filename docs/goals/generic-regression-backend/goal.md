# Generic Regression Row Schema and Native Optimizer Backend

## Objective

Execute Task K: make regression generic around target rows, parameter maps, and native optimizer loops with real Jacobians.

## Original Request

Run Goal Prep for Task K, Generic regression row schema and native optimizer backend, using the local live GoalBuddy board in Codex.

## Intake Summary

- Input shape: `existing_plan`
- Audience: `epcsaft` package maintainers and future roadmap agents
- Authority: `requested`
- Proof type: `artifact`
- Completion proof: generic regression row compilation, parameter-map handling, native optimizer integration, and the requested validation/PR lifecycle are complete without application-specific public APIs or finite differences.
- Likely misfire: Adding application-specific regression APIs, Python-owned production objective loops, or finite-difference-backed derivatives instead of a generic native backend.
- Blind spots considered:
  - Dependencies B, C, and E are already merged into `origin/main`.
  - The branch must stay on `codex/generic-regression-backend`.
  - `auto_start_after_gate` is true, so the task should not stop at a waiting checkpoint.
  - No finite-difference derivative route is allowed.
  - No application-specific public APIs are allowed.
- Existing plan facts:
  - Task key: `K`.
  - Title: `Generic regression row schema and native optimizer backend`.
  - Assigned branch: `codex/generic-regression-backend`.
  - Dependencies: `B, C, E`.
  - Prompt file: `docs/roadmaps/agent_prompts/K_generic_regression_row_schema_and_native_optimizer_backend.md`.
  - Validation commands: `uv run python scripts/build_epcsaft.py --enable-ceres --enable-cppad`; `uv run python run_pytest.py tests/native/test_ceres_pure_regression.py tests/native/test_ceres_binary_regression.py tests/api/test_regression_api.py -q`; `uv run python scripts/validate_project.py quick`.
  - Do not add `fit_lithium_extraction_parameters` or `fit_mea_absorption` APIs.
  - Do not use finite difference or Python-owned production objective loops.

## Goal Kind

`existing_plan`

## Current Tranche

The current tranche is bounded watcher auto-run plus Task K implementation. The dependency gate passed on the live branch, so implementation starts immediately without asking.

## Non-Negotiable Constraints

- Stay on branch `codex/generic-regression-backend` before writing or executing Task K work.
- Do not create local `.worktrees/`.
- Do not add application-specific public APIs.
- Keep `epcsaft` general-purpose.
- No finite difference.
- Do not use Python-owned production objective loops.
- Ceres is preferred when it owns the native loop.
- Other native optimizers are allowed only with analytic, CppAD, or implicit derivatives.
- Use generic concepts such as `TargetDataset`, `RegressionProblem`, `EquilibriumProblem`, `ReactionSet`, `PhaseSpec`, and `ParameterSet`.

## Canonical Board

Machine truth lives at:

`docs/goals/generic-regression-backend/state.yaml`

If this charter and `state.yaml` disagree, `state.yaml` wins for task status, active task, receipts, verification freshness, and completion truth.

## Run Command

```text
/goal Follow docs/goals/generic-regression-backend/goal.md.
```

## PM Loop

On every `/goal` continuation:

1. Read this charter.
2. Read `state.yaml`.
3. Read `notes/dependency_gate.yaml`.
4. Re-check dependency B / C / E merge state before implementation.
5. Work only on the active board task.
6. If the gate fails or the branch cannot remain aligned with `origin/main`, write a compact receipt and stop with `BLOCKED_DEPENDENCY_OR_REBASE`.
7. If the gate passes, continue implementation automatically because `auto_start_after_gate` is true.
