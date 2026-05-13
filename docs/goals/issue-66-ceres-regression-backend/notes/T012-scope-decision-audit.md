# T012 Scope Decision Audit

## Objective Restated

Complete GitHub issue #66 by implementing native Ceres-backed pure and binary ePC-SAFT regression, validating it locally, publishing a focused PR with `Closes #66`, waiting for green checks, squash-merging, deleting the branch, and fast-forwarding the root checkout local `main`.

## Prompt-To-Artifact Checklist

| Requirement | Evidence | Status |
| --- | --- | --- |
| Work on branch `codex/issue-66-ceres-regression-backend-2` | `git branch --show-current` returned `codex/issue-66-ceres-regression-backend-2`. | Satisfied |
| Work on issue #66 only | Current GoalBuddy artifacts and implementation receipts are scoped to issue #66. | Satisfied so far |
| Do not create `.worktrees/` | No task receipt records creating a worktree; work is in the Codex app worktree. | Satisfied so far |
| Confirm prerequisite issues #59, #60, #61, #62, #65 merged | Goal prep setup evidence records closed issues with merged PRs #73, #75, #76, #77, #74. | Satisfied during setup |
| Read issue #66 completely including addenda | Goal prep and T001 receipts record complete issue body/addenda read; latest `gh issue view 66` confirms the same required scope. | Satisfied |
| Ceres owns optimizer loop for supported production rows | Pure-neutral native Ceres path uses `ceres::CostFunction` and `ceres::Solve` in `src/epcsaft/native/epcsaft_regression.cpp`. | Satisfied for pure-neutral slice |
| Python validates and serializes only | Pure-neutral path calls native `_fit_pure_neutral_native_ceres`; test asserts `python_objective_used is False`. | Satisfied for pure-neutral slice |
| No finite-difference fallback for Ceres production rows | Pure-neutral Ceres uses existing native residual/Jacobian evaluator; binary/liquid Ceres requests return `backend_unavailable`. | Satisfied for current supported/gated rows |
| Required pure rows: density, vapor pressure | Pure-neutral Ceres test covers density and vapor-pressure residuals. | Satisfied for pure-neutral slice |
| Required relative permittivity where needed | Current issue slice does not implement production relative-permittivity Ceres rows beyond gating. | Incomplete |
| Required binary `k_ij` target rows | `tests/native/test_ceres_binary_regression.py` asserts explicit Ceres binary request returns `backend_unavailable`. | Incomplete |
| Required liquid-electrolyte `d_born`/`f_solv` rows | `tests/native/test_ceres_liquid_electrolyte_regression.py` asserts Ceres liquid-electrolyte fitting returns `backend_unavailable`. | Incomplete |
| Association parameters | T011 records association parameter production support remains blocked without implicit sensitivities. | Incomplete |
| Result contract fields | `FitResult` and native result conversion include optimizer/derivative/status/iteration/objective/residual/gradient/step/map/diagnostics fields for current paths. | Partially satisfied |
| Issue validation command: Ceres build | Previously passed: `uv run python scripts/build_epcsaft.py --clean --enable-ceres --enable-cppad`. | Satisfied for current diff |
| Issue validation command: Ceres tests | Previously passed with current tests, but binary/liquid tests cover truthful gates, not production support. | Weak for full issue |
| Issue validation command: quick validation | Previously passed: `uv run python scripts/validate_project.py quick`. | Satisfied for current diff |
| Open focused PR linked to #66 with `Closes #66` | `gh pr list --head codex/issue-66-ceres-regression-backend-2` returned `[]`. | Missing |
| GitHub checks green | No PR exists. | Missing |
| Squash-merge PR | No PR exists. | Missing |
| Delete branch | Branch still active. | Missing |
| Fast-forward root checkout local `main` | No merged PR exists to fast-forward from. | Missing |

## Current Scope Decision

Full issue #66 is not complete. The current diff implements and validates a useful pure-neutral Ceres slice plus truthful binary/liquid capability gates, but issue #66 still requires production support for binary `k_ij`, liquid-electrolyte `d_born`/`f_solv`, and association-parameter handling or explicitly accepted issue narrowing.

T011 judged the next production step as a broad native derivative-subsystem implementation, not a bounded wrapper task. T012 therefore remains blocked until one of these scope decisions is provided:

1. Authorize the larger native derivative-subsystem implementation for full issue #66.
2. Split or narrow issue #66 so the current pure-neutral plus capability-gate slice can be published without `Closes #66`.
3. Pause the goal with this board state as the handoff.

