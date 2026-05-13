# Issue #66 Ceres Pure And Binary Regression Backend

## Objective

Prepare and execute GitHub issue #66: implement native Ceres-backed regression for pure and binary ePC-SAFT parameter fitting, with Ceres owning the optimizer loop and Python limited to validation and serialization.

## Original Request

Initiate Goal Prep for issue #66, using the local live GoalBuddy board in Codex, creating only the GoalBuddy setup artifacts and preserving the execution constraints for the later `/goal` run.

## Intake Summary

- Input shape: `existing_plan`
- Audience: ePC-SAFT maintainers and downstream users who need production native parameter regression.
- Authority: `requested`
- Proof type: `test`
- Completion proof: issue #66 is implemented on a focused PR with `Closes #66`, the issue-listed validation commands and GitHub checks pass, the PR is squash-merged, the branch is deleted, and the root checkout local `main` is fast-forwarded from `origin/main`.
- Likely misfire: GoalBuddy could implement a Python-owned fitting loop, finite-difference fallback, or broad reactive pressure/speciation fit and still appear to have "regression" working.
- Blind spots considered: derivative coverage may be incomplete after #62; Ceres rows must capability-gate unsupported derivative families; association and solved-state sensitivities need implicit derivatives or truthful `backend_unavailable`; PR #56 is reference-only and must not become the base.
- Existing plan facts: issue #66 and its detailed addenda define files to inspect/create, required Ceres behavior, target rows, result fields, tests, validation commands, stop conditions, derivative policy, and post-#62 capability gating.

## Goal Kind

`existing_plan`

## Current Tranche

Carry issue #66 from verified setup through a focused implementation PR and post-merge cleanup. The first pure-neutral Ceres foundation slice is implemented locally. After issue comment `https://github.com/tannerpolley/ePC-SAFT/issues/66#issuecomment-4436520857`, the active path is to continue on `codex/issue-66-ceres-regression-backend-2`, keep the pure-neutral foundation logically isolated, and implement the native binary `k_ij` derivative subsystem next. Do not close #66 until pure + binary Ceres support is validated, or until #66 is explicitly split/narrowed.

## Non-Negotiable Constraints

- Work only in the existing Codex app worktree; do not create a local git worktree under `.worktrees/`.
- Current branch must be `codex/issue-66-ceres-regression-backend-2`.
- Local main is assumed to have been fast-forwarded from `origin/main` before this app worktree was created; if evidence contradicts that, stop and report it.
- Work on GitHub issue #66 only.
- Confirm issues #59, #60, #61, #62, and #65 remain merged into `main` before implementation.
- Read issue #66 completely, including the Detailed implementation addendum and later derivative/capability-gating addenda, before planning implementation or editing product/source files.
- Do not use PR #56 as a base branch. Do not continue or merge PR #56. Use PR #56 only as reference if needed.
- Implement Ceres-backed regression for pure and binary ePC-SAFT parameter fitting.
- Ceres must own the optimizer loop.
- Python validates and serializes only.
- No finite-difference fallback. For the resumed binary derivative work, no finite differences whatsoever are allowed in the Ceres production derivative path.
- No missing Jacobian columns.
- No Python-owned production objective loop.
- No optimistic production capability flags.
- Do not start with fully coupled reactive pressure/speciation fits.
- Unsupported derivative families must return `backend_unavailable` and must not be advertised as production-supported rows.
- Continue on the current branch rather than opening the optional partial foundation PR unless the board is explicitly redirected.
- Keep the existing pure-neutral Ceres foundation logically separate from binary derivative-subsystem edits in receipts, tests, and PR description.
- Run the validation commands listed in issue #66:

```powershell
uv run python scripts/build_epcsaft.py --clean --enable-ceres --enable-cppad
uv run python run_pytest.py tests/native/test_ceres_pure_regression.py tests/native/test_ceres_binary_regression.py tests/native/test_ceres_liquid_electrolyte_regression.py -q
uv run python scripts/validate_project.py quick
```

- Open a focused PR linked to issue #66.
- The PR body must include `Closes #66`.
- If validation and GitHub checks pass, squash-merge the PR, delete the branch, then fast-forward the root checkout's local `main` from `origin/main` before reporting completion.

## Verified Setup Facts From Goal Prep

- Current branch was verified as `codex/issue-66-ceres-regression-backend-2`.
- GitHub CLI authentication was verified for account `tannerpolley`.
- Issues #59, #60, #61, #62, and #65 were closed and had merged PRs into `main`.
- Issue #66 was open and its body/addenda were read during Goal Prep.
- GoalBuddy Scout/Judge/Worker agent configs were found installed under the user Codex agents directory.

## Stop Rule

Stop only when a final audit proves the full original outcome is complete.

Do not stop after planning, discovery, or Judge selection if a safe Worker task can be activated.

Do not stop after a single verified Worker slice while issue #66 still has required local follow-up work. After each slice audit, advance the board to the next highest-leverage safe Worker task and continue.

Stop and write a handoff if Ceres integration depends on unfinished APIs from issues #59-#62, if finite-difference fallback appears necessary, if optimizer ownership drifts back into Python, if the branch changes away from `codex/issue-66-ceres-regression-backend-2`, or if implementation requires PR #56 as a base.

## Canonical Board

Machine truth lives at:

`docs/goals/issue-66-ceres-regression-backend/state.yaml`

If this charter and `state.yaml` disagree, `state.yaml` wins for task status, active task, receipts, verification freshness, and completion truth.

## Run Command

```text
/goal Follow docs/goals/issue-66-ceres-regression-backend/goal.md.
```

## PM Loop

On every `/goal` continuation:

1. Read this charter.
2. Read `state.yaml`.
3. Run the bundled GoalBuddy update checker when available and mention a newer version without blocking.
4. Re-check the intake, current branch, prerequisite issue merge status, original request, proof, blind spots, existing plan facts, and likely misfire.
5. Work only on the active board task.
6. Assign Scout, Judge, Worker, or PM according to the task.
7. Write a compact task receipt.
8. Update the board.
9. If Judge selected a safe Worker task with `allowed_files`, `verify`, and `stop_if`, activate it and continue unless blocked.
10. Treat a slice audit as a checkpoint, not completion, unless it explicitly proves the full original issue #66 outcome is complete.
11. Finish only with a Judge/PM audit receipt that maps receipts and verification back to issue #66 and records `full_outcome_complete: true`.
