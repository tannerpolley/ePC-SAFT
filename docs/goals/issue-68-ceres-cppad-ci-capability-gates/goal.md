# Issue #68 Ceres/CppAD CI And Capability Gates

## Objective

Complete GitHub issue #68 only: add the requested Ceres/CppAD CI and capability gates, literature-backed regression validation using only existing repo assets, required validation evidence, and a focused PR linked to and closing #68.

## Original Request

Initiate Goal Prep for Issue #68 Ceres/CppAD CI and capability gates using the local live GoalBuddy board. Do not start `/goal` execution yet. Preserve the execution constraints, prerequisites, validation, PR, merge, branch cleanup, and issue closure requirements.

## Intake Summary

- Input shape: `existing_plan`
- Audience: ePC-SAFT maintainers, CI users, and downstream modeling agents that need reliable Ceres/CppAD capability detection and regression validation.
- Authority: `requested`
- Proof type: `test | artifact | review`
- Completion proof: issue #68 is implemented on the correct branch after issue #67 is merged and local main is updated; the validation commands listed in issue #68 pass or have explicit blocker evidence; a focused PR with `Closes #68` passes checks; the PR is squash-merged; the branch is deleted; root checkout local `main` is fast-forwarded to match `origin/main`; issue #68 is closed and labeled when labels are available.
- Likely misfire: starting implementation from the wrong base, continuing PR #56, skipping the #67 merge gate, fetching new internet data, adding finite-difference derivative tests, or validating vapor Born derivatives despite the user forbidding those paths.
- Blind spots considered: issue #68 may have a detailed addendum that changes task order; #67/its PR may still block the correct base; validation command names must come from issue #68, not memory; GitHub permissions may block merge, branch deletion, labeling, or issue closure.
- Existing plan facts: preserve every user-supplied execution constraint in this charter and `state.yaml`; do not use PR #56 as a base branch; use PR #56 only as optional reference; stop after issue #68.

## Goal Kind

`existing_plan`

## Current Tranche

Prepare for and then execute issue #68 end to end after the prerequisite #67 gate is satisfied. The first safe slice for `/goal` is a PM prerequisite gate: confirm #67 and its attached PR are merged into `origin/main`, confirm the root checkout local `main` has fast-forwarded to `origin/main`, ensure the worktree branch is `codex/issue-68-ceres-cppad-ci-capabilities`, and rebase the branch onto current `origin/main` or local `main` before reading and implementing issue #68.

## Non-Negotiable Constraints

- Goal Prep must not start `/goal` execution.
- Work on GitHub issue #68 only.
- Do not plan implementation details beyond what is needed to create the GoalBuddy board during Goal Prep.
- Do not edit product/source files during Goal Prep.
- You are already in a Codex app worktree created from updated local `main`; local `main` should have been fast-forwarded from `origin/main` immediately before this worktree was created.
- Then checkout into branch `codex/issue-68-ceres-cppad-ci-capabilities`.
- Before starting issue #68 work, confirm issue #67 and its attached PR have merged to `origin/main`, then confirm that merge has been pulled to local root `main`.
- Set up an automation watcher script for the #67 prerequisite gate if #67 is not yet merged; the watcher must wait until #67/its PR are confirmed merged to `origin/main`, confirm local root `main` was pulled, then rebase this branch onto current `origin/main` or local `main`.
- If issue #67 is not merged when checked, stop and report that #67 is blocked.
- Read issue #68 completely, including the Detailed implementation addendum, before planning or editing.
- Do not use PR #56 as a base branch.
- Do not continue or merge PR #56.
- Use PR #56 only as reference if needed.
- Add literature-backed regression validation using only existing repo assets: `docs/papers`, `analyses`, `data/reference`, Figiel 2025, and MIAC outputs.
- Do not fetch validation data from the internet.
- No finite-difference derivatives.
- Do not test vapor Born derivatives.
- Run the validation commands listed in issue #68.
- Open a focused PR linked to issue #68.
- The PR body must include `Closes #68`.
- If validation and GitHub checks pass, squash-merge the PR, delete the branch, fast-forward the root checkout local `main`, confirm issue #68 closed, close it manually if needed with a final comment, apply a completion label if labels are available, and report the PR URL, merge commit, validation summary, label status, and final issue state.
- If merge or issue closure is blocked, leave the PR ready and report the PR URL plus exactly what blocks merge or issue closure.
- Stop after this issue.

## Required Post-Merge Root Checkout Fast-Forward

Before continuing or reporting completion after a successful merge, run from any shell:

```powershell
git -C C:\Users\Tanner\Documents\git\ePC-SAFT pull --ff-only origin main
git -C C:\Users\Tanner\Documents\git\ePC-SAFT rev-parse --short=12 main
git -C C:\Users\Tanner\Documents\git\ePC-SAFT rev-parse --short=12 origin/main
```

The `main` and `origin/main` hashes must match and include the PR merge commit. If this fails, report that future worktrees from `main` are blocked until local `main` is updated.

## Stop Rule

Stop only when a final audit proves the full issue #68 outcome is complete, or when a user-specified prerequisite or GitHub permission gate blocks safe continuation.

Do not continue into issue #68 implementation if issue #67 is not confirmed merged into `origin/main` and local root `main` is not confirmed updated.

Do not stop after reading issue #68, planning, or choosing a Worker task if the prerequisite gate is satisfied and safe implementation tasks remain.

## Canonical Board

Machine truth lives at:

`docs/goals/issue-68-ceres-cppad-ci-capability-gates/state.yaml`

If this charter and `state.yaml` disagree, `state.yaml` wins for task status, active task, receipts, verification freshness, and completion truth.

## Run Command

```text
/goal Follow docs/goals/issue-68-ceres-cppad-ci-capability-gates/goal.md.
```

## PM Loop

On every `/goal` continuation:

1. Read this charter.
2. Read `state.yaml`.
3. Run the bundled GoalBuddy update checker when available and mention a newer version without blocking.
4. Work only on the active board task.
5. Enforce the #67 prerequisite gate before any issue #68 planning or implementation.
6. Read issue #68 completely, including the Detailed implementation addendum, before planning or editing.
7. Preserve the user constraints around PR #56, validation assets, no internet data, no finite differences, and no vapor Born derivative tests.
8. Assign Scout, Judge, Worker, or PM according to the task.
9. Write a compact task receipt and update the board.
10. Treat a slice audit as a checkpoint, not completion, unless it explicitly proves every issue #68 completion requirement.
11. Finish only with a Judge/PM audit receipt that maps receipts and validation back to issue #68 and records `full_outcome_complete: true`.
