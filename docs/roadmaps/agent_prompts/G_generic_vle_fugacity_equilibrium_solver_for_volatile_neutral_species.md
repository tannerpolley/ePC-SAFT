# Task G: Generic VLE/fugacity-equilibrium solver for volatile neutral species

## Goal Prep invocation

```text
[$goal-prep](C:\Users\Tanner\.codex\skills\goalbuddy\SKILL.md)

Initiate Goal Prep for this exact goal:

Task G: Generic VLE/fugacity-equilibrium solver for volatile neutral species

Use the local live GoalBuddy board in Codex.

Before Goal Prep writes files, run the branch bootstrap below. Do not write Goal Prep files on `main`.

```powershell
git fetch origin --prune
$current = (git branch --show-current).Trim()
if ($current -ne "codex/generic-vle-fugacity-equilibrium") {
    git ls-remote --exit-code --heads origin codex/generic-vle-fugacity-equilibrium
    if ($LASTEXITCODE -eq 0) {
        git switch codex/generic-vle-fugacity-equilibrium 2>$null
        if ($LASTEXITCODE -ne 0) {
            git switch -c codex/generic-vle-fugacity-equilibrium --track origin/codex/generic-vle-fugacity-equilibrium
        }
    } else {
        git switch -c codex/generic-vle-fugacity-equilibrium origin/main
    }
}
$current = (git branch --show-current).Trim()
if ($current -ne "codex/generic-vle-fugacity-equilibrium") {
    Write-Error "BRANCH_BOOTSTRAP_FAILED expected=codex/generic-vle-fugacity-equilibrium actual=$current"
    exit 1
}
git status --short
```

Do not start `/goal` execution yet.
Do not edit product/source files during Goal Prep.
Do not create local `.worktrees/`.
Do not ask repeated confirmation questions.

Create:
- docs/goals/<slug>/goal.md
- docs/goals/<slug>/state.yaml
- docs/goals/<slug>/notes/
- docs/goals/<slug>/notes/dependency_gate.yaml
- docs/goals/<slug>/notes/watch_dependency.ps1

After branch bootstrap, current branch must be:

codex/generic-vle-fugacity-equilibrium

Before planning or editing, read:

docs/roadmaps/general_reactive_electrolyte_equilibrium_readiness.md
docs/roadmaps/agent_dependency_plan.md
docs/roadmaps/agent_prompts/index.yaml

Dependencies:
- C

After Goal Prep, use bounded watcher auto-run:
1. Create or update `docs/goals/<slug>/notes/dependency_gate.yaml`.
2. Create or update `docs/goals/<slug>/notes/watch_dependency.ps1` during Goal Prep, even if dependencies appear satisfied, so the GoalBuddy board and file-based watcher state are complete before `/goal` execution.
3. Set `watcher_mode: bounded`, `auto_start_after_gate: true`, `poll_interval_seconds: 120`, and `max_wait_minutes: 480` in local GoalBuddy/dependency files.
4. If dependencies are satisfied, start implementation immediately without asking.
5. If dependencies are not satisfied, run the bounded watcher, poll every 120 seconds, stop after 480 minutes with `PREPARED_WAITING` if still blocked, and start implementation without asking when the gate passes.

Do not stop at `PREPARED_READY` merely because `auto_start_after_gate` is false. That manual checkpoint mode is retired unless the user explicitly asks for it in a future issue.
Do not ask the user what to do next.
Do not ask repeated confirmation questions.
```

## Non-interaction rule

Do not ask the user to choose among options.
Do not ask repeated confirmation questions.
Use the defaults written in this prompt.
If a required value is missing, choose the safest conservative default and record it in `docs/goals/<slug>/notes/assumptions.md`.
If the missing value prevents safe work, stop with status `BLOCKED_MISSING_INPUT`.

## Package-scope rule

The `epcsaft` package must stay general-purpose.

Do not add public APIs named after MEA, lithium extraction, absorption columns, extraction efficiency, distribution coefficient, selectivity, or any application-specific workflow.

Use generic concepts:

```python
equilibrium(...)
regress_parameters(...)
ReactionSet(...)
EquilibriumProblem(...)
RegressionProblem(...)
TargetDataset(...)
PhaseSpec(...)
ParameterSet(...)
```

## Derivative rule

No finite difference.

Allowed derivative backends:

```text
analytic
cppad
analytic_implicit
cppad_implicit
legacy_eigen_forward only for validated legacy/local paths
backend_unavailable only for explicitly out-of-scope workflows
```

CppAD is the default for explicit algebraic derivatives.
Solved states use analytic_implicit or cppad_implicit sensitivities.
Do not tape iterative solver loops as production derivatives.

## GitHub issue

#90

## Task summary

Implement or harden generic fugacity-equilibrium routes for volatile neutral species.

## Scope

- direct volatile partial pressure from liquid fugacity where valid
- bubble pressure
- bubble temperature
- dew pressure
- dew temperature
- TP flash
- diagnostics report route used

## Do not do

- do not create MEA bubble-pressure-specific public APIs
- do not assume ions distribute to vapor unless problem explicitly models it

## Dependency gate before implementation

Before implementation, read the GitHub issue and run the dependency gate again:

```powershell
git fetch origin --prune
gh issue view 90 --repo tannerpolley/ePC-SAFT --json number,title,state,body,url
git branch --show-current
git status --short
git rebase origin/main
```

Read the GitHub issue before implementation and include its current scope in `docs/goals/<slug>/notes/issue_scope.md`. Verify prerequisite issues are closed and PRs are merged. If any dependency is missing or the rebase conflicts, stop with status `BLOCKED_DEPENDENCY_OR_REBASE`.

## Implementation, PR, self-review, and merge automation

This prompt pre-authorizes the agent to do the full issue lifecycle without asking for another yes when the gates below pass.

Required sequence:

1. Inspect GitHub issue #90 and implement only that assigned roadmap task on the assigned branch.
2. Run the task-specific validation commands and the repo-level validation named in this prompt.
3. Run `git diff --check`.
4. Rebase or fast-forward against `origin/main`, then review the branch against `origin/main`:
   - `git fetch origin --prune`
   - `git status --short`
   - `git diff --stat origin/main...HEAD`
   - inspect the changed files and confirm they match this task scope.
5. Open a focused draft PR with the GitHub CLI if no PR exists for the branch, and include Closes #90 in the PR body.
6. Self-review the PR against `origin/main` before marking it ready:
   - no application-specific public APIs
   - no finite-difference derivative route
   - no unrelated files
   - all in-scope items are classified as `implemented`, `already_supported_with_tests`, `blocker_requires_followup`, or `out_of_scope_by_roadmap`
   - no silent narrowing of the prompt scope
7. If self-review passes, mark the PR ready for review.
8. Wait for GitHub checks to finish.
9. If checks pass, the PR is mergeable, and the final GoalBuddy audit says `full_outcome_complete: true`, merge the PR without asking for additional user confirmation.
10. After merge, confirm issue #90 is closed; if it is still open, close it with a comment that names the merged PR and merge commit.
11. Delete the remote branch and local branch used by this task without asking for additional user confirmation.
12. Record the issue URL, issue close status, PR URL, merge commit, validation commands, remote branch deletion, local branch deletion, and any cleanup blocker in `state.yaml`.

Suggested GitHub CLI commands:

```powershell
gh pr list --head <branch> --state open --json number,url,isDraft,mergeable,statusCheckRollup
gh pr create --draft --base main --head <branch> --title "<task title>" --body "<summary, tests, limitations, dependencies>

Closes #90"
gh pr diff <number> --name-only
gh pr ready <number>
gh pr checks <number> --watch --fail-fast
gh pr merge <number> --merge --delete-branch
gh issue view 90 --repo tannerpolley/ePC-SAFT --json state,url
gh issue close 90 --repo tannerpolley/ePC-SAFT --comment "Completed by merged PR <number> at <merge commit>."  # only if still open
git fetch origin main --prune
git switch --detach origin/main
git branch -d codex/generic-vle-fugacity-equilibrium
```

Do not merge if any of these are true:

- GitHub issue scope was not inspected or recorded
- required validation failed or was not run
- GitHub checks failed, are pending, or are unavailable
- PR is not mergeable
- branch is not current with `origin/main`
- final GoalBuddy audit is missing or says `full_outcome_complete: false`
- task scope was narrowed silently
- changed files are outside the assigned task scope
- review finds application-specific APIs or finite differences
- credentials, network, or GitHub policy prevent the merge
- local branch cleanup would discard uncommitted work or delete the wrong branch

If blocked, do not ask repeated confirmation questions. Stop with a precise status such as `BLOCKED_ISSUE_SCOPE_UNREAD`, `BLOCKED_ISSUE_CLOSE_FAILED`, `BLOCKED_CHECKS_FAILED`, `BLOCKED_REBASE_CONFLICT`, `BLOCKED_SCOPE_GAP`, `BLOCKED_GITHUB_POLICY`, `BLOCKED_MERGE_CONFLICT`, `BLOCKED_REMOTE_BRANCH_DELETE`, or `BLOCKED_LOCAL_BRANCH_DELETE`, and write the exact blocker and next command to `state.yaml`.

## Validation

```powershell
uv run python run_pytest.py tests/equilibrium/test_derivative_policy.py tests/native/test_cppad_bubble_derivatives.py -q
```
```powershell
uv run python scripts/validate_project.py quick
```

## Branch cleanup after merge

After the PR is merged, delete both the remote and local task branch without asking again.

Use this cleanup sequence:

```powershell
git fetch origin main --prune
git switch --detach origin/main
git branch -d codex/generic-vle-fugacity-equilibrium
git push origin --delete codex/generic-vle-fugacity-equilibrium  # only if gh pr merge --delete-branch did not already delete it
git -C C:\Users\Tanner\Documents\git\ePC-SAFT fetch origin main --prune
git -C C:\Users\Tanner\Documents\git\ePC-SAFT pull --ff-only origin main
```

If local branch deletion fails because the branch is still checked out or has unmerged local work, stop with `BLOCKED_LOCAL_BRANCH_DELETE` and record the exact checkout path, branch name, and safe next command. Do not force-delete a branch with unmerged work.

Report:

```text
issue URL
issue closed yes/no
PR URL
merge commit
remote branch deleted yes/no
local branch deleted yes/no
local main updated yes/no
next unblocked issues
```
