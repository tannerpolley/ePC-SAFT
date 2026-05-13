# Task E: Generic target-row and dataset schema

## Goal Prep invocation

```text
[$goal-prep](C:\Users\Tanner\.codex\skills\goalbuddy\SKILL.md)

Initiate Goal Prep for this exact goal:

Task E: Generic target-row and dataset schema

Use the local live GoalBuddy board in Codex.

Before Goal Prep writes files, run the branch bootstrap below. Do not write Goal Prep files on `main`.

```powershell
git fetch origin --prune
$current = (git branch --show-current).Trim()
if ($current -ne "codex/generic-target-row-schema") {
    git ls-remote --exit-code --heads origin codex/generic-target-row-schema
    if ($LASTEXITCODE -eq 0) {
        git switch codex/generic-target-row-schema 2>$null
        if ($LASTEXITCODE -ne 0) {
            git switch -c codex/generic-target-row-schema --track origin/codex/generic-target-row-schema
        }
    } else {
        git switch -c codex/generic-target-row-schema origin/main
    }
}
$current = (git branch --show-current).Trim()
if ($current -ne "codex/generic-target-row-schema") {
    Write-Error "BRANCH_BOOTSTRAP_FAILED expected=codex/generic-target-row-schema actual=$current"
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

After branch bootstrap, current branch must be:

codex/generic-target-row-schema

Before planning or editing, read:

docs/roadmaps/general_reactive_electrolyte_equilibrium_readiness.md
docs/roadmaps/agent_dependency_plan.md
docs/roadmaps/agent_prompts/index.yaml

Dependencies:
- A

After Goal Prep, use bounded watcher auto-run:
1. Create or update `docs/goals/<slug>/notes/dependency_gate.yaml`.
2. Create or update `docs/goals/<slug>/notes/watch_dependency.ps1` when dependencies are not already satisfied.
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

## Task summary

Create application-neutral regression and validation target-row schemas for future generic regression/equilibrium workflows.

## Scope

- pure density
- pure vapor pressure
- P-rho-T
- binary VLE
- binary LLE
- osmotic coefficient
- MIAC
- relative permittivity
- activity/fugacity
- speciation
- VLE partial pressure
- LLE phase composition
- regularization

## Do not do

- do not add lithium-specific or MEA-specific public APIs
- do not implement optimizer internals here

## Dependency gate before implementation

Before implementation, run the dependency gate again:

```powershell
git fetch origin --prune
git branch --show-current
git status --short
git rebase origin/main
```

Verify prerequisite issues are closed and PRs are merged. If any dependency is missing or the rebase conflicts, stop with status `BLOCKED_DEPENDENCY_OR_REBASE`.

## Implementation expectations

- One issue per branch.
- One PR for this issue.
- Open a focused PR.
- PR body includes issue link, summary, tests, limitations, and next dependencies.
- Do not close unrelated issues.

## Validation

```powershell
uv run python run_pytest.py tests/api/test_regression_problem_schema.py tests/api/test_regression_api.py -q
```
```powershell
uv run python scripts/validate_project.py docs
```

## Branch cleanup after merge

After the PR is merged, run or report the equivalent:

```powershell
git -C C:\Users\Tanner\Documents\git\ePC-SAFT fetch origin main --prune
git -C C:\Users\Tanner\Documents\git\ePC-SAFT pull --ff-only origin main
git push origin --delete codex/generic-target-row-schema
```

Report:

```text
issue URL
PR URL
merge commit
branch deleted yes/no
local main updated yes/no
next unblocked issues
```
