# Task J: Generic reactive LLE and chemical phase equilibrium

## Goal Prep invocation

```text
[$goal-prep](C:\Users\Tanner\.codex\skills\goalbuddy\SKILL.md)

Initiate Goal Prep for this exact goal:

Task J: Generic reactive LLE and chemical phase equilibrium

Use the local live GoalBuddy board in Codex.

Before Goal Prep writes files, run the branch bootstrap below. Do not write Goal Prep files on `main`.

```powershell
git fetch origin --prune
$current = (git branch --show-current).Trim()
if ($current -ne "codex/generic-reactive-lle") {
    git ls-remote --exit-code --heads origin codex/generic-reactive-lle
    if ($LASTEXITCODE -eq 0) {
        git switch codex/generic-reactive-lle 2>$null
        if ($LASTEXITCODE -ne 0) {
            git switch -c codex/generic-reactive-lle --track origin/codex/generic-reactive-lle
        }
    } else {
        git switch -c codex/generic-reactive-lle origin/main
    }
}
$current = (git branch --show-current).Trim()
if ($current -ne "codex/generic-reactive-lle") {
    Write-Error "BRANCH_BOOTSTRAP_FAILED expected=codex/generic-reactive-lle actual=$current"
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
- docs/goals/<slug>/dependency_gate.yaml

After branch bootstrap, current branch must be:

codex/generic-reactive-lle

Before planning or editing, read:

docs/roadmaps/general_reactive_electrolyte_equilibrium_readiness.md
docs/roadmaps/agent_dependency_plan.md
docs/roadmaps/agent_prompts/index.yaml

Dependencies:
- I
- D
- F

After Goal Prep, do exactly one of:
1. If dependencies are not satisfied: write PREPARED_WAITING and stop.
2. If dependencies are satisfied and auto_start_after_gate is true: begin following goal.md directly without asking.
3. If dependencies are satisfied and auto_start_after_gate is false: write PREPARED_READY and stop.

Do not ask the user what to do next.
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

Combine reaction/speciation and LLE in a generic chemical-phase-equilibrium framework.

## Scope

- reaction coordinates
- element balances
- phase split
- fugacity equality
- reaction equilibrium
- nonnegativity
- Ascani 2023 benchmark attempt

## Do not do

- do not create extraction-specific or MEA-specific public APIs
- do not start with a full constrained NLP unless the issue explicitly decides that path

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
- Do not modify PR #56.
- Do not base work on PR #56.

## Validation

```powershell
uv run python run_pytest.py tests/equilibrium/test_reactive_lle.py tests/api/test_reactive_staged_workflow_contract.py -q
```
```powershell
uv run python scripts/validate_project.py quick
```

## Branch cleanup after merge

After the PR is merged, run or report the equivalent:

```powershell
git -C C:\Users\Tanner\Documents\git\ePC-SAFT fetch origin main --prune
git -C C:\Users\Tanner\Documents\git\ePC-SAFT pull --ff-only origin main
git push origin --delete codex/generic-reactive-lle
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
