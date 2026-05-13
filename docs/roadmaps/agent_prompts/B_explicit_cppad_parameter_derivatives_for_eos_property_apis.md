# Task B: Explicit CppAD parameter derivatives for EOS/property APIs

## Goal Prep invocation

```text
[$goal-prep](C:\Users\Tanner\.codex\skills\goalbuddy\SKILL.md)

Initiate Goal Prep for this exact goal:

Task B: Explicit CppAD parameter derivatives for EOS/property APIs

Use the local live GoalBuddy board in Codex.

Before Goal Prep writes files, run the branch bootstrap below. Do not write Goal Prep files on `main`.

```powershell
git fetch origin --prune
$current = (git branch --show-current).Trim()
if ($current -ne "codex/cppad-explicit-parameter-derivatives") {
    git ls-remote --exit-code --heads origin codex/cppad-explicit-parameter-derivatives
    if ($LASTEXITCODE -eq 0) {
        git switch codex/cppad-explicit-parameter-derivatives 2>$null
        if ($LASTEXITCODE -ne 0) {
            git switch -c codex/cppad-explicit-parameter-derivatives --track origin/codex/cppad-explicit-parameter-derivatives
        }
    } else {
        git switch -c codex/cppad-explicit-parameter-derivatives origin/main
    }
}
$current = (git branch --show-current).Trim()
if ($current -ne "codex/cppad-explicit-parameter-derivatives") {
    Write-Error "BRANCH_BOOTSTRAP_FAILED expected=codex/cppad-explicit-parameter-derivatives actual=$current"
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

codex/cppad-explicit-parameter-derivatives

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

Implement or verify CppAD/analytic derivatives for explicit algebraic parameter effects used by properties and regression.

## Scope

- m, sigma, epsilon explicit effects
- k_ij, l_ij, k_hb_ij direct effects where algebraic
- d_born, f_solv, dielectric parameters
- pressure/fugacity/activity/chemical-potential parameter derivatives where explicit
- regression APIs request derivatives by parameter name

## Do not do

- do not implement implicit solved-state sensitivities in this issue
- do not tape solver loops

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
uv run python scripts/build_epcsaft.py --enable-cppad
```
```powershell
uv run python run_pytest.py tests/native/test_cppad_pressure_derivatives.py tests/native/test_cppad_fugacity_derivatives.py tests/native/test_cppad_activity_derivatives.py tests/native/test_cppad_relative_permittivity_derivatives.py -q
```
```powershell
uv run python scripts/validate_project.py quick
```

## Branch cleanup after merge

After the PR is merged, run or report the equivalent:

```powershell
git -C C:\Users\Tanner\Documents\git\ePC-SAFT fetch origin main --prune
git -C C:\Users\Tanner\Documents\git\ePC-SAFT pull --ff-only origin main
git push origin --delete codex/cppad-explicit-parameter-derivatives
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
