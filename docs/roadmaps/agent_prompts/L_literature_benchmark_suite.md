# Task L: Literature benchmark suite

## Goal Prep invocation

```text
[$goal-prep](C:\Users\Tanner\.codex\skills\goalbuddy\SKILL.md)

Initiate Goal Prep for this exact goal:

Task L: Literature benchmark suite

Use the local live GoalBuddy board in Codex.

Before Goal Prep writes files, run the branch bootstrap below. Do not write Goal Prep files on `main`.

```powershell
git fetch origin --prune
$current = (git branch --show-current).Trim()
if ($current -ne "codex/literature-benchmark-suite") {
    git ls-remote --exit-code --heads origin codex/literature-benchmark-suite
    if ($LASTEXITCODE -eq 0) {
        git switch codex/literature-benchmark-suite 2>$null
        if ($LASTEXITCODE -ne 0) {
            git switch -c codex/literature-benchmark-suite --track origin/codex/literature-benchmark-suite
        }
    } else {
        git switch -c codex/literature-benchmark-suite origin/main
    }
}
$current = (git branch --show-current).Trim()
if ($current -ne "codex/literature-benchmark-suite") {
    Write-Error "BRANCH_BOOTSTRAP_FAILED expected=codex/literature-benchmark-suite actual=$current"
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

codex/literature-benchmark-suite

Before planning or editing, read:

docs/roadmaps/general_reactive_electrolyte_equilibrium_readiness.md
docs/roadmaps/agent_dependency_plan.md
docs/roadmaps/agent_prompts/index.yaml

Dependencies:
- none

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

Inventory and then build generic literature benchmarks for package-level confidence.

## Scope

- fixture inventory can start early
- implementation waits on relevant solver/regression issues
- MEA simple workflow benchmark
- MDEA ePC-SAFT benchmark
- Figiel 2025 SSM+DS Born benchmark
- Held 2014 revised ePC-SAFT benchmark
- non-electrolyte LLE benchmark
- Ascani 2022 electrolyte LLE benchmark
- Ascani 2023 reactive LLE benchmark
- Khudaida salting-out LLE benchmark
- Hubach/Yu lithium-related equilibrium benchmark

## Do not do

- do not implement benchmark tests that depend on missing APIs
- do not require downstream repo access for package-level benchmarks

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
uv run python scripts/validate_project.py docs
```
```powershell
uv run python run_pytest.py tests/regression/test_literature_pure_parameter_regression.py tests/regression/test_literature_binary_kij_regression.py tests/regression/test_figiel_2025_born_parameter_parity.py -q
```

## Branch cleanup after merge

After the PR is merged, run or report the equivalent:

```powershell
git -C C:\Users\Tanner\Documents\git\ePC-SAFT fetch origin main --prune
git -C C:\Users\Tanner\Documents\git\ePC-SAFT pull --ff-only origin main
git push origin --delete codex/literature-benchmark-suite
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
