# Copy/Paste Prompts for Agents A-M

Use these after the setup PR has merged and the prompt registry is present on `origin/main`.

Before each Codex task, create the Codex worktree from updated `main` or `origin/main`. The prompt will instruct the agent to switch/create the assigned branch before Goal Prep writes files.

Each prompt points the agent to the repo prompt registry and invokes Goal Prep explicitly.

---

## Agent A: Derivative backend completion audit and coverage matrix hard gate

Branch:

```text
codex/backend-coverage-hard-gate
```

Copy/paste prompt:

```text
[$goal-prep](C:\Users\Tanner\.codex\skills\goalbuddy\SKILL.md)

Task key: A

Read:
docs/roadmaps/agent_prompts/index.yaml

Load the prompt file for task key A.

Run Goal Prep for this exact task:
Derivative backend completion audit and coverage matrix hard gate

Use the local live GoalBuddy board in Codex.

Do not start implementation until the dependency gate in the task prompt passes.

The Codex worktree may have started on `main`. Before Goal Prep writes files, switch/create the assigned branch.

Assigned branch:
codex/backend-coverage-hard-gate

Branch bootstrap:

```powershell
git fetch origin --prune
$current = (git branch --show-current).Trim()
if ($current -ne "codex/backend-coverage-hard-gate") {
    git ls-remote --exit-code --heads origin codex/backend-coverage-hard-gate
    if ($LASTEXITCODE -eq 0) {
        git switch codex/backend-coverage-hard-gate 2>$null
        if ($LASTEXITCODE -ne 0) {
            git switch -c codex/backend-coverage-hard-gate --track origin/codex/backend-coverage-hard-gate
        }
    } else {
        git switch -c codex/backend-coverage-hard-gate origin/main
    }
}
$current = (git branch --show-current).Trim()
if ($current -ne "codex/backend-coverage-hard-gate") {
    Write-Error "BRANCH_BOOTSTRAP_FAILED expected=codex/backend-coverage-hard-gate actual=$current"
    exit 1
}
git status --short
```

After branch bootstrap, current branch must be:
codex/backend-coverage-hard-gate

Dependencies:
none

Non-interaction rule:
Do not ask me to confirm one of several options.
Do not ask repeated confirmation questions.
Use the defaults written in the task prompt.
Use bounded watcher mode for dependency waits.
Set `watcher_mode: bounded`, `auto_start_after_gate: true`, `poll_interval_seconds: 120`, and `max_wait_minutes: 480` in the local GoalBuddy/dependency files.
If dependencies are not satisfied, create the GoalBuddy board, create/update the per-goal watcher, run the bounded watcher, write PREPARED_WAITING only on timeout/blocker, and stop.
If dependencies are satisfied or the watcher returns GATE_PASS, continue implementation without asking.
Do not stop at PREPARED_READY merely because auto_start_after_gate is false; that manual checkpoint mode is retired unless explicitly requested.
If a required value is missing, choose the safest conservative default and record it in docs/goals/<slug>/notes/assumptions.md.
If the missing value prevents safe work, stop with BLOCKED_MISSING_INPUT.

Do not create local .worktrees/.
Do not use PR #56 as a base.
Do not add application-specific public APIs.
No finite difference.

After Goal Prep, follow the behavior defined in the loaded task prompt.
```
---

## Agent B: Explicit CppAD parameter derivatives for EOS/property APIs

Branch:

```text
codex/cppad-explicit-parameter-derivatives
```

Copy/paste prompt:

```text
[$goal-prep](C:\Users\Tanner\.codex\skills\goalbuddy\SKILL.md)

Task key: B

Read:
docs/roadmaps/agent_prompts/index.yaml

Load the prompt file for task key B.

Run Goal Prep for this exact task:
Explicit CppAD parameter derivatives for EOS/property APIs

Use the local live GoalBuddy board in Codex.

Do not start implementation until the dependency gate in the task prompt passes.

The Codex worktree may have started on `main`. Before Goal Prep writes files, switch/create the assigned branch.

Assigned branch:
codex/cppad-explicit-parameter-derivatives

Branch bootstrap:

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

After branch bootstrap, current branch must be:
codex/cppad-explicit-parameter-derivatives

Dependencies:
A

Non-interaction rule:
Do not ask me to confirm one of several options.
Do not ask repeated confirmation questions.
Use the defaults written in the task prompt.
Use bounded watcher mode for dependency waits.
Set `watcher_mode: bounded`, `auto_start_after_gate: true`, `poll_interval_seconds: 120`, and `max_wait_minutes: 480` in the local GoalBuddy/dependency files.
If dependencies are not satisfied, create the GoalBuddy board, create/update the per-goal watcher, run the bounded watcher, write PREPARED_WAITING only on timeout/blocker, and stop.
If dependencies are satisfied or the watcher returns GATE_PASS, continue implementation without asking.
Do not stop at PREPARED_READY merely because auto_start_after_gate is false; that manual checkpoint mode is retired unless explicitly requested.
If a required value is missing, choose the safest conservative default and record it in docs/goals/<slug>/notes/assumptions.md.
If the missing value prevents safe work, stop with BLOCKED_MISSING_INPUT.

Do not create local .worktrees/.
Do not use PR #56 as a base.
Do not add application-specific public APIs.
No finite difference.

After Goal Prep, follow the behavior defined in the loaded task prompt.
```
---

## Agent C: Generic implicit sensitivity framework for solved states

Branch:

```text
codex/generic-implicit-sensitivity-framework
```

Copy/paste prompt:

```text
[$goal-prep](C:\Users\Tanner\.codex\skills\goalbuddy\SKILL.md)

Task key: C

Read:
docs/roadmaps/agent_prompts/index.yaml

Load the prompt file for task key C.

Run Goal Prep for this exact task:
Generic implicit sensitivity framework for solved states

Use the local live GoalBuddy board in Codex.

Do not start implementation until the dependency gate in the task prompt passes.

The Codex worktree may have started on `main`. Before Goal Prep writes files, switch/create the assigned branch.

Assigned branch:
codex/generic-implicit-sensitivity-framework

Branch bootstrap:

```powershell
git fetch origin --prune
$current = (git branch --show-current).Trim()
if ($current -ne "codex/generic-implicit-sensitivity-framework") {
    git ls-remote --exit-code --heads origin codex/generic-implicit-sensitivity-framework
    if ($LASTEXITCODE -eq 0) {
        git switch codex/generic-implicit-sensitivity-framework 2>$null
        if ($LASTEXITCODE -ne 0) {
            git switch -c codex/generic-implicit-sensitivity-framework --track origin/codex/generic-implicit-sensitivity-framework
        }
    } else {
        git switch -c codex/generic-implicit-sensitivity-framework origin/main
    }
}
$current = (git branch --show-current).Trim()
if ($current -ne "codex/generic-implicit-sensitivity-framework") {
    Write-Error "BRANCH_BOOTSTRAP_FAILED expected=codex/generic-implicit-sensitivity-framework actual=$current"
    exit 1
}
git status --short
```

After branch bootstrap, current branch must be:
codex/generic-implicit-sensitivity-framework

Dependencies:
A, B

Non-interaction rule:
Do not ask me to confirm one of several options.
Do not ask repeated confirmation questions.
Use the defaults written in the task prompt.
Use bounded watcher mode for dependency waits.
Set `watcher_mode: bounded`, `auto_start_after_gate: true`, `poll_interval_seconds: 120`, and `max_wait_minutes: 480` in the local GoalBuddy/dependency files.
If dependencies are not satisfied, create the GoalBuddy board, create/update the per-goal watcher, run the bounded watcher, write PREPARED_WAITING only on timeout/blocker, and stop.
If dependencies are satisfied or the watcher returns GATE_PASS, continue implementation without asking.
Do not stop at PREPARED_READY merely because auto_start_after_gate is false; that manual checkpoint mode is retired unless explicitly requested.
If a required value is missing, choose the safest conservative default and record it in docs/goals/<slug>/notes/assumptions.md.
If the missing value prevents safe work, stop with BLOCKED_MISSING_INPUT.

Do not create local .worktrees/.
Do not use PR #56 as a base.
Do not add application-specific public APIs.
No finite difference.

After Goal Prep, follow the behavior defined in the loaded task prompt.
```
---

## Agent D: General reaction and equilibrium-constant convention layer

Branch:

```text
codex/reaction-constant-conventions
```

Copy/paste prompt:

```text
[$goal-prep](C:\Users\Tanner\.codex\skills\goalbuddy\SKILL.md)

Task key: D

Read:
docs/roadmaps/agent_prompts/index.yaml

Load the prompt file for task key D.

Run Goal Prep for this exact task:
General reaction and equilibrium-constant convention layer

Use the local live GoalBuddy board in Codex.

Do not start implementation until the dependency gate in the task prompt passes.

The Codex worktree may have started on `main`. Before Goal Prep writes files, switch/create the assigned branch.

Assigned branch:
codex/reaction-constant-conventions

Branch bootstrap:

```powershell
git fetch origin --prune
$current = (git branch --show-current).Trim()
if ($current -ne "codex/reaction-constant-conventions") {
    git ls-remote --exit-code --heads origin codex/reaction-constant-conventions
    if ($LASTEXITCODE -eq 0) {
        git switch codex/reaction-constant-conventions 2>$null
        if ($LASTEXITCODE -ne 0) {
            git switch -c codex/reaction-constant-conventions --track origin/codex/reaction-constant-conventions
        }
    } else {
        git switch -c codex/reaction-constant-conventions origin/main
    }
}
$current = (git branch --show-current).Trim()
if ($current -ne "codex/reaction-constant-conventions") {
    Write-Error "BRANCH_BOOTSTRAP_FAILED expected=codex/reaction-constant-conventions actual=$current"
    exit 1
}
git status --short
```

After branch bootstrap, current branch must be:
codex/reaction-constant-conventions

Dependencies:
A

Non-interaction rule:
Do not ask me to confirm one of several options.
Do not ask repeated confirmation questions.
Use the defaults written in the task prompt.
Use bounded watcher mode for dependency waits.
Set `watcher_mode: bounded`, `auto_start_after_gate: true`, `poll_interval_seconds: 120`, and `max_wait_minutes: 480` in the local GoalBuddy/dependency files.
If dependencies are not satisfied, create the GoalBuddy board, create/update the per-goal watcher, run the bounded watcher, write PREPARED_WAITING only on timeout/blocker, and stop.
If dependencies are satisfied or the watcher returns GATE_PASS, continue implementation without asking.
Do not stop at PREPARED_READY merely because auto_start_after_gate is false; that manual checkpoint mode is retired unless explicitly requested.
If a required value is missing, choose the safest conservative default and record it in docs/goals/<slug>/notes/assumptions.md.
If the missing value prevents safe work, stop with BLOCKED_MISSING_INPUT.

Do not create local .worktrees/.
Do not use PR #56 as a base.
Do not add application-specific public APIs.
No finite difference.

After Goal Prep, follow the behavior defined in the loaded task prompt.
```
---

## Agent E: Generic target-row and dataset schema

Branch:

```text
codex/generic-target-row-schema
```

Copy/paste prompt:

```text
[$goal-prep](C:\Users\Tanner\.codex\skills\goalbuddy\SKILL.md)

Task key: E

Read:
docs/roadmaps/agent_prompts/index.yaml

Load the prompt file for task key E.

Run Goal Prep for this exact task:
Generic target-row and dataset schema

Use the local live GoalBuddy board in Codex.

Do not start implementation until the dependency gate in the task prompt passes.

The Codex worktree may have started on `main`. Before Goal Prep writes files, switch/create the assigned branch.

Assigned branch:
codex/generic-target-row-schema

Branch bootstrap:

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

After branch bootstrap, current branch must be:
codex/generic-target-row-schema

Dependencies:
A

Non-interaction rule:
Do not ask me to confirm one of several options.
Do not ask repeated confirmation questions.
Use the defaults written in the task prompt.
Use bounded watcher mode for dependency waits.
Set `watcher_mode: bounded`, `auto_start_after_gate: true`, `poll_interval_seconds: 120`, and `max_wait_minutes: 480` in the local GoalBuddy/dependency files.
If dependencies are not satisfied, create the GoalBuddy board, create/update the per-goal watcher, run the bounded watcher, write PREPARED_WAITING only on timeout/blocker, and stop.
If dependencies are satisfied or the watcher returns GATE_PASS, continue implementation without asking.
Do not stop at PREPARED_READY merely because auto_start_after_gate is false; that manual checkpoint mode is retired unless explicitly requested.
If a required value is missing, choose the safest conservative default and record it in docs/goals/<slug>/notes/assumptions.md.
If the missing value prevents safe work, stop with BLOCKED_MISSING_INPUT.

Do not create local .worktrees/.
Do not use PR #56 as a base.
Do not add application-specific public APIs.
No finite difference.

After Goal Prep, follow the behavior defined in the loaded task prompt.
```
---

## Agent F: Generic speciation solver using ePC-SAFT activities

Branch:

```text
codex/generic-activity-speciation
```

Copy/paste prompt:

```text
[$goal-prep](C:\Users\Tanner\.codex\skills\goalbuddy\SKILL.md)

Task key: F

Read:
docs/roadmaps/agent_prompts/index.yaml

Load the prompt file for task key F.

Run Goal Prep for this exact task:
Generic speciation solver using ePC-SAFT activities

Use the local live GoalBuddy board in Codex.

Do not start implementation until the dependency gate in the task prompt passes.

The Codex worktree may have started on `main`. Before Goal Prep writes files, switch/create the assigned branch.

Assigned branch:
codex/generic-activity-speciation

Branch bootstrap:

```powershell
git fetch origin --prune
$current = (git branch --show-current).Trim()
if ($current -ne "codex/generic-activity-speciation") {
    git ls-remote --exit-code --heads origin codex/generic-activity-speciation
    if ($LASTEXITCODE -eq 0) {
        git switch codex/generic-activity-speciation 2>$null
        if ($LASTEXITCODE -ne 0) {
            git switch -c codex/generic-activity-speciation --track origin/codex/generic-activity-speciation
        }
    } else {
        git switch -c codex/generic-activity-speciation origin/main
    }
}
$current = (git branch --show-current).Trim()
if ($current -ne "codex/generic-activity-speciation") {
    Write-Error "BRANCH_BOOTSTRAP_FAILED expected=codex/generic-activity-speciation actual=$current"
    exit 1
}
git status --short
```

After branch bootstrap, current branch must be:
codex/generic-activity-speciation

Dependencies:
C, D

Non-interaction rule:
Do not ask me to confirm one of several options.
Do not ask repeated confirmation questions.
Use the defaults written in the task prompt.
Use bounded watcher mode for dependency waits.
Set `watcher_mode: bounded`, `auto_start_after_gate: true`, `poll_interval_seconds: 120`, and `max_wait_minutes: 480` in the local GoalBuddy/dependency files.
If dependencies are not satisfied, create the GoalBuddy board, create/update the per-goal watcher, run the bounded watcher, write PREPARED_WAITING only on timeout/blocker, and stop.
If dependencies are satisfied or the watcher returns GATE_PASS, continue implementation without asking.
Do not stop at PREPARED_READY merely because auto_start_after_gate is false; that manual checkpoint mode is retired unless explicitly requested.
If a required value is missing, choose the safest conservative default and record it in docs/goals/<slug>/notes/assumptions.md.
If the missing value prevents safe work, stop with BLOCKED_MISSING_INPUT.

Do not create local .worktrees/.
Do not use PR #56 as a base.
Do not add application-specific public APIs.
No finite difference.

Task completion automation:
- Open a focused draft PR if one does not already exist for the task branch.
- Review the PR against `origin/main` for task scope, no finite differences, no application-specific public APIs, no PR #56 dependency, and no silent scope narrowing.
- If validation passes, GitHub checks pass, the PR is mergeable, and the final GoalBuddy audit says `full_outcome_complete: true`, mark the PR ready and merge without another yes.
- After merge, delete the remote branch and switch off/delete the local task branch; record `BLOCKED_REMOTE_BRANCH_DELETE` or `BLOCKED_LOCAL_BRANCH_DELETE` if cleanup cannot be completed safely.
After Goal Prep, follow the behavior defined in the loaded task prompt.
```
---

## Agent G: Generic VLE/fugacity-equilibrium solver for volatile neutral species

Branch:

```text
codex/generic-vle-fugacity-equilibrium
```

Copy/paste prompt:

```text
[$goal-prep](C:\Users\Tanner\.codex\skills\goalbuddy\SKILL.md)

Task key: G

Read:
docs/roadmaps/agent_prompts/index.yaml

Load the prompt file for task key G.

Run Goal Prep for this exact task:
Generic VLE/fugacity-equilibrium solver for volatile neutral species

Use the local live GoalBuddy board in Codex.

Do not start implementation until the dependency gate in the task prompt passes.

The Codex worktree may have started on `main`. Before Goal Prep writes files, switch/create the assigned branch.

Assigned branch:
codex/generic-vle-fugacity-equilibrium

Branch bootstrap:

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

After branch bootstrap, current branch must be:
codex/generic-vle-fugacity-equilibrium

Dependencies:
C

Non-interaction rule:
Do not ask me to confirm one of several options.
Do not ask repeated confirmation questions.
Use the defaults written in the task prompt.
Use bounded watcher mode for dependency waits.
Set `watcher_mode: bounded`, `auto_start_after_gate: true`, `poll_interval_seconds: 120`, and `max_wait_minutes: 480` in the local GoalBuddy/dependency files.
If dependencies are not satisfied, create the GoalBuddy board, create/update the per-goal watcher, run the bounded watcher, write PREPARED_WAITING only on timeout/blocker, and stop.
If dependencies are satisfied or the watcher returns GATE_PASS, continue implementation without asking.
Do not stop at PREPARED_READY merely because auto_start_after_gate is false; that manual checkpoint mode is retired unless explicitly requested.
If a required value is missing, choose the safest conservative default and record it in docs/goals/<slug>/notes/assumptions.md.
If the missing value prevents safe work, stop with BLOCKED_MISSING_INPUT.

Do not create local .worktrees/.
Do not use PR #56 as a base.
Do not add application-specific public APIs.
No finite difference.

Task completion automation:
- Open a focused draft PR if one does not already exist for the task branch.
- Review the PR against `origin/main` for task scope, no finite differences, no application-specific public APIs, no PR #56 dependency, and no silent scope narrowing.
- If validation passes, GitHub checks pass, the PR is mergeable, and the final GoalBuddy audit says `full_outcome_complete: true`, mark the PR ready and merge without another yes.
- After merge, delete the remote branch and switch off/delete the local task branch; record `BLOCKED_REMOTE_BRANCH_DELETE` or `BLOCKED_LOCAL_BRANCH_DELETE` if cleanup cannot be completed safely.
After Goal Prep, follow the behavior defined in the loaded task prompt.
```
---

## Agent H: Generic non-electrolyte LLE benchmark and solver hardening

Branch:

```text
codex/generic-non-electrolyte-lle
```

Copy/paste prompt:

```text
[$goal-prep](C:\Users\Tanner\.codex\skills\goalbuddy\SKILL.md)

Task key: H

Read:
docs/roadmaps/agent_prompts/index.yaml

Load the prompt file for task key H.

Run Goal Prep for this exact task:
Generic non-electrolyte LLE benchmark and solver hardening

Use the local live GoalBuddy board in Codex.

Do not start implementation until the dependency gate in the task prompt passes.

The Codex worktree may have started on `main`. Before Goal Prep writes files, switch/create the assigned branch.

Assigned branch:
codex/generic-non-electrolyte-lle

Branch bootstrap:

```powershell
git fetch origin --prune
$current = (git branch --show-current).Trim()
if ($current -ne "codex/generic-non-electrolyte-lle") {
    git ls-remote --exit-code --heads origin codex/generic-non-electrolyte-lle
    if ($LASTEXITCODE -eq 0) {
        git switch codex/generic-non-electrolyte-lle 2>$null
        if ($LASTEXITCODE -ne 0) {
            git switch -c codex/generic-non-electrolyte-lle --track origin/codex/generic-non-electrolyte-lle
        }
    } else {
        git switch -c codex/generic-non-electrolyte-lle origin/main
    }
}
$current = (git branch --show-current).Trim()
if ($current -ne "codex/generic-non-electrolyte-lle") {
    Write-Error "BRANCH_BOOTSTRAP_FAILED expected=codex/generic-non-electrolyte-lle actual=$current"
    exit 1
}
git status --short
```

After branch bootstrap, current branch must be:
codex/generic-non-electrolyte-lle

Dependencies:
C

Non-interaction rule:
Do not ask me to confirm one of several options.
Do not ask repeated confirmation questions.
Use the defaults written in the task prompt.
Use bounded watcher mode for dependency waits.
Set `watcher_mode: bounded`, `auto_start_after_gate: true`, `poll_interval_seconds: 120`, and `max_wait_minutes: 480` in the local GoalBuddy/dependency files.
If dependencies are not satisfied, create the GoalBuddy board, create/update the per-goal watcher, run the bounded watcher, write PREPARED_WAITING only on timeout/blocker, and stop.
If dependencies are satisfied or the watcher returns GATE_PASS, continue implementation without asking.
Do not stop at PREPARED_READY merely because auto_start_after_gate is false; that manual checkpoint mode is retired unless explicitly requested.
If a required value is missing, choose the safest conservative default and record it in docs/goals/<slug>/notes/assumptions.md.
If the missing value prevents safe work, stop with BLOCKED_MISSING_INPUT.

Do not create local .worktrees/.
Do not use PR #56 as a base.
Do not add application-specific public APIs.
No finite difference.

Task completion automation:
- Open a focused draft PR if one does not already exist for the task branch.
- Review the PR against `origin/main` for task scope, no finite differences, no application-specific public APIs, no PR #56 dependency, and no silent scope narrowing.
- If validation passes, GitHub checks pass, the PR is mergeable, and the final GoalBuddy audit says `full_outcome_complete: true`, mark the PR ready and merge without another yes.
- After merge, delete the remote branch and switch off/delete the local task branch; record `BLOCKED_REMOTE_BRANCH_DELETE` or `BLOCKED_LOCAL_BRANCH_DELETE` if cleanup cannot be completed safely.
After Goal Prep, follow the behavior defined in the loaded task prompt.
```
---

## Agent I: Generic electrolyte LLE with distributed ions

Branch:

```text
codex/generic-electrolyte-lle
```

Copy/paste prompt:

```text
[$goal-prep](C:\Users\Tanner\.codex\skills\goalbuddy\SKILL.md)

Task key: I

Read:
docs/roadmaps/agent_prompts/index.yaml

Load the prompt file for task key I.

Run Goal Prep for this exact task:
Generic electrolyte LLE with distributed ions

Use the local live GoalBuddy board in Codex.

Do not start implementation until the dependency gate in the task prompt passes.

The Codex worktree may have started on `main`. Before Goal Prep writes files, switch/create the assigned branch.

Assigned branch:
codex/generic-electrolyte-lle

Branch bootstrap:

```powershell
git fetch origin --prune
$current = (git branch --show-current).Trim()
if ($current -ne "codex/generic-electrolyte-lle") {
    git ls-remote --exit-code --heads origin codex/generic-electrolyte-lle
    if ($LASTEXITCODE -eq 0) {
        git switch codex/generic-electrolyte-lle 2>$null
        if ($LASTEXITCODE -ne 0) {
            git switch -c codex/generic-electrolyte-lle --track origin/codex/generic-electrolyte-lle
        }
    } else {
        git switch -c codex/generic-electrolyte-lle origin/main
    }
}
$current = (git branch --show-current).Trim()
if ($current -ne "codex/generic-electrolyte-lle") {
    Write-Error "BRANCH_BOOTSTRAP_FAILED expected=codex/generic-electrolyte-lle actual=$current"
    exit 1
}
git status --short
```

After branch bootstrap, current branch must be:
codex/generic-electrolyte-lle

Dependencies:
H, C

Non-interaction rule:
Do not ask me to confirm one of several options.
Do not ask repeated confirmation questions.
Use the defaults written in the task prompt.
Use bounded watcher mode for dependency waits.
Set `watcher_mode: bounded`, `auto_start_after_gate: true`, `poll_interval_seconds: 120`, and `max_wait_minutes: 480` in the local GoalBuddy/dependency files.
If dependencies are not satisfied, create the GoalBuddy board, create/update the per-goal watcher, run the bounded watcher, write PREPARED_WAITING only on timeout/blocker, and stop.
If dependencies are satisfied or the watcher returns GATE_PASS, continue implementation without asking.
Do not stop at PREPARED_READY merely because auto_start_after_gate is false; that manual checkpoint mode is retired unless explicitly requested.
If a required value is missing, choose the safest conservative default and record it in docs/goals/<slug>/notes/assumptions.md.
If the missing value prevents safe work, stop with BLOCKED_MISSING_INPUT.

Do not create local .worktrees/.
Do not use PR #56 as a base.
Do not add application-specific public APIs.
No finite difference.

Task completion automation:
- Open a focused draft PR if one does not already exist for the task branch.
- Review the PR against `origin/main` for task scope, no finite differences, no application-specific public APIs, no PR #56 dependency, and no silent scope narrowing.
- If validation passes, GitHub checks pass, the PR is mergeable, and the final GoalBuddy audit says `full_outcome_complete: true`, mark the PR ready and merge without another yes.
- After merge, delete the remote branch and switch off/delete the local task branch; record `BLOCKED_REMOTE_BRANCH_DELETE` or `BLOCKED_LOCAL_BRANCH_DELETE` if cleanup cannot be completed safely.
After Goal Prep, follow the behavior defined in the loaded task prompt.
```
---

## Agent J: Generic reactive LLE and chemical phase equilibrium

Branch:

```text
codex/generic-reactive-lle
```

Copy/paste prompt:

```text
[$goal-prep](C:\Users\Tanner\.codex\skills\goalbuddy\SKILL.md)

Task key: J

Read:
docs/roadmaps/agent_prompts/index.yaml

Load the prompt file for task key J.

Run Goal Prep for this exact task:
Generic reactive LLE and chemical phase equilibrium

Use the local live GoalBuddy board in Codex.

Do not start implementation until the dependency gate in the task prompt passes.

The Codex worktree may have started on `main`. Before Goal Prep writes files, switch/create the assigned branch.

Assigned branch:
codex/generic-reactive-lle

Branch bootstrap:

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

After branch bootstrap, current branch must be:
codex/generic-reactive-lle

Dependencies:
I, D, F

Non-interaction rule:
Do not ask me to confirm one of several options.
Do not ask repeated confirmation questions.
Use the defaults written in the task prompt.
Use bounded watcher mode for dependency waits.
Set `watcher_mode: bounded`, `auto_start_after_gate: true`, `poll_interval_seconds: 120`, and `max_wait_minutes: 480` in the local GoalBuddy/dependency files.
If dependencies are not satisfied, create the GoalBuddy board, create/update the per-goal watcher, run the bounded watcher, write PREPARED_WAITING only on timeout/blocker, and stop.
If dependencies are satisfied or the watcher returns GATE_PASS, continue implementation without asking.
Do not stop at PREPARED_READY merely because auto_start_after_gate is false; that manual checkpoint mode is retired unless explicitly requested.
If a required value is missing, choose the safest conservative default and record it in docs/goals/<slug>/notes/assumptions.md.
If the missing value prevents safe work, stop with BLOCKED_MISSING_INPUT.

Do not create local .worktrees/.
Do not use PR #56 as a base.
Do not add application-specific public APIs.
No finite difference.

Task completion automation:
- Open a focused draft PR if one does not already exist for the task branch.
- Review the PR against `origin/main` for task scope, no finite differences, no application-specific public APIs, no PR #56 dependency, and no silent scope narrowing.
- If validation passes, GitHub checks pass, the PR is mergeable, and the final GoalBuddy audit says `full_outcome_complete: true`, mark the PR ready and merge without another yes.
- After merge, delete the remote branch and switch off/delete the local task branch; record `BLOCKED_REMOTE_BRANCH_DELETE` or `BLOCKED_LOCAL_BRANCH_DELETE` if cleanup cannot be completed safely.
After Goal Prep, follow the behavior defined in the loaded task prompt.
```
---

## Agent K: Generic regression row schema and native optimizer backend

Branch:

```text
codex/generic-regression-backend
```

Copy/paste prompt:

```text
[$goal-prep](C:\Users\Tanner\.codex\skills\goalbuddy\SKILL.md)

Task key: K

Read:
docs/roadmaps/agent_prompts/index.yaml

Load the prompt file for task key K.

Run Goal Prep for this exact task:
Generic regression row schema and native optimizer backend

Use the local live GoalBuddy board in Codex.

Do not start implementation until the dependency gate in the task prompt passes.

The Codex worktree may have started on `main`. Before Goal Prep writes files, switch/create the assigned branch.

Assigned branch:
codex/generic-regression-backend

Branch bootstrap:

```powershell
git fetch origin --prune
$current = (git branch --show-current).Trim()
if ($current -ne "codex/generic-regression-backend") {
    git ls-remote --exit-code --heads origin codex/generic-regression-backend
    if ($LASTEXITCODE -eq 0) {
        git switch codex/generic-regression-backend 2>$null
        if ($LASTEXITCODE -ne 0) {
            git switch -c codex/generic-regression-backend --track origin/codex/generic-regression-backend
        }
    } else {
        git switch -c codex/generic-regression-backend origin/main
    }
}
$current = (git branch --show-current).Trim()
if ($current -ne "codex/generic-regression-backend") {
    Write-Error "BRANCH_BOOTSTRAP_FAILED expected=codex/generic-regression-backend actual=$current"
    exit 1
}
git status --short
```

After branch bootstrap, current branch must be:
codex/generic-regression-backend

Dependencies:
B, C, E

Non-interaction rule:
Do not ask me to confirm one of several options.
Do not ask repeated confirmation questions.
Use the defaults written in the task prompt.
Use bounded watcher mode for dependency waits.
Set `watcher_mode: bounded`, `auto_start_after_gate: true`, `poll_interval_seconds: 120`, and `max_wait_minutes: 480` in the local GoalBuddy/dependency files.
If dependencies are not satisfied, create the GoalBuddy board, create/update the per-goal watcher, run the bounded watcher, write PREPARED_WAITING only on timeout/blocker, and stop.
If dependencies are satisfied or the watcher returns GATE_PASS, continue implementation without asking.
Do not stop at PREPARED_READY merely because auto_start_after_gate is false; that manual checkpoint mode is retired unless explicitly requested.
If a required value is missing, choose the safest conservative default and record it in docs/goals/<slug>/notes/assumptions.md.
If the missing value prevents safe work, stop with BLOCKED_MISSING_INPUT.

Do not create local .worktrees/.
Do not use PR #56 as a base.
Do not add application-specific public APIs.
No finite difference.

Task completion automation:
- Open a focused draft PR if one does not already exist for the task branch.
- Review the PR against `origin/main` for task scope, no finite differences, no application-specific public APIs, no PR #56 dependency, and no silent scope narrowing.
- If validation passes, GitHub checks pass, the PR is mergeable, and the final GoalBuddy audit says `full_outcome_complete: true`, mark the PR ready and merge without another yes.
- After merge, delete the remote branch and switch off/delete the local task branch; record `BLOCKED_REMOTE_BRANCH_DELETE` or `BLOCKED_LOCAL_BRANCH_DELETE` if cleanup cannot be completed safely.
After Goal Prep, follow the behavior defined in the loaded task prompt.
```
---

## Agent L: Literature benchmark suite

Branch:

```text
codex/literature-benchmark-suite
```

Copy/paste prompt:

```text
[$goal-prep](C:\Users\Tanner\.codex\skills\goalbuddy\SKILL.md)

Task key: L

Read:
docs/roadmaps/agent_prompts/index.yaml

Load the prompt file for task key L.

Run Goal Prep for this exact task:
Literature benchmark suite

Use the local live GoalBuddy board in Codex.

Do not start implementation until the dependency gate in the task prompt passes.

The Codex worktree may have started on `main`. Before Goal Prep writes files, switch/create the assigned branch.

Assigned branch:
codex/literature-benchmark-suite

Branch bootstrap:

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

After branch bootstrap, current branch must be:
codex/literature-benchmark-suite

Dependencies:
none

Non-interaction rule:
Do not ask me to confirm one of several options.
Do not ask repeated confirmation questions.
Use the defaults written in the task prompt.
Use bounded watcher mode for dependency waits.
Set `watcher_mode: bounded`, `auto_start_after_gate: true`, `poll_interval_seconds: 120`, and `max_wait_minutes: 480` in the local GoalBuddy/dependency files.
If dependencies are not satisfied, create the GoalBuddy board, create/update the per-goal watcher, run the bounded watcher, write PREPARED_WAITING only on timeout/blocker, and stop.
If dependencies are satisfied or the watcher returns GATE_PASS, continue implementation without asking.
Do not stop at PREPARED_READY merely because auto_start_after_gate is false; that manual checkpoint mode is retired unless explicitly requested.
If a required value is missing, choose the safest conservative default and record it in docs/goals/<slug>/notes/assumptions.md.
If the missing value prevents safe work, stop with BLOCKED_MISSING_INPUT.

Do not create local .worktrees/.
Do not use PR #56 as a base.
Do not add application-specific public APIs.
No finite difference.

Task completion automation:
- Open a focused draft PR if one does not already exist for the task branch.
- Review the PR against `origin/main` for task scope, no finite differences, no application-specific public APIs, no PR #56 dependency, and no silent scope narrowing.
- If validation passes, GitHub checks pass, the PR is mergeable, and the final GoalBuddy audit says `full_outcome_complete: true`, mark the PR ready and merge without another yes.
- After merge, delete the remote branch and switch off/delete the local task branch; record `BLOCKED_REMOTE_BRANCH_DELETE` or `BLOCKED_LOCAL_BRANCH_DELETE` if cleanup cannot be completed safely.
After Goal Prep, follow the behavior defined in the loaded task prompt.
```
---

## Agent M: Downstream integration smoke tests

Branch:

```text
codex/downstream-integration-smokes
```

Copy/paste prompt:

```text
[$goal-prep](C:\Users\Tanner\.codex\skills\goalbuddy\SKILL.md)

Task key: M

Read:
docs/roadmaps/agent_prompts/index.yaml

Load the prompt file for task key M.

Run Goal Prep for this exact task:
Downstream integration smoke tests

Use the local live GoalBuddy board in Codex.

Do not start implementation until the dependency gate in the task prompt passes.

The Codex worktree may have started on `main`. Before Goal Prep writes files, switch/create the assigned branch.

Assigned branch:
codex/downstream-integration-smokes

Branch bootstrap:

```powershell
git fetch origin --prune
$current = (git branch --show-current).Trim()
if ($current -ne "codex/downstream-integration-smokes") {
    git ls-remote --exit-code --heads origin codex/downstream-integration-smokes
    if ($LASTEXITCODE -eq 0) {
        git switch codex/downstream-integration-smokes 2>$null
        if ($LASTEXITCODE -ne 0) {
            git switch -c codex/downstream-integration-smokes --track origin/codex/downstream-integration-smokes
        }
    } else {
        git switch -c codex/downstream-integration-smokes origin/main
    }
}
$current = (git branch --show-current).Trim()
if ($current -ne "codex/downstream-integration-smokes") {
    Write-Error "BRANCH_BOOTSTRAP_FAILED expected=codex/downstream-integration-smokes actual=$current"
    exit 1
}
git status --short
```

After branch bootstrap, current branch must be:
codex/downstream-integration-smokes

Dependencies:
F, G, I, J, K

Non-interaction rule:
Do not ask me to confirm one of several options.
Do not ask repeated confirmation questions.
Use the defaults written in the task prompt.
Use bounded watcher mode for dependency waits.
Set `watcher_mode: bounded`, `auto_start_after_gate: true`, `poll_interval_seconds: 120`, and `max_wait_minutes: 480` in the local GoalBuddy/dependency files.
If dependencies are not satisfied, create the GoalBuddy board, create/update the per-goal watcher, run the bounded watcher, write PREPARED_WAITING only on timeout/blocker, and stop.
If dependencies are satisfied or the watcher returns GATE_PASS, continue implementation without asking.
Do not stop at PREPARED_READY merely because auto_start_after_gate is false; that manual checkpoint mode is retired unless explicitly requested.
If a required value is missing, choose the safest conservative default and record it in docs/goals/<slug>/notes/assumptions.md.
If the missing value prevents safe work, stop with BLOCKED_MISSING_INPUT.

Do not create local .worktrees/.
Do not use PR #56 as a base.
Do not add application-specific public APIs.
No finite difference.

Task completion automation:
- Open a focused draft PR if one does not already exist for the task branch.
- Review the PR against `origin/main` for task scope, no finite differences, no application-specific public APIs, no PR #56 dependency, and no silent scope narrowing.
- If validation passes, GitHub checks pass, the PR is mergeable, and the final GoalBuddy audit says `full_outcome_complete: true`, mark the PR ready and merge without another yes.
- After merge, delete the remote branch and switch off/delete the local task branch; record `BLOCKED_REMOTE_BRANCH_DELETE` or `BLOCKED_LOCAL_BRANCH_DELETE` if cleanup cannot be completed safely.
After Goal Prep, follow the behavior defined in the loaded task prompt.
```
