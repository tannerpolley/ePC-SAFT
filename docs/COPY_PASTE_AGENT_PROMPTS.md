# Copy/Paste Prompts for Agents F-M

Use these after the setup and prerequisite PRs have merged and the prompt registry is present on `origin/main`.

Before each Codex task, create the Codex worktree from updated `main` or `origin/main`. The prompt will instruct the agent to switch/create the assigned branch before Goal Prep writes files.

Each prompt points the agent to the repo prompt registry and invokes Goal Prep explicitly.

---
## Agent F: Generic speciation solver using ePC-SAFT activities

Branch:

```text
codex/generic-activity-speciation
```

Copy/paste prompt:

````text
[$goal-prep](C:\Users\Tanner\.codex\skills\goalbuddy\SKILL.md)

Task key: F
GitHub issue: #89

Read:
docs/roadmaps/agent_prompts/index.yaml

Load the prompt file for task key F.
Read GitHub issue #89 and record the current issue scope in docs/goals/<slug>/notes/issue_scope.md during /goal execution.

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
Create/update docs/goals/<slug>/notes/watch_dependency.ps1 during Goal Prep. If dependencies are not satisfied, create the GoalBuddy board, run the bounded watcher, write PREPARED_WAITING only on timeout/blocker, and stop.
If dependencies are satisfied or the watcher returns GATE_PASS, continue implementation without asking.
Do not stop at PREPARED_READY merely because auto_start_after_gate is false; that manual checkpoint mode is retired unless explicitly requested.
If a required value is missing, choose the safest conservative default and record it in docs/goals/<slug>/notes/assumptions.md.
If the missing value prevents safe work, stop with BLOCKED_MISSING_INPUT.

Do not create local .worktrees/.

Do not add application-specific public APIs.
No finite difference.

Task completion automation:
- Open a focused draft PR if one does not already exist for the task branch, and include `Closes #89` in the PR body.
- Review the PR against `origin/main` for task scope, no finite differences, no application-specific public APIs, and no silent scope narrowing.
- If validation passes, GitHub checks pass, the PR is mergeable, and the final GoalBuddy audit says `full_outcome_complete: true`, mark the PR ready and merge without another yes.
- After merge, confirm the corresponding GitHub issue is closed; if it is still open, close it with a comment that names the merged PR and merge commit.
- Delete the remote branch and switch off/delete the local task branch; record `BLOCKED_ISSUE_CLOSE_FAILED`, `BLOCKED_REMOTE_BRANCH_DELETE`, or `BLOCKED_LOCAL_BRANCH_DELETE` if cleanup cannot be completed safely.
After Goal Prep, follow the behavior defined in the loaded task prompt.
````
---

## Agent G: Generic VLE/fugacity-equilibrium solver for volatile neutral species

Branch:

```text
codex/generic-vle-fugacity-equilibrium
```

Copy/paste prompt:

````text
[$goal-prep](C:\Users\Tanner\.codex\skills\goalbuddy\SKILL.md)

Task key: G
GitHub issue: #90

Read:
docs/roadmaps/agent_prompts/index.yaml

Load the prompt file for task key G.
Read GitHub issue #90 and record the current issue scope in docs/goals/<slug>/notes/issue_scope.md during /goal execution.

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
Create/update docs/goals/<slug>/notes/watch_dependency.ps1 during Goal Prep. If dependencies are not satisfied, create the GoalBuddy board, run the bounded watcher, write PREPARED_WAITING only on timeout/blocker, and stop.
If dependencies are satisfied or the watcher returns GATE_PASS, continue implementation without asking.
Do not stop at PREPARED_READY merely because auto_start_after_gate is false; that manual checkpoint mode is retired unless explicitly requested.
If a required value is missing, choose the safest conservative default and record it in docs/goals/<slug>/notes/assumptions.md.
If the missing value prevents safe work, stop with BLOCKED_MISSING_INPUT.

Do not create local .worktrees/.

Do not add application-specific public APIs.
No finite difference.

Task completion automation:
- Open a focused draft PR if one does not already exist for the task branch, and include `Closes #90` in the PR body.
- Review the PR against `origin/main` for task scope, no finite differences, no application-specific public APIs, and no silent scope narrowing.
- If validation passes, GitHub checks pass, the PR is mergeable, and the final GoalBuddy audit says `full_outcome_complete: true`, mark the PR ready and merge without another yes.
- After merge, confirm the corresponding GitHub issue is closed; if it is still open, close it with a comment that names the merged PR and merge commit.
- Delete the remote branch and switch off/delete the local task branch; record `BLOCKED_ISSUE_CLOSE_FAILED`, `BLOCKED_REMOTE_BRANCH_DELETE`, or `BLOCKED_LOCAL_BRANCH_DELETE` if cleanup cannot be completed safely.
After Goal Prep, follow the behavior defined in the loaded task prompt.
````
---

## Agent H: Generic non-electrolyte LLE benchmark and solver hardening

Branch:

```text
codex/generic-non-electrolyte-lle
```

Copy/paste prompt:

````text
[$goal-prep](C:\Users\Tanner\.codex\skills\goalbuddy\SKILL.md)

Task key: H
GitHub issue: #91

Read:
docs/roadmaps/agent_prompts/index.yaml

Load the prompt file for task key H.
Read GitHub issue #91 and record the current issue scope in docs/goals/<slug>/notes/issue_scope.md during /goal execution.

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
Create/update docs/goals/<slug>/notes/watch_dependency.ps1 during Goal Prep. If dependencies are not satisfied, create the GoalBuddy board, run the bounded watcher, write PREPARED_WAITING only on timeout/blocker, and stop.
If dependencies are satisfied or the watcher returns GATE_PASS, continue implementation without asking.
Do not stop at PREPARED_READY merely because auto_start_after_gate is false; that manual checkpoint mode is retired unless explicitly requested.
If a required value is missing, choose the safest conservative default and record it in docs/goals/<slug>/notes/assumptions.md.
If the missing value prevents safe work, stop with BLOCKED_MISSING_INPUT.

Do not create local .worktrees/.

Do not add application-specific public APIs.
No finite difference.

Task completion automation:
- Open a focused draft PR if one does not already exist for the task branch, and include `Closes #91` in the PR body.
- Review the PR against `origin/main` for task scope, no finite differences, no application-specific public APIs, and no silent scope narrowing.
- If validation passes, GitHub checks pass, the PR is mergeable, and the final GoalBuddy audit says `full_outcome_complete: true`, mark the PR ready and merge without another yes.
- After merge, confirm the corresponding GitHub issue is closed; if it is still open, close it with a comment that names the merged PR and merge commit.
- Delete the remote branch and switch off/delete the local task branch; record `BLOCKED_ISSUE_CLOSE_FAILED`, `BLOCKED_REMOTE_BRANCH_DELETE`, or `BLOCKED_LOCAL_BRANCH_DELETE` if cleanup cannot be completed safely.
After Goal Prep, follow the behavior defined in the loaded task prompt.
````
---

## Agent I: Generic electrolyte LLE with distributed ions

Branch:

```text
codex/generic-electrolyte-lle
```

Copy/paste prompt:

````text
[$goal-prep](C:\Users\Tanner\.codex\skills\goalbuddy\SKILL.md)

Task key: I
GitHub issue: #92

Read:
docs/roadmaps/agent_prompts/index.yaml

Load the prompt file for task key I.
Read GitHub issue #92 and record the current issue scope in docs/goals/<slug>/notes/issue_scope.md during /goal execution.

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
Create/update docs/goals/<slug>/notes/watch_dependency.ps1 during Goal Prep. If dependencies are not satisfied, create the GoalBuddy board, run the bounded watcher, write PREPARED_WAITING only on timeout/blocker, and stop.
If dependencies are satisfied or the watcher returns GATE_PASS, continue implementation without asking.
Do not stop at PREPARED_READY merely because auto_start_after_gate is false; that manual checkpoint mode is retired unless explicitly requested.
If a required value is missing, choose the safest conservative default and record it in docs/goals/<slug>/notes/assumptions.md.
If the missing value prevents safe work, stop with BLOCKED_MISSING_INPUT.

Do not create local .worktrees/.

Do not add application-specific public APIs.
No finite difference.

Task completion automation:
- Open a focused draft PR if one does not already exist for the task branch, and include `Closes #92` in the PR body.
- Review the PR against `origin/main` for task scope, no finite differences, no application-specific public APIs, and no silent scope narrowing.
- If validation passes, GitHub checks pass, the PR is mergeable, and the final GoalBuddy audit says `full_outcome_complete: true`, mark the PR ready and merge without another yes.
- After merge, confirm the corresponding GitHub issue is closed; if it is still open, close it with a comment that names the merged PR and merge commit.
- Delete the remote branch and switch off/delete the local task branch; record `BLOCKED_ISSUE_CLOSE_FAILED`, `BLOCKED_REMOTE_BRANCH_DELETE`, or `BLOCKED_LOCAL_BRANCH_DELETE` if cleanup cannot be completed safely.
After Goal Prep, follow the behavior defined in the loaded task prompt.
````
---

## Agent J: Generic reactive LLE and chemical phase equilibrium

Branch:

```text
codex/generic-reactive-lle
```

Copy/paste prompt:

````text
[$goal-prep](C:\Users\Tanner\.codex\skills\goalbuddy\SKILL.md)

Task key: J
GitHub issue: #93

Read:
docs/roadmaps/agent_prompts/index.yaml

Load the prompt file for task key J.
Read GitHub issue #93 and record the current issue scope in docs/goals/<slug>/notes/issue_scope.md during /goal execution.

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
Create/update docs/goals/<slug>/notes/watch_dependency.ps1 during Goal Prep. If dependencies are not satisfied, create the GoalBuddy board, run the bounded watcher, write PREPARED_WAITING only on timeout/blocker, and stop.
If dependencies are satisfied or the watcher returns GATE_PASS, continue implementation without asking.
Do not stop at PREPARED_READY merely because auto_start_after_gate is false; that manual checkpoint mode is retired unless explicitly requested.
If a required value is missing, choose the safest conservative default and record it in docs/goals/<slug>/notes/assumptions.md.
If the missing value prevents safe work, stop with BLOCKED_MISSING_INPUT.

Do not create local .worktrees/.

Do not add application-specific public APIs.
No finite difference.

Task completion automation:
- Open a focused draft PR if one does not already exist for the task branch, and include `Closes #93` in the PR body.
- Review the PR against `origin/main` for task scope, no finite differences, no application-specific public APIs, and no silent scope narrowing.
- If validation passes, GitHub checks pass, the PR is mergeable, and the final GoalBuddy audit says `full_outcome_complete: true`, mark the PR ready and merge without another yes.
- After merge, confirm the corresponding GitHub issue is closed; if it is still open, close it with a comment that names the merged PR and merge commit.
- Delete the remote branch and switch off/delete the local task branch; record `BLOCKED_ISSUE_CLOSE_FAILED`, `BLOCKED_REMOTE_BRANCH_DELETE`, or `BLOCKED_LOCAL_BRANCH_DELETE` if cleanup cannot be completed safely.
After Goal Prep, follow the behavior defined in the loaded task prompt.
````
---

## Agent K: Generic regression row schema and native optimizer backend

Branch:

```text
codex/generic-regression-backend
```

Copy/paste prompt:

````text
[$goal-prep](C:\Users\Tanner\.codex\skills\goalbuddy\SKILL.md)

Task key: K
GitHub issue: #94

Read:
docs/roadmaps/agent_prompts/index.yaml

Load the prompt file for task key K.
Read GitHub issue #94 and record the current issue scope in docs/goals/<slug>/notes/issue_scope.md during /goal execution.

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
Create/update docs/goals/<slug>/notes/watch_dependency.ps1 during Goal Prep. If dependencies are not satisfied, create the GoalBuddy board, run the bounded watcher, write PREPARED_WAITING only on timeout/blocker, and stop.
If dependencies are satisfied or the watcher returns GATE_PASS, continue implementation without asking.
Do not stop at PREPARED_READY merely because auto_start_after_gate is false; that manual checkpoint mode is retired unless explicitly requested.
If a required value is missing, choose the safest conservative default and record it in docs/goals/<slug>/notes/assumptions.md.
If the missing value prevents safe work, stop with BLOCKED_MISSING_INPUT.

Do not create local .worktrees/.

Do not add application-specific public APIs.
No finite difference.

Task completion automation:
- Open a focused draft PR if one does not already exist for the task branch, and include `Closes #94` in the PR body.
- Review the PR against `origin/main` for task scope, no finite differences, no application-specific public APIs, and no silent scope narrowing.
- If validation passes, GitHub checks pass, the PR is mergeable, and the final GoalBuddy audit says `full_outcome_complete: true`, mark the PR ready and merge without another yes.
- After merge, confirm the corresponding GitHub issue is closed; if it is still open, close it with a comment that names the merged PR and merge commit.
- Delete the remote branch and switch off/delete the local task branch; record `BLOCKED_ISSUE_CLOSE_FAILED`, `BLOCKED_REMOTE_BRANCH_DELETE`, or `BLOCKED_LOCAL_BRANCH_DELETE` if cleanup cannot be completed safely.
After Goal Prep, follow the behavior defined in the loaded task prompt.
````
---

## Agent L: Literature benchmark suite

Branch:

```text
codex/literature-benchmark-suite
```

Copy/paste prompt:

````text
[$goal-prep](C:\Users\Tanner\.codex\skills\goalbuddy\SKILL.md)

Task key: L
GitHub issue: #95

Read:
docs/roadmaps/agent_prompts/index.yaml

Load the prompt file for task key L.
Read GitHub issue #95 and record the current issue scope in docs/goals/<slug>/notes/issue_scope.md during /goal execution.

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
Create/update docs/goals/<slug>/notes/watch_dependency.ps1 during Goal Prep. If dependencies are not satisfied, create the GoalBuddy board, run the bounded watcher, write PREPARED_WAITING only on timeout/blocker, and stop.
If dependencies are satisfied or the watcher returns GATE_PASS, continue implementation without asking.
Do not stop at PREPARED_READY merely because auto_start_after_gate is false; that manual checkpoint mode is retired unless explicitly requested.
If a required value is missing, choose the safest conservative default and record it in docs/goals/<slug>/notes/assumptions.md.
If the missing value prevents safe work, stop with BLOCKED_MISSING_INPUT.

Do not create local .worktrees/.

Do not add application-specific public APIs.
No finite difference.

Task completion automation:
- Open a focused draft PR if one does not already exist for the task branch, and include `Closes #95` in the PR body.
- Review the PR against `origin/main` for task scope, no finite differences, no application-specific public APIs, and no silent scope narrowing.
- If validation passes, GitHub checks pass, the PR is mergeable, and the final GoalBuddy audit says `full_outcome_complete: true`, mark the PR ready and merge without another yes.
- After merge, confirm the corresponding GitHub issue is closed; if it is still open, close it with a comment that names the merged PR and merge commit.
- Delete the remote branch and switch off/delete the local task branch; record `BLOCKED_ISSUE_CLOSE_FAILED`, `BLOCKED_REMOTE_BRANCH_DELETE`, or `BLOCKED_LOCAL_BRANCH_DELETE` if cleanup cannot be completed safely.
After Goal Prep, follow the behavior defined in the loaded task prompt.
````
---

## Agent M: Downstream integration smoke tests

Branch:

```text
codex/downstream-integration-smokes
```

Copy/paste prompt:

````text
[$goal-prep](C:\Users\Tanner\.codex\skills\goalbuddy\SKILL.md)

Task key: M
GitHub issue: #96

Read:
docs/roadmaps/agent_prompts/index.yaml

Load the prompt file for task key M.
Read GitHub issue #96 and record the current issue scope in docs/goals/<slug>/notes/issue_scope.md during /goal execution.

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
Create/update docs/goals/<slug>/notes/watch_dependency.ps1 during Goal Prep. If dependencies are not satisfied, create the GoalBuddy board, run the bounded watcher, write PREPARED_WAITING only on timeout/blocker, and stop.
If dependencies are satisfied or the watcher returns GATE_PASS, continue implementation without asking.
Do not stop at PREPARED_READY merely because auto_start_after_gate is false; that manual checkpoint mode is retired unless explicitly requested.
If a required value is missing, choose the safest conservative default and record it in docs/goals/<slug>/notes/assumptions.md.
If the missing value prevents safe work, stop with BLOCKED_MISSING_INPUT.

Do not create local .worktrees/.

Do not add application-specific public APIs.
No finite difference.

Task completion automation:
- Open a focused draft PR if one does not already exist for the task branch, and include `Closes #96` in the PR body.
- Review the PR against `origin/main` for task scope, no finite differences, no application-specific public APIs, and no silent scope narrowing.
- If validation passes, GitHub checks pass, the PR is mergeable, and the final GoalBuddy audit says `full_outcome_complete: true`, mark the PR ready and merge without another yes.
- After merge, confirm the corresponding GitHub issue is closed; if it is still open, close it with a comment that names the merged PR and merge commit.
- Delete the remote branch and switch off/delete the local task branch; record `BLOCKED_ISSUE_CLOSE_FAILED`, `BLOCKED_REMOTE_BRANCH_DELETE`, or `BLOCKED_LOCAL_BRANCH_DELETE` if cleanup cannot be completed safely.
After Goal Prep, follow the behavior defined in the loaded task prompt.
````
