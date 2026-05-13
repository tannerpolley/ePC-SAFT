# Setup Codex Agent Prompt

Use this prompt in a new Codex app worktree checked out from updated `main`.

## Before creating the Codex task

Run:

```powershell
git -C C:\Users\Tanner\Documents\git\ePC-SAFT fetch origin main --prune
git -C C:\Users\Tanner\Documents\git\ePC-SAFT pull --ff-only origin main
git -C C:\Users\Tanner\Documents\git\ePC-SAFT rev-parse --short=12 main
git -C C:\Users\Tanner\Documents\git\ePC-SAFT rev-parse --short=12 origin/main
```

The hashes must match.

Create the Codex worktree from updated `main` or `origin/main`. The Codex app may initially place the worktree on `main`; that is okay.

Assigned branch:

```text
codex/general-reactive-electrolyte-roadmap-setup
```

Before Goal Prep writes any files, the agent must switch/create this assigned branch using the branch bootstrap below:

```powershell
git fetch origin --prune
$current = (git branch --show-current).Trim()
if ($current -ne "codex/general-reactive-electrolyte-roadmap-setup") {
    git ls-remote --exit-code --heads origin codex/general-reactive-electrolyte-roadmap-setup
    if ($LASTEXITCODE -eq 0) {
        git switch codex/general-reactive-electrolyte-roadmap-setup 2>$null
        if ($LASTEXITCODE -ne 0) {
            git switch -c codex/general-reactive-electrolyte-roadmap-setup --track origin/codex/general-reactive-electrolyte-roadmap-setup
        }
    } else {
        git switch -c codex/general-reactive-electrolyte-roadmap-setup origin/main
    }
}
$current = (git branch --show-current).Trim()
if ($current -ne "codex/general-reactive-electrolyte-roadmap-setup") {
    Write-Error "BRANCH_BOOTSTRAP_FAILED expected=codex/general-reactive-electrolyte-roadmap-setup actual=$current"
    exit 1
}
git status --short
```


## First message: Goal Prep only

```text
[$goal-prep](C:\Users\Tanner\.codex\skills\goalbuddy\SKILL.md)

Initiate Goal Prep for this exact goal:

General reactive/electrolyte roadmap setup for ePC-SAFT agents

Use the local live GoalBuddy board in Codex.

Before Goal Prep writes files, run the branch bootstrap from this prompt. Do not write Goal Prep files on `main`.

Do not start `/goal` execution yet.
Do not edit package source files during Goal Prep.
Do not create local .worktrees/.
Do not ask me to choose among options.
Do not ask repeated confirmation questions.
Use the defaults in this prompt.

If a required value is missing, choose the safest conservative default and record it in docs/goals/<slug>/notes/assumptions.md.
If the missing value prevents safe work, stop with status BLOCKED_MISSING_INPUT.

If GoalBuddy or local board setup is unavailable, stop immediately and report setup is blocked.

Create:
- docs/goals/<slug>/goal.md
- docs/goals/<slug>/state.yaml
- docs/goals/<slug>/notes/

After branch bootstrap, current branch must be:

codex/general-reactive-electrolyte-roadmap-setup

This setup task is docs/issues/project/branch/prompt setup only.

Do not implement source code.

After Goal Prep, print the exact `/goal Follow docs/goals/<slug>/goal.md.` command once and stop.
Do not ask what to do next.
```

## Second message: Start goal

Use the exact `/goal Follow ...` path printed by Goal Prep, then provide:

```text
You are the setup agent for the next ePC-SAFT roadmap.

Your job is to create shared planning infrastructure. Do not implement package source code.

Use the attached setup bundle as the source of truth. Place files into the repository at the paths indicated inside the bundle.

Create or update these repo files:

docs/roadmaps/general_reactive_electrolyte_equilibrium_readiness.md
docs/roadmaps/agent_dependency_plan.md
docs/roadmaps/github_project_dashboard_plan.md
docs/roadmaps/branch_setup_commands.md
docs/roadmaps/agent_prompts/index.yaml
docs/roadmaps/agent_prompts/README.md
docs/roadmaps/agent_prompts/*.md
docs/roadmaps/issue_drafts/*.md
docs/roadmaps/watcher_templates/dependency_gate_template.yaml
docs/roadmaps/watcher_templates/bounded_dependency_watcher.ps1

Create GitHub issues A-M from docs/roadmaps/issue_drafts if GitHub issue tools or gh CLI are available.

If possible, create a GitHub Project named:

ePC-SAFT — General Reactive/Electrolyte Roadmap

with statuses:
Draft
Goal Prep Ready
Blocked on Dependency
Ready to Start
In Progress
PR Open
Needs Review
Merged
Stopped / Needs Handoff
Deferred

If GitHub Project creation is unavailable, do not improvise. Keep the project plan in docs and report that the project must be created manually.

Do not create placeholder implementation branches until this setup PR is merged into origin/main.

If you can merge the setup PR after validation:
1. Merge setup PR.
2. Fast-forward local main from origin/main.
3. Create placeholder branches from updated origin/main using docs/roadmaps/branch_setup_commands.md.
4. Push those branches.

If you cannot merge:
1. Open the setup PR.
2. Do not create implementation branches yet.
3. Report that branches must be created after setup PR merge.

Validation:
- uv run python scripts/validate_project.py docs
- git status --short
- verify files exist

PR body:
- Summary
- Roadmap files added
- Issue drafts or issue URLs
- Project board URL if created
- Branches created, or why not
- Validation
- Confirmation that no source code was modified

Do not close parent roadmap issues unless explicitly instructed.
```
