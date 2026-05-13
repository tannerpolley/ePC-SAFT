# Running Agent Repair Prompts

Use these prompts only for already-running roadmap agents that are stuck on stale GoalBuddy automation, root-level dependency gate files, or silent scope narrowing.

## Generic Bounded Watcher Repair

```text
Correction: repair this local GoalBuddy task for bounded watcher auto-run and checker compatibility.

Do not restart the task unless you are blocked.
Do not ask repeated confirmation questions.
Do not ask me to choose among options.

Move or copy any root-level dependency gate:

docs/goals/<slug>/notes/dependency_gate.yaml

to:

docs/goals/<slug>/notes/dependency_gate.yaml

Create or update:

docs/goals/<slug>/notes/watch_dependency.ps1

Update local goal.md and state.yaml references to use:

docs/goals/<slug>/notes/dependency_gate.yaml
docs/goals/<slug>/notes/watch_dependency.ps1

Set:

watcher_mode: bounded
auto_start_after_gate: true
poll_interval_seconds: 120
max_wait_minutes: 480

Do not stop at PREPARED_READY merely because auto_start_after_gate is false.

New behavior:
- If dependencies are satisfied, update/fast-forward or rebase branch onto origin/main and start implementation without asking.
- If dependencies are not satisfied, run the bounded watcher.
- If timeout, write PREPARED_WAITING and stop.
- If rebase/merge conflicts occur, stop with BLOCKED_REBASE_CONFLICT.
- After implementation is complete, open a focused draft PR if needed, self-review it against origin/main, mark ready only after validation passes, wait for checks, and merge without another yes only when checks pass, the PR is mergeable, and the final GoalBuddy audit says full_outcome_complete: true.
- After merge, delete both the remote branch and the local task branch. If cleanup cannot be completed safely, stop with BLOCKED_REMOTE_BRANCH_DELETE or BLOCKED_LOCAL_BRANCH_DELETE and record the exact branch/path/next command.

If git rebase is blocked by tool approval policy and this branch has no unique commits, use git merge --ff-only origin/main.

Do not create local .worktrees/.
Do not write files on main.
Do not use PR #56 as a base.
Do not add application-specific public APIs.
No finite difference.
```

## Direct Task B Scope Repair

```text
Correction: Task B narrowed scope too aggressively and must repair its scope contract.

Do not close Task B and do not stop with only a recorded limitation.

The statement:

"k_hb_ij remains a recorded limitation because it is associating-specific and needs a separate validated association CppAD path rather than being mislabeled as explicit support. A PR has not been opened yet."

is not an acceptable task ending.

Technical correction:
- k_hb_ij is a real ePC-SAFT parameter family.
- Full active-association property/regression derivatives through k_hb_ij may require implicit association site-fraction sensitivities, which are owned by Task C.
- But Task B still owns explicit parameter derivative contracts, coverage/capability rows, and tests preventing silent omission or overclaiming.

Required Task B output:
1. Implement or verify explicit algebraic parameter derivative support/contracts for m, sigma, epsilon, k_ij, l_ij, d_born, f_solv, dielectric/relative-permittivity parameters, explicit property APIs where expressible through available blocks, and regression API derivative lookup by parameter name.
2. Include k_hb_ij in coverage/capability.
3. If full k_hb_ij property/regression derivatives require association implicit sensitivities, classify as blocker_requires_implicit_association_sensitivity, add/update coverage row, add/update capability row, add/update a no-silent-omission/no-overclaim test, and identify Task C as owner.
4. Open a focused PR.

Do not use finite difference.
Do not add application-specific public APIs.
Do not create local .worktrees/.
Do not write files on main.
Do not leave without PR unless you hit BLOCKED_SCOPE_GAP.

If blocked, stop with:

BLOCKED_SCOPE_GAP

and list exact missing derivative/function, file/function, parameter family, why it cannot be completed now, and future owner.
```

