# ePC-SAFT Roadmap Worktree Repair Status

Generated: 2026-05-13 America/Denver

## Root Checkout

- Root checkout: `C:\Users\Tanner\Documents\git\ePC-SAFT`
- Branch: `codex/worktree-coordinator-repair`
- Root `HEAD`: `660ca5884821e9b693478bb9ca9bbdd8abe934b0`
- `origin/main`: `660ca5884821e9b693478bb9ca9bbdd8abe934b0`
- PR #98: merged into `main` at `660ca5884821e9b693478bb9ca9bbdd8abe934b0`
- PR #99: open, mergeable, checks green, truthfulness acceptable for convention-layer scope

## Worktrees Found

| Worktree | Repo path | Task | Branch | Status class | PR |
| --- | --- | --- | --- | --- | --- |
| `4744` | `C:\Users\Tanner\.codex\worktrees\4744` | unknown | no checkout | `UNKNOWN_NEEDS_HUMAN` | n/a |
| `4a52` | `C:\Users\Tanner\.codex\worktrees\4a52\ePC-SAFT` | A | `codex/backend-coverage-hard-gate` | `MERGED_CLEANUP_NEEDED` | PR #98 merged |
| `964e` | `C:\Users\Tanner\.codex\worktrees\964e\ePC-SAFT` | B | `codex/cppad-explicit-parameter-derivatives` | `NEEDS_SCOPE_REPAIR` | none found |
| `b752` | `C:\Users\Tanner\.codex\worktrees\b752\ePC-SAFT` | C | `codex/generic-implicit-sensitivity-framework` | `OK_WAITING_WATCHER` | none found |
| `6443` | `C:\Users\Tanner\.codex\worktrees\6443\ePC-SAFT` | D | `codex/reaction-constant-conventions` | `HAS_OPEN_PR_REVIEW` | PR #99 open |
| `e188` | `C:\Users\Tanner\.codex\worktrees\e188\ePC-SAFT` | E | `codex/generic-target-row-schema` | `OK_RUNNING` locally complete, PR needed | none found |
| `f5b0` | `C:\Users\Tanner\.codex\worktrees\f5b0\ePC-SAFT` | policy repair | `codex/autorun-watcher-prompt-fix` | local policy patch copied into coordinator branch | none found |

No Codex transcript/session/jsonl files were found by the targeted worktree scan. Status above is based on Git state, GoalBuddy control files, local notes, and live GitHub PR data.

## Task Mapping A-M

| Task | Branch | Issue | Current local state |
| --- | --- | --- | --- |
| A | `codex/backend-coverage-hard-gate` | #84 | PR #98 merged; worktree can be cleaned after dependent agents are stable. |
| B | `codex/cppad-explicit-parameter-derivatives` | #85 | Active local implementation, no PR. Needs PR or `BLOCKED_SCOPE_GAP`; `k_hb_ij` must be represented in coverage/capability. |
| C | `codex/generic-implicit-sensitivity-framework` | #86 | Waiting on B with bounded watcher; local GoalBuddy files repaired. |
| D | `codex/reaction-constant-conventions` | #87 | PR #99 open and mergeable; truthfulness review passed for convention-layer scope. |
| E | `codex/generic-target-row-schema` | #88 | Local board says complete and validation passed; no PR found, so owning agent should open PR. |
| F | `codex/generic-activity-speciation` | #89 | Future task; repo prompt policy now uses bounded watcher auto-run. |
| G | `codex/generic-vle-fugacity-equilibrium` | #90 | Future task; repo prompt policy now uses bounded watcher auto-run. |
| H | `codex/generic-non-electrolyte-lle` | #91 | Future task; repo prompt policy now uses bounded watcher auto-run. |
| I | `codex/generic-electrolyte-lle` | #92 | Future task; repo prompt policy now uses bounded watcher auto-run. |
| J | `codex/generic-reactive-lle` | #93 | Future task; repo prompt policy now uses bounded watcher auto-run. |
| K | `codex/generic-regression-backend` | #94 | Future task; repo prompt policy now uses bounded watcher auto-run. |
| L | `codex/literature-benchmark-suite` | #95 | Future task; repo prompt policy now uses bounded watcher auto-run. |
| M | `codex/downstream-integration-smokes` | #96 | Future task; repo prompt policy now uses bounded watcher auto-run. |

## Current B/C/D/E Status

### Task B

Status: `NEEDS_SCOPE_REPAIR`.

The Task B worktree has implementation edits and passing local validation recorded in its GoalBuddy state, but no PR exists for `codex/cppad-explicit-parameter-derivatives`. The local board is now repaired so active task `T009` must finish with either:

- `PR_OPENED`, with truthful coverage/capability language, or
- `BLOCKED_SCOPE_GAP`, with exact file/function, derivative, parameter family, reason, and future owner.

`k_hb_ij` must not be hidden as a vague limitation. If full active-association derivatives require implicit association site-fraction sensitivity, Task B must mark that as `blocker_requires_implicit_association_sensitivity`, include coverage/capability rows, add a no-silent-omission/no-overclaim test, and identify Task C as owner.

### Task C

Status: `OK_WAITING_WATCHER`.

Task C remains blocked until Task B merges. Its local GoalBuddy root-file issue was repaired: dependency gate and watcher script now live under `notes/`, and the state checker passes. The watcher should continue waiting for Task B, then update from `origin/main` and start implementation automatically.

### Task D

Status: `HAS_OPEN_PR_REVIEW`.

PR #99 is open and mergeable: `https://github.com/tannerpolley/ePC-SAFT/pull/99`.

Truthfulness review passed for the stated scope:

- `ReactionConstantConvention` is generic.
- No application-specific API was found.
- No finite-difference route was found.
- Fitted constants are explicit, not default.
- Molality/apparent constants are defined at the convention layer and documented as `backend_unavailable` for unsupported native residual routes.
- Checks are green.

Recommendation: safe to merge after normal Tanner/repo authorization. This coordinator did not auto-merge it.

### Task E

Status: `OK_RUNNING` locally complete, but PR is missing.

The Task E board says the tranche completed and docs validation passed. No PR was found for `codex/generic-target-row-schema`. The owning Task E agent should open a focused PR from that branch, or write `BLOCKED_SCOPE_GAP` if it cannot.

## Local GoalBuddy Files Repaired

Moved root-level gates into `notes/`, created or moved watcher scripts, updated path references, and ran the GoalBuddy checker successfully for:

- `C:\Users\Tanner\.codex\worktrees\4a52\ePC-SAFT\docs\goals\derivative-backend-coverage-hard-gate`
- `C:\Users\Tanner\.codex\worktrees\964e\ePC-SAFT\docs\goals\explicit-cppad-parameter-derivatives`
- `C:\Users\Tanner\.codex\worktrees\b752\ePC-SAFT\docs\goals\generic-implicit-sensitivity-framework`
- `C:\Users\Tanner\.codex\worktrees\6443\ePC-SAFT\docs\goals\reaction-constant-conventions`
- `C:\Users\Tanner\.codex\worktrees\e188\ePC-SAFT\docs\goals\generic-target-row-schema`

Task B also received `notes/repair_note.md` and a state/goal update requiring `PR_OPENED` or `BLOCKED_SCOPE_GAP`.

## GoalBuddy UI And Checker Health

The Codex app UI was inspected from the visible Codex window. No separate Goal tab was visible above the prompt window. The right-side panel showed normal branch details and an artifact link for `goal.md`, but not a Goal tab. Treat this as a UI integration limitation unless the app later exposes the tab.

The dependency gate remains file-based and does not depend on the UI tab:

- dependency issue closed
- dependency PR merged
- `origin/main` contains the dependency merge commit
- assigned branch is current
- branch updates cleanly from `origin/main`
- `docs/goals/<slug>/notes/dependency_gate.yaml` passes

| Worktree | Task | Root entries after repair | Checker | Local board URL | `/goal` state activity | Checker conflict |
| --- | --- | --- | --- | --- | --- | --- |
| `4a52` | A | `goal.md`, `state.yaml`, `notes/` | pass | stale/not live | done, no active task | none |
| `964e` | B | `goal.md`, `state.yaml`, `notes/` | pass | stale/not live | active task `T009` | none |
| `b752` | C | `goal.md`, `state.yaml`, `notes/` | pass | stale/not live | active task `T001`, waiting on B | none |
| `6443` | D | `goal.md`, `state.yaml`, `notes/` | pass | not configured | done, PR #99 open | none |
| `e188` | E | `goal.md`, `state.yaml`, `notes/` | pass | not configured | done, PR missing | none |
| `f5b0` | policy | `goal.md`, `state.yaml`, `notes/` | pass | stale/not live | done, policy patch copied forward | none |

The coordinator `/goal` run is changing state normally in `docs/goals/worktree-coordinator-repair/state.yaml`; the missing UI tab did not block file-based GoalBuddy execution.

## Repo Files Patched

- `docs/COPY_PASTE_AGENT_PROMPTS.md`
- `docs/roadmaps/agent_prompts/index.yaml`
- `docs/roadmaps/agent_prompts/README.md`
- `docs/roadmaps/agent_prompts/A_derivative_backend_completion_audit_and_coverage_matrix_hard_gate.md`
- `docs/roadmaps/agent_prompts/B_explicit_cppad_parameter_derivatives_for_eos_property_apis.md`
- `docs/roadmaps/agent_prompts/C_generic_implicit_sensitivity_framework_for_solved_states.md`
- `docs/roadmaps/agent_prompts/D_general_reaction_and_equilibrium_constant_convention_layer.md`
- `docs/roadmaps/agent_prompts/E_generic_target_row_and_dataset_schema.md`
- `docs/roadmaps/agent_prompts/F_generic_speciation_solver_using_epc_saft_activities.md`
- `docs/roadmaps/agent_prompts/G_generic_vle_fugacity_equilibrium_solver_for_volatile_neutral_species.md`
- `docs/roadmaps/agent_prompts/H_generic_non_electrolyte_lle_benchmark_and_solver_hardening.md`
- `docs/roadmaps/agent_prompts/I_generic_electrolyte_lle_with_distributed_ions.md`
- `docs/roadmaps/agent_prompts/J_generic_reactive_lle_and_chemical_phase_equilibrium.md`
- `docs/roadmaps/agent_prompts/K_generic_regression_row_schema_and_native_optimizer_backend.md`
- `docs/roadmaps/agent_prompts/L_literature_benchmark_suite.md`
- `docs/roadmaps/agent_prompts/M_downstream_integration_smoke_tests.md`
- `docs/roadmaps/watcher_templates/dependency_gate_template.yaml`
- `docs/roadmaps/watcher_templates/bounded_dependency_watcher.ps1`
- `docs/roadmaps/autorun_watcher_policy.md`
- `docs/roadmaps/agent_scope_contract_policy.md`
- `docs/roadmaps/running_agent_repair_prompts.md`
- `docs/roadmaps/worktree_repair_status.md`

Validation:

```text
uv run python scripts/validate_project.py docs
```

Result: pass.

## Next Actions For Tanner

1. Merge PR #99 if you accept the convention-layer scope and `backend_unavailable` limitation for molality/apparent native residual routes.
2. Restart or message Task B with the direct repair prompt in `docs/roadmaps/running_agent_repair_prompts.md`; it must open a PR or stop with `BLOCKED_SCOPE_GAP`.
3. Leave Task C watcher waiting until Task B merges.
4. Ask the Task E agent to open a focused PR from `codex/generic-target-row-schema`.
5. After B, D, and E are merged, clean up merged/stale worktrees including A and any empty/unknown `4744` directory if no active Codex thread owns it.
