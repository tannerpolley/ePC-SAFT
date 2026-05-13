# T001 Worktree Inventory

Timestamp: 2026-05-13 America/Denver

## Root State

- Root checkout: `C:\Users\Tanner\Documents\git\ePC-SAFT`
- `HEAD`: `660ca5884821e9b693478bb9ca9bbdd8abe934b0`
- `origin/main`: `660ca5884821e9b693478bb9ca9bbdd8abe934b0`
- PR #98: merged, merge commit `660ca5884821e9b693478bb9ca9bbdd8abe934b0`
- PR #99: open, mergeable, head `codex/reaction-constant-conventions`

## Worktrees Found

| ID | Path | Task | Branch | Goal slug | Current status |
| --- | --- | --- | --- | --- | --- |
| 4744 | `C:\Users\Tanner\.codex\worktrees\4744` | unknown | no repo checkout | n/a | `UNKNOWN_NEEDS_HUMAN` |
| 4a52 | `C:\Users\Tanner\.codex\worktrees\4a52\ePC-SAFT` | A | `codex/backend-coverage-hard-gate` | `derivative-backend-coverage-hard-gate` | `MERGED_CLEANUP_NEEDED` |
| 964e | `C:\Users\Tanner\.codex\worktrees\964e\ePC-SAFT` | B | `codex/cppad-explicit-parameter-derivatives` | `explicit-cppad-parameter-derivatives` | `NEEDS_SCOPE_REPAIR` |
| b752 | `C:\Users\Tanner\.codex\worktrees\b752\ePC-SAFT` | C | `codex/generic-implicit-sensitivity-framework` | `generic-implicit-sensitivity-framework` | `NEEDS_LOCAL_GOAL_REPAIR` and waiting on B |
| 6443 | `C:\Users\Tanner\.codex\worktrees\6443\ePC-SAFT` | D | `codex/reaction-constant-conventions` | `reaction-constant-conventions` | `HAS_OPEN_PR_REVIEW` |
| e188 | `C:\Users\Tanner\.codex\worktrees\e188\ePC-SAFT` | E | `codex/generic-target-row-schema` | `generic-target-row-schema` | `OK_RUNNING` locally complete but no open PR found |
| f5b0 | `C:\Users\Tanner\.codex\worktrees\f5b0\ePC-SAFT` | policy | `codex/autorun-watcher-prompt-fix` | `fix-a-m-roadmap-prompts-bounded-watcher-auto-run` | docs policy patch exists locally, no open PR found |

## GoalBuddy Repair Findings

- Active roadmap goal directories with root-level `dependency_gate.yaml`: A, B, C, D, E.
- Active roadmap goal directories with `notes/dependency_gate.yaml`: none.
- Active roadmap goal directories with root-level `watch_dependency.ps1`: C only.
- Active roadmap goal directories with `notes/watch_dependency.ps1`: none.
- B/C/D/E goal text or state still mention `PREPARED_READY`; C specifically records the checker rejecting root `dependency_gate.yaml` and `watch_dependency.ps1`.
- Targeted transcript/session scan found no transcript/session/jsonl files in these worktree directories. No transcript content was used.

## Task-Specific Facts

- Task A: PR #98 is merged; worktree has only GoalBuddy state/board artifacts dirty.
- Task B: implementation edits are present; `state.yaml` active task is `T009`. The final audit says the current ending is not complete because no PR was opened and `k_hb_ij` remains unavailable. This matches the handoff's unacceptable failure mode.
- Task C: active watcher task is waiting for B. It has bounded watcher logic but root-level gate/script paths break GoalBuddy v2 checker compatibility.
- Task D: PR #99 is open and mergeable. PR body explicitly says molality/apparent constants raise `backend_unavailable` rather than fake support. Tests/checks are green.
- Task E: state says complete after accepting PR #98 as dependency; package source edits are present in the E worktree. No open PR for `codex/generic-target-row-schema` was found.
- Policy worktree: f5b0 contains the repo prompt-policy patch requested by the handoff, but the main checkout still has stale prompt policy values.

## Candidate Repairs

- Repair local GoalBuddy control files for B/C/D/E and optionally A: move/copy root `dependency_gate.yaml` into `notes/`, move/copy `watch_dependency.ps1` into `notes/`, and update references to `notes/` paths.
- For Task B, write an explicit local repair note and state entry requiring `PR_OPENED` or `BLOCKED_SCOPE_GAP`, with `k_hb_ij` represented in coverage/capability and Task C named as implicit-association owner.
- Patch main repo prompt policy from the f5b0 docs-only worktree or equivalent edits on a focused branch.
- Create `docs/roadmaps/worktree_repair_status.md`.
