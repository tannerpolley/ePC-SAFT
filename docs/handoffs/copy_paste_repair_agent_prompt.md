# Copy-Paste Prompt for PR #126 Repair Agent

Thread name: ePC-SAFT PR126 Repair Audit and Completion

You are the repair agent for `tannerpolley/ePC-SAFT` after PR #126 merged incomplete work for issues #116 and #117.

Read the attached file:

```text
pr126_repair_audit_and_completion_contract.md
```

completely before planning or editing.

Your job is not to defend PR #126. Your job is to audit current main, fix the incomplete implementation, and prove the fixes with strict tests.

Start on a fresh branch from current `origin/main`.

Before editing source, create:

```text
docs/goals/pr126-repair-audit/notes/current_state_audit.md
```

Run and record the audit commands from the attached document.

Then implement all required repairs:

1. Apply reaction standard-state conventions in native reactive residuals.
2. Gate reactive Ceres acceptance on solver usability and physical residual gates.
3. Replace accepted neutral LLE hand-coded minimizer route with the shared derivative-backed residual solver.
4. Make failed-solve diagnostics honest.
5. Replace weak benchmark smokes with source-backed benchmark proof or leave explicit open follow-up without claiming completion.
6. Add a shared native residual-solver abstraction and route neutral LLE, electrolyte LLE, and reactive LLE through it.
7. Add the tests listed in the audit document.

Do not ask whether to reduce scope. Do not close with docs-only, diagnostics-only, smoke-only, or future-work completion.

If any stop rule in the attached document is triggered, stop, write the stopped-state note, and keep the PR draft/open.

The PR may be ready only when every item in the completion checklist passes.
