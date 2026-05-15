# T028 PR #126 Post-Merge Review Reopen

PR #126 merged at `869e3354ddc0b52075ddc9efe687b34d6aa98316` and closed issues #116 and #117, but the post-merge review identified unresolved completion gaps. This note reopens the GoalBuddy audit without erasing the earlier contradictory receipts.

## Reopened gates

1. `reaction_standard_states` must control the native coupled reactive phase-equilibrium reaction residual values and Jacobians. The merged residual path currently uses activity-form terms for all standard states.
2. The native reactive Ceres route must not return an accepted production result unless Ceres produced a usable converged solution and the final thermodynamic residual checks pass.
3. The reactive electrolyte benchmark must include charged species in the reaction itself. Ions that only participate in the phase split are not sufficient proof for issue #117.
4. The neutral reactive LLE benchmark must become source-backed Ascani-style esterification evidence, not a synthetic two-component model-consistent fixture.
5. The previous `T999` not-complete receipt and later `T027` complete receipt are contradictory. The final reopened audit must state which receipt is superseded by current evidence and must not claim completion until all reopened gates pass.

## Current status

`T027` is treated as superseded by this post-merge review. The board is active again, with `T029` through `T031` queued for source/test repair, benchmark replacement, and final audit.
