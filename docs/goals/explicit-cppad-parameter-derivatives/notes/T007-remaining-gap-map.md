T007 remaining gap map

Dependency gate:
- Task A / PR #98 is merged into origin/main.
- Current branch is codex/cppad-explicit-parameter-derivatives.
- origin/main is an ancestor of HEAD.
- Corrected watcher settings require bounded auto-start, so implementation continued.

Implemented before this scout:
- Neutral binary k_ij pressure and ln-fugacity parameter derivatives through state/property result APIs.
- Pure-neutral m/sigma/epsilon pressure and ln-fugacity parameter derivatives.

Next safe implementation slice:
- Add residual-chemical-potential parameter derivatives for the same explicit CppAD paths.
- Generalize neutral binary pair-parameter CppAD support from k_ij to l_ij for nonassociating neutral binaries.

Not selected:
- k_hb_ij is associating-specific. It should remain unavailable in this tranche unless the association CppAD path is expanded and validated separately, because Task B must not tape solver loops or mislabel implicit association behavior as explicit coverage.
