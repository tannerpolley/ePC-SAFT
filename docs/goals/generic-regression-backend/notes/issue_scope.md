# Issue Scope

Issue `#94` is the live work item for Task K: Generic regression row schema and native optimizer backend.

Current scope recorded from the issue:

- generic target-row compilation
- generic parameter maps
- Ceres preferred when it owns the native loop
- other native optimizers allowed only with analytic, CppAD, or implicit derivatives
- row diagnostics
- source summaries
- active bounds
- parameter movement

Out of scope:

- `fit_lithium_extraction_parameters`
- `fit_mea_absorption`
- finite difference
- Python-owned production objective loops

Current dependency state:

- Task B / PR `#102` is merged.
- Task C / PR `#104` is merged.
- Task E / PR `#101` is merged.

The issue itself remains open while the task branch proceeds through implementation and review.
