# Phase 8 Ambiguous Wording Audit

## Package-Facing Rewrites

- Rephrased IPOPT docs as an explicit opt-in residual-minimization bridge, not a default production solver or full constrained thermodynamic NLP.
- Rephrased reactive staged workflow docs as sequential coupling where prose could imply a full coupled reactive flash.
- Replaced broad fallback wording in user-facing guidance with explicit route, substitute, conservative-path, or compatibility-shim language where behavior was not being renamed.
- Updated release, installation, architecture, cookbook, downstream install, diagnostics, presentation, and current-state planning text to keep the same capability claims with sharper role labels.

## Remaining Search Matches

- Public API names and diagnostic keys were intentionally left unchanged during this audit, including staged reactive-equilibrium route names, neutral candidate-selection flags, density closure candidate-selection flags, and `experimental_coupled_density_lle`.
- Later pruning work removed the neutral fallback flag and renamed density closure diagnostics to `density_best_candidate_refinement_used`, `density_best_candidate_rejection_reason`, and `density_warm_start_rejections`.
- Source-literature and validation-data labels named experimental remain unchanged because they refer to paper or digitized data, not package maturity.
- Git staging language in the LaTeX mirror helper remains unchanged because it refers to Git index state.
- FULL_ROADMAP warnings remain unchanged where they explicitly forbid closure on staged-only or diagnostic-only evidence.
- Compatibility-shim wording remains where the text explicitly describes a retained legacy import or dispatcher path.

## Boundary

No behavior, public API names, diagnostic keys, solver routes, or capability
truth values were changed in this wording phase.
