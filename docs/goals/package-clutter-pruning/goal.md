# Package Clutter Pruning Goal

Source authority: GitHub issue #120, "Prune generated artifacts, runtime clutter, and staged/optional API noise before vertical implementation work"

Issue URL: https://github.com/tannerpolley/ePC-SAFT/issues/120

## Outcome

Make the package tree easier for future vertical implementation agents to reason about before #114-#119 work starts. This goal is cleanup, relocation, de-duplication, and guardrail hardening only.

## Non-Goals

- Do not change EOS equations, association equations, Born / SSM / DS equations, Debye-Huckel equations, mixing rules, parameter values, Ceres objective math, equilibrium residual math, regression targets, benchmark tolerances, or public scientific outputs.
- Do not implement #114, #115, #116, #117, #118, or #119.
- Do not close by inventory, diagnostics, staged workflows, synthetic proof, or limitation text.
- Do not broaden public capability claims.

## Required Phase Order

1. Phase 0 - Intake and safety baseline.
2. Phase 1 - Remove generated and compiled artifacts from source tracking.
3. Phase 2 - Move benchmark execution out of the runtime package.
4. Phase 3 - Move analysis and confidence scoring out of runtime diagnostics.
5. Phase 4 - Clarify native-regression scaffolding by removing placeholders or making the folder real.
6. Phase 5 - Quarantine staged and optional top-level helpers.
7. Phase 6 - Create the module ownership and safe split plan.
8. Phase 7 - Audit and reduce broad exception handling where safe.
9. Phase 8 - Rewrite ambiguous compatibility, staged, and debug-route language.
10. Phase 9 - Audit and tighten the public API surface.
11. Final audit - Run required validation, prepare PR evidence, and keep the issue open or draft if any phase is incomplete.

## Completion Proof

The goal is complete only when all #120 deliverables exist, required validation passes, the generated-artifact scan reports zero tracked generated or compiled artifacts or justified exceptions, public/runtime ownership is explicit, and the final PR evidence confirms no scientific equations, solver algorithms, regression math, benchmark targets, or scientific outputs changed.

## Starter Command

`/goal Follow docs/goals/package-clutter-pruning/goal.md.`
