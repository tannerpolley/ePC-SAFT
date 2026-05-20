# Architecture Deepening Fullest

## Objective

Implement the five architecture-deepening candidates from the latest `improve-codebase-architecture` pass as far as the current repo can safely support: native route bridge, equilibrium route recipes, regression problem assembly, parameter source resolution, and EOS Harness state views.

## Original Request

`/goal Implement all the candidates to their fullest`

## Intake Summary

- Input shape: `existing_plan`
- Audience: ePC-SAFT package maintainers and downstream consumers.
- Authority: `requested`
- Proof type: `test`
- Completion proof: each candidate is implemented with code, tests, and docs where useful; any blocked portion has a concrete receipt explaining the blocker; focused validation and cleanup pass; tracked changes are committed locally unless the user redirects.
- Likely misfire: repeating the completed `architecture-deepening` and `architecture-deepening-followup` work, or adding shallow wrappers while leaving route, regression, parameter, and EOS behavior split across broad modules.
- Blind spots considered: native `_core` rebuild coordination, public behavior compatibility, no unsupported capability claims, derivative-path correctness, and avoiding broad rewrites that are not backed by executable evidence.

## Candidate Lanes

1. Native Route Bridge Module
2. Equilibrium Route Recipe Module
3. Regression Problem Assembly Module
4. Parameter Source Resolution Module
5. EOS Harness State View Module

## Run Command

```text
/goal Follow docs/goals/architecture-deepening-fullest/goal.md.
```

## Board Truth

Machine truth lives at `docs/goals/architecture-deepening-fullest/state.yaml`.
