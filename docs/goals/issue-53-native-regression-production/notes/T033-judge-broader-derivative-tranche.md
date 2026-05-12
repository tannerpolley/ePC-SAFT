# T033 broader derivative tranche decision

Date: 2026-05-11

## Result

Approved with a new architecture-first sequence.

## Decision

Prioritize the Born/speciation derivative substrate before reactive-electrolyte
bubble-pressure derivatives.

Reason:

- it matches the remaining Issue #53 regression test scope the user clarified:
  `d_born`, `f_solv`, and generic binary parameters such as `k_ij`;
- it is smaller than the bubble-pressure path, which still needs a coupled
  log-pressure plus vapor-composition implicit sensitivity system;
- it directly serves the next honest native thermodynamic row expansion for
  reactive-speciation regression.

## First bounded implementation slice

Land the shared derivative-policy scaffolding that the broader tranche needs:

- move native backend-unavailable debug gating to a shared helper under
  `src/epcsaft/native/autodiff/`;
- start scalar-templating cleanup in the Born path where double and autodiff
  logic are currently duplicated;
- do not widen production capability claims yet.

This slice is structural, not celebratory. It is allowed because it removes
duplication and narrows the next implementation surface without pretending that
`r_theta` support already exists for `d_born`, `f_solv`, or `k_ij`.

## Next implementation target after this slice

Add the first real non-logK native parameter-sensitivity path for
`reactive_speciation` rows:

1. define parameter-to-thermo sensitivity ownership for `d_born`, `f_solv`, and
   `k_ij`;
2. thread those sensitivities into the residual/activity/fugacity path;
3. extend implicit speciation `R_theta` assembly;
4. keep bubble-pressure rows gated as `backend_unavailable` until their own
   implicit derivative system exists.

