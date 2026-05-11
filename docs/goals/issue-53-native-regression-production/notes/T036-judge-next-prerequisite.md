# T036 next prerequisite decision

Date: 2026-05-11

## Result

Approved.

## Decision

The next honest prerequisite after the T035 blocker is
**concentration-standard-state native chemical-equilibrium derivatives**, not
full activity-coupled standard states.

## Why concentration comes first

It is the smaller mathematically valid next slice:

- the reaction residual still avoids activity coefficients;
- the only added coupling is through molar density in
  `species_activity = x_i * rho`;
- that means the derivative problem is:
  - implicit composition sensitivity through the inner chemical-equilibrium
    Jacobian, plus
  - density sensitivity with respect to composition and fitted parameters.

By contrast, activity-coupled standard states require all of that **plus**
activity-coefficient sensitivities, which is the broader Born/SSM+DS path the
user ultimately wants but is not the smallest next step.

## Immediate implication for issue #53

The next worker should:

1. add native concentration-standard-state `R_theta` support to chemical
   equilibrium;
2. keep activity-coupled standard states gated as `backend_unavailable`;
3. only after concentration works, reconsider whether the next derivative slice
   is direct activity-coupled standard states or a narrower Born-specific
   activity sensitivity path.
