# Issue Scope: Generic electrolyte LLE with distributed ions

Source: [GitHub issue #92](https://github.com/tannerpolley/ePC-SAFT/issues/92)

## Current scope

- Build electrolyte LLE on the generic LLE foundation.
- Use ion-based notation and distributed ions.
- Keep `epcsaft` general-purpose.
- Avoid application-specific public APIs.

## In scope

- ion-based species
- phase electroneutrality
- distributed ions
- mixed solvents
- mixed electrolytes
- charge-balance diagnostics
- Ascani 2022 benchmark attempt

## Out of scope

- lithium-extraction-specific public APIs
- forcing Ascani 2022 Case Study 2 into a fake pass if it is inconsistent

## Required policy

- No finite difference.
- CppAD for explicit algebraic derivatives.
- Analytic formulas where exact and validated.
- Implicit sensitivities for solved states.
- No backend_unavailable for required workflows.
- backend_unavailable only for explicitly out-of-scope workflows.

## Current gate status

- Dependency C / issue #86: closed.
- Dependency H / issue #91: closed.
- The bounded dependency watcher reported `GATE_PASS` on May 13, 2026.
- Implementation is unblocked on `codex/generic-electrolyte-lle`.

## Notes for the next `/goal`

- Start from the generic electrolyte LLE foundation, not an application-specific API layer.
- Preserve the issue scope above exactly.
- If Ascani 2022 Case Study 2 remains inconsistent, document it rather than forcing a pass.
