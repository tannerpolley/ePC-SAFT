# T025 Judge: Next Repair Decision

## Decision

repair_next

## Reason

The blocker is localized enough for one more bounded Worker task before declaring a broader design handoff. The neutral phase-state sensitivity path matches projected perturbations, while the active association/ionic path does not. That points the next slice at active association plus ionic/Born terms in `phase_state_ln_fugacity_composition_sensitivity_cpp(...)`, not at the Ceres route itself.

## Next Worker Scope

Allowed files:

- `src/epcsaft/native/epcsaft_ares.cpp`
- `src/epcsaft/native/epcsaft_core_internal.h`
- `src/epcsaft/native/epcsaft_electrolyte.h`
- `src/epcsaft/native/epcsaft_equilibrium.cpp`
- `tests/native/cppad/test_phase_state_sensitivities.py`
- `tests/native/equilibrium/test_electrolyte_lle_ceres_solver.py`
- `tests/native/equilibrium/test_electrolyte_lle_residual_jacobian.py`
- `docs/goals/native-lle-reactive-production-solvers/notes/**`

Required proof:

- projected solved-state perturbation test for a neutral pressure state;
- projected solved-state perturbation test for the active association/ionic Ascani distributed-ion state;
- accepted Ceres solve no longer reports `local_residual_slope` if the analytic repair succeeds.

Stop condition:

- if the active association/ionic mismatch requires redesigning residual chemical potential or fugacity coefficient derivative ownership across contribution files, stop with a broader design blocker instead of widening the Worker slice.
