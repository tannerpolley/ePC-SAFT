# MEA_CO2_H2O_phase2

This is the single Phase 2 parameter artifact for the MEA-CO2-H2O true-species activity-based ePC-SAFT workflow.

## Source

The initial artifact is copied from:

`data/reference/epcsaft_datasets/MEA_CO2_H2O_ionic_fit/`

It is promoted into a Phase 2-owned directory so all Phase 2 pressure, speciation, residual, and comparison artifacts can name one parameter source.

## Policy

- Use this artifact for Phase 2 calculations and summaries.
- Do not update this artifact from failed or nonpromoted regression candidates.
- Keep MEA `f_solv = 1` unless direct MEA-H2O dielectric or activity evidence supports a different policy.
- Keep `HCO3- d_born = 3.0` and `CO3^2- d_born = 3.0` unless coupled pressure/speciation evidence and approval gates support movement.
- Do not use this artifact to claim a final Phase 3 globally regressed parameter set.
