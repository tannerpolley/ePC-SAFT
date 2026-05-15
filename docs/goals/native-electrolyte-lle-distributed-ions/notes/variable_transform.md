# Issue #116 T004 Variable Transform Note

Date: 2026-05-15

## Scope Completed In T004

T004 records native ownership of the electrolyte LLE basis and transformed feasible variables. This slice does not claim the Ceres residual solve or production Jacobian from T005.

## Native Basis Model

- The native solver now reports `basis_model = charge_neutral_salt_pair_coordinates`.
- Public species remain explicit in result phases and diagnostics.
- The formula basis is a coordinate system only. It reports:
  - neutral species indices
  - cation species indices
  - anion species indices
  - charged species indices
  - species charge vector
  - formula feed
  - salt-pair cation and anion indices
  - salt-pair stoichiometry
  - row-major salt-pair basis vectors in public species space

## Transform Guarantees Reported

- `phase_charge_enforced_by_basis = true`
- `material_balance_enforced_by_formula_transform = true`
- `formula_phase_positivity_enforced_by_transform = true`
- `explicit_public_species_reported = true`
- Accepted results report `accepted_transformed_variables_feasible = true`.

## Solved-State Diagnostics Added

Accepted electrolyte LLE results now report:

- `phase_charge_balance_feed`
- `phase_charge_balance_aq`
- `phase_charge_balance_org`
- `phase_charge_balance_max_abs`
- `accepted_beta_formula`
- `accepted_beta_org`
- `accepted_aq_formula`
- `accepted_org_formula`
- `material_balance_residual`

Failure diagnostics can report the same basis and best-candidate fields without marking them as accepted.

## Remaining Issue #116 Work

T005 must still replace the accepted production route with a Ceres trust-region residual solve and a real Jacobian path. The old predictive solver labels remain visible until that slice is complete.
