# Rezaee 2026 calibrated native Ipopt attempt

This diagnostic uses the local copy of the downstream paper-constants refit calibration.
It does not promote the Rezaee lane to direct closure.

## Source-phase residuals

- rows: 26
- median abs ln residual: 334.451
- mean abs ln residual: 333.201
- max abs ln residual: 345.513

## Public route attempt

- experiment_no: 1
- accepted: False
- error_type: SolutionError
- solver_backend: None
- derivative_backend: None
- density_backend: None
- source_charge: -1.499996e-07
- neutralized_charge: 2.168404e-19
- neutralization: 1.499996e-07 added to H+

### Error

('Native reactive LLE route was rejected.', {'accepted': False, 'rejection_reason': '', 'derivative_backend': 'cppad_implicit', 'density_backend': 'liquid_pressure_root', 'phase_count': 0, 'species_count': 0, 'balance_row_count': 0, 'reaction_count': 0, 'conserved_balance_norm': 0.0, 'charge_balance_norm': 0.0, 'pressure_consistency_norm': 0.0, 'phase_equilibrium_norm': 0.0, 'reaction_stationarity_norm': 0.0, 'phase_distance': 0.0, 'objective': 0.0, 'standard_mu_rt': [], 'constraints': [], 'reaction_stationarity_residuals': [], 'phase_equilibrium_residuals': [], 'phase_charge_residuals': [], 'phase_amount_totals': [], 'phase_volumes': [], 'phase_densities': [], 'phase_compositions': [], 'phase_ln_fugacity_coefficients': [], 'route_status': 'solver_rejected', 'solver_status': 'error_in_step_computation'})
