# Rezaee 2026 calibrated native Ipopt attempt

This diagnostic uses `current_package_paper_constants_refit_candidate`.
It does not promote the Rezaee lane to direct closure.

## Source-phase residuals

Separate-phase replay convention:
- rows: 26
- median abs ln residual: 1.13388
- mean abs ln residual: 1.12034
- max abs ln residual: 2.58537

Native combined-mixture phase-tagged residual:
- rows: 26
- median abs ln residual: 367.962
- mean abs ln residual: 368.01
- max abs ln residual: 385.856

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

('Native reactive LLE route was rejected.', {'accepted': False, 'rejection_reason': 'charge_balance', 'derivative_backend': 'cppad_implicit', 'density_backend': 'liquid_pressure_root', 'phase_count': 2, 'species_count': 11, 'balance_row_count': 8, 'reaction_count': 2, 'conserved_balance_norm': 1.1435297153639112e-14, 'charge_balance_norm': 0.4796437339413225, 'pressure_consistency_norm': 0.0, 'phase_equilibrium_norm': 0.0, 'reaction_stationarity_norm': 185.15270534675156, 'phase_distance': 0.9600822411811756, 'objective': 28126.470507010275, 'standard_mu_rt': [], 'constraints': [1.1435297153639112e-14, 1.0195837230053684e-15, 1.1969591984239969e-15, 0.0, 9.75781955236954e-19, 2.9947908641359856e-18, 1.27675647831893e-15, -6.938893903907228e-18, 148.22680816467727, 185.15270534675156, -0.4796437339413225, 6.245373013030353e-14], 'reaction_stationarity_residuals': [148.22680816467727, 185.15270534675156], 'phase_equilibrium_residuals': [], 'phase_charge_residuals': [-0.4796437339413225, 6.245373013030353e-14], 'phase_amount_totals': [0.5099498774871067, 0.49999992500029167], 'phase_volumes': [6.40762383930233e-05, 0.00012950314825630187], 'phase_densities': [7958.486488536299, 3860.9094198291873], 'phase_compositions': [[1.9609764051860798e-14, 0.006543089475139049, 0.02164367925817333, 0.010202933025535453, 0.0015274984970027974, 0.9600822411811756, 5.585629541977556e-07, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8811999524799697, 0.11880004751999035, 1.9999996571518496e-14, 1.99999965672325e-14]], 'phase_ln_fugacity_coefficients': [[-103.53511631410706, -180.3612225537372, -181.2581951297914, -184.21107699874224, -170.94272850223172, -169.04982330253858, -143.4961857147545, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -53.793320192238916, -41.91297841941582, -113.06306001938843, -83.80132106126852]], 'route_status': 'postsolve_rejected', 'solver_status': 'success'})
