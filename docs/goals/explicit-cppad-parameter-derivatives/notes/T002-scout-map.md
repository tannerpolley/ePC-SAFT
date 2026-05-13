# T002 Scout Map

## Verified Facts

- `src/epcsaft/epcsaft.py` owns public state derivative result methods:
  - `pressure_parameter_derivative_result` is still `backend_unavailable`.
  - `ln_fugacity_parameter_derivative_result` supports analytic Born SSM+DS derivatives when `born_ssmds_liquid_derivatives()` is supported, otherwise it reports unavailable.
  - `activity_parameter_derivative_result` follows the same Born SSM+DS path for natural-log activity coefficients.
  - `relative_permittivity_parameter_derivative_result` already supports analytic linear mole-fraction dielectric parameter derivatives.
- `ePCSAFTState.derivative_coverage_matrix()` currently omits explicit pressure/fugacity/activity parameter rows except for the Born SSM+DS row.
- `src/epcsaft/native/epcsaft_ares.cpp` already contains `neutral_binary_kij_phase_derivatives_cpp`, a CppAD-backed evaluator that computes:
  - `dpdk` for pressure with respect to a binary `k_ij` entry at fixed density.
  - `dlnphi_dk_fixed_rho` and `dlnphi_dk_total`.
  - The routine is intentionally limited to exactly two neutral, nonassociating components with a dense `k_ij` matrix.
- `src/epcsaft/native/epcsaft_regression.cpp` already consumes `neutral_binary_kij_phase_derivatives_cpp` for Ceres binary VLE `k_ij` residual Jacobians, so the derivative path is not speculative.
- Existing targeted tests already cover result-shape contracts and no-finite-difference assertions:
  - `tests/native/test_cppad_pressure_derivatives.py`
  - `tests/native/test_cppad_fugacity_derivatives.py`
  - `tests/native/test_cppad_activity_derivatives.py`
  - `tests/native/test_cppad_relative_permittivity_derivatives.py`
  - `tests/native/test_derivative_coverage_matrix.py`

## Coverage Gaps

- Pressure parameter derivatives are not wired through the state/property API.
- Fugacity parameter derivatives do not expose the existing neutral binary `k_ij` path.
- Runtime coverage rows do not enumerate pressure/fugacity/activity/relative-permittivity parameter result rows from the state matrix.
- General `m/sigma/epsilon`, `l_ij`, and `k_hb_ij` state/property parameter derivatives remain broader than the smallest existing verified native path.

## Candidate Worker Slice

Expose the existing neutral binary `k_ij` CppAD derivative path through generic state/property result methods:

- Add a pybind helper for neutral binary `k_ij` property derivatives.
- Wire `pressure_parameter_derivative_result()` and `ln_fugacity_parameter_derivative_result()` to return one generic parameter column for `k_ij:<left>:<right>` when the state is an eligible neutral binary case.
- Keep all other cases as honest `backend_unavailable`.
- Add targeted tests for pressure and fugacity parameter derivatives, and add coverage-matrix rows for parameter derivatives.
