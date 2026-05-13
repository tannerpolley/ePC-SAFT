# T005 Binary Derivative Blocker

## Blocker

Binary `k_ij` Ceres regression cannot be truthfully marked production in this slice without implementing missing native derivative coverage.

## Evidence

- `src/epcsaft/native/epcsaft_regression.cpp` still routes generic binary residual optimization through `GenericLeastSquaresFunctor::df`.
- That `df` implementation computes central finite differences, which issue #66 forbids for Ceres regression.
- `evaluate_generic_residuals_with_jacobian_cpp` currently throws `backend_unavailable: generic regression sensitivities are not implemented.`
- Binary VLE residuals depend on liquid and vapor `ln_fugacity_coefficient` values at density-root states. A real Ceres Jacobian for `k_ij` therefore needs supported parameter sensitivities through fugacity and solved density states, or an explicit analytic/CppAD/implicit implementation for that row family.

## Safe Next Options

1. Add a derivative coverage/capability gate that keeps binary Ceres rows out of production capability claims and returns `backend_unavailable`.
2. Implement the missing binary `k_ij` VLE derivative path with real analytic, CppAD, or implicit sensitivities before wiring Ceres to optimize it.
3. Keep the existing generic least-squares path as legacy/reference only; do not use it as issue #66 Ceres production support because its Jacobian is finite-difference based.
