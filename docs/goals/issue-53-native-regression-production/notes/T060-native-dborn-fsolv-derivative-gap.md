# T060 Native d_born/f_solv Derivative Gap

Date: 2026-05-11  
Goal lane: `issue-53-native-regression-production`
Thread note: add backend-unavailable eradication contract for active tranche.

## Hard Instruction Update

This tranche now has a hard non-negotiable contract:

- Do not keep, add, or preserve any backend-unavailable derivative path in core regression, Eq/thermo derivative routing, or debug surfaces for this issue-53 lane.
- The lane is an absolute zero-tolerance policy for backend-unavailable usage as compatibility/fallback/debug/benchmark-emitted behavior for supported paths. Use `backend_unavailable` when unsupported.
- No debug gate should be the mechanism for accepting Backend unavailable in this lane.

## Immediate scope for this tranche

- `d_born` and `f_solv` native generic regression derivative ownership for Figiel residuals.
- Any existing `Backend_unavailable` checks tied to these paths should be reworked to CppAD/analytic/implicit derivatives or explicit unsupported-path diagnostics.
- Tests in this tranche should verify absence of backend-unavailable assertions as required behavior for solved native slices.

## Post-audit blocker list (active scan 2026-05-11)

`rg -l "finite[_-]difference|Backend_unavailable|Backend unavailable" --glob '!docs/goals/**' src tests docs/pages`
currently returns remaining execution-surface hits in:

- `src/epcsaft/bindings.cpp`
- `src/epcsaft/epcsaft.py`
- `src/epcsaft/equilibrium.py`
- `src/epcsaft/ipopt_backend.py`
- `src/epcsaft/native/epcsaft_ares.cpp`
- `src/epcsaft/native/epcsaft_activity.cpp`
- `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp`
- `src/epcsaft/native/epcsaft_chemical_equilibrium.h`
- `src/epcsaft/native/epcsaft_equilibrium.cpp`
- `src/epcsaft/native/epcsaft_electrolyte.h`
- `src/epcsaft/native/epcsaft_fugcoef.cpp`
- `src/epcsaft/native/epcsaft_parameter_setup.cpp`
- `src/epcsaft/native/epcsaft_regression.cpp`
- `src/epcsaft/native/epcsaft_Z.cpp`
- `src/epcsaft/native/regression/regression_types.cpp`
- `src/epcsaft/native/regression/regression_types.h`
- `src/epcsaft/native/regression/thermo_regression.cpp`
- `src/epcsaft/native/autodiff/ad_derivative_checks.h`
- `src/epcsaft/native/autodiff/debug_gate.h`
- `src/epcsaft/native_regression.py`
- `src/epcsaft/parameters.py`
- `src/epcsaft/reactive_regression.py`
- `src/epcsaft/reactive_speciation.py`
- `src/epcsaft/benchmarks/native_ceres_thermo_regression.py`
- `src/epcsaft/benchmarks/native_regression.py`
- `src/epcsaft/benchmarks/neutral_equilibrium.py`
- `src/epcsaft/benchmarks/reactive_regression.py`
- `src/epcsaft/regression.py`
- `src/epcsaft/runtime.py`
- `tests/api/test_Backend_unavailable_debug_policy.py`
- `tests/api/test_reactive_regression.py`
- `tests/api/test_regression_api.py`
- `tests/api/test_reactive_speciation.py`
- `tests/api/test_runtime.py`
- `tests/equilibrium/test_lle.py`
- `tests/native/test_cppad_bubble_derivatives.py`
- `tests/native/test_cppad_eos_derivatives.py`
- `tests/native/test_cppad_reactive_speciation_derivatives.py`
- `tests/native/test_chemical_equilibrium_native.py`
- `tests/native/test_equation_registry.py`
- `tests/native/test_equilibrium_native_contracts.py`
- `tests/native/test_native_reactive_regression.py`
- `tests/native/test_native_regression_autodiff.py`
- `tests/native/test_native_regression_types.py`
- `tests/native/test_runtime_contracts.py`
- `tests/workflows/test_benchmark_native_regression.py`
- `docs/pages/diagnostics.rst`
- `docs/pages/development_workflows.rst`
- `docs/pages/electrolyte_vle_reactive_workflow.rst`
- `docs/pages/equilibrium_cookbook.rst`
- `docs/pages/parameter_regression.rst`
- `docs/pages/user_options.rst`

## Notes

- Keep this as a blocker-level constraint until the thread proves backend-unavailable-free behavior for the active Figiel-supported residual slice.

