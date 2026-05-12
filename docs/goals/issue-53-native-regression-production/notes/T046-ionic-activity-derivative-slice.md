T046 receipt: ionic component-activity derivative slice
======================================================

Status
------

T046 is complete.

This worker tranche now covers both halves of the scoped activity slice:

- backend-unavailable-free state-level derivative support for the supported
  nonassociating ionic `epcsaft_component_activity` path
- native chemical-equilibrium and native thermodynamic regression support for
  activity-coupled `reactive_speciation` `logK` rows on that same slice

What changed
------------

1. Added state-level `ln(phi)` derivative substrates:

   - fixed-`rho` composition derivatives
   - fixed-composition density derivatives

2. Added a native MIAC log-derivative payload in `epcsaft_activity.cpp`:

   - computes `d log(gamma) / d log(n)` for the supported ionic,
     nonassociating, fixed-pressure component-activity slice
   - includes both current-state and reference-state pressure-closure terms

3. Extended native chemical-equilibrium autodiff routing:

   - `mole_fraction_activity` + `epcsaft_component_activity` now routes to
     `chemical_equilibrium:mole_fraction_activity:log_amounts:component_activity_cppad`
   - ideal and concentration slices remain unchanged

4. Extended native thermodynamic regression support:

   - activity-coupled `reactive_speciation` `logK` rows are now accepted by the
     Ceres implicit derivative gate
   - mixed concentration+activity rows remain honestly unsupported

5. Added activity-coupled benchmark coverage:

   - `reactive_speciation_activity_logk_implicit`

What is verified
----------------

- native activity-coupled chemical-equilibrium autodiff Jacobian matches the
  existing backend-unavailable path on the supported salt case
- public `solve_reactive_speciation(..., jacobian_backend="auto")` routes the
  supported activity-standard-state case to autodiff
- native Ceres thermodynamic regression now accepts the activity-standard-state
  `logK` row and reports `derivative_backend="cppad_implicit"`
- activity benchmark runs the native hot loop with `backend=ceres` and
  `derivative=cppad_implicit`

Validation evidence
-------------------

- `uv run python scripts/build_epcsaft.py --build-only --parallel 10`
  - passed
- `uv run python run_pytest.py tests/native/test_chemical_equilibrium_native.py tests/api/test_reactive_speciation.py tests/native/test_runtime_contracts.py tests/native/test_cppad_reactive_speciation_derivatives.py tests/native/test_native_ceres_thermodynamic_regression.py tests/workflows/test_benchmark_native_regression.py -q`
  - `84 passed, 1 skipped in 7.16s`
- `uv run python scripts/benchmark_native_ceres_thermo_regression.py --case reactive_speciation_activity_logk_implicit --warmup 1 --repeat 1`
  - `reactive_speciation_activity_logk_implicit 1891.322 converged ceres cppad_implicit True 0.020167 6.34303e-20`

Key files
---------

- `src/epcsaft/native/epcsaft_activity.cpp`
- `src/epcsaft/native/epcsaft_fugcoef.cpp`
- `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp`
- `src/epcsaft/native/regression/thermo_regression.cpp`
- `src/epcsaft/native/epcsaft_core_internal.h`
- `src/epcsaft/native/epcsaft_electrolyte.h`
- `src/epcsaft/benchmarks/native_ceres_thermo_regression.py`
- `tests/native/test_chemical_equilibrium_native.py`
- `tests/native/test_native_ceres_thermodynamic_regression.py`
- `tests/api/test_reactive_speciation.py`

Remaining issue #53 scope after T046
------------------------------------

This does not complete issue #53.

The next honest unresolved surfaces are still:

- bubble-pressure differentiation for native thermodynamic regression
- non-`logK` thermodynamic parameter sensitivities such as `d_born`, `f_solv`,
  and generic binary parameters through the production native regression path
- broader activity/concentration mixed standard-state combinations

