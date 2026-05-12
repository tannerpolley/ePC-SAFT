Diagnostics
===========

Runtime Capabilities
--------------------

``epcsaft.capabilities()`` reports installed runtime features, solver backends,
and package workflow status. The important backend labels are:

* ``native``: production native runtime path.
* ``native_state_fugacity_with_python_scalar_root``: Python owns scalar
  iteration while native states and fugacity calculations own thermodynamics.
* ``python_batched_native_solvers``: Python batches rows and manages seeds;
  native solvers evaluate thermodynamics and residuals.
* ``legacy wrapper``: compatibility path retained for older callers.
* ``experimental IPOPT``: explicit opt-in residual-minimization refinement.

For reactive electrolyte regression, inspect
``capabilities()["regression"]["reactive_electrolyte_batch_context"]["fit_status_contract"]``.
It lists the public fit statuses, the top-level convergence fields, and the
empty ``public_placeholder_statuses`` list that downstream agents can use to
confirm the package is not returning provisional labels. The sibling
``bounded_mixed_pressure_speciation_regression`` capability advertises the
supported target families and the current production-readiness gate. The public
default fit backend now calls the native residual-record solve boundary, while
``backend="python_compat"`` is the explicit compatibility-only Python
Gauss-Newton path. The capability reports
``status="partial_native_ceres_thermodynamic_slice"`` and
``issue53_native_production_ready=False`` until full Ceres parameter iteration
and all production parameter sensitivities are wired. The first supported Ceres
thermodynamic slice is native reactive speciation with reaction logK
parameters, ideal mole-fraction standard states, and speciation targets.

For reactive speciation, ``capabilities()["equilibrium"]["reactive_speciation"]``
now reports the implemented default Jacobian truth:

* ``jacobian_auto_policy="cppad_supported_else_debug_fd_or_backend_unavailable"``
* auto/default supports ``ideal_mole_fraction``, ``concentration``, and
  ``mole_fraction_activity`` on the supported native CppAD slice
* without CppAD, the ideal mole-fraction path still falls back to the exact
  analytic Jacobian
* public runtime ``state.dadx()`` auto/default remains analytic; that runtime
  policy is separate from the reactive-speciation Jacobian selection

Contribution Maps
-----------------

State objects expose contribution-map helpers:

* ``state.helmholtz_contributions()``
* ``state.residual_helmholtz_contributions()``
* ``state.pressure_contributions()``
* ``state.chemical_potential_contributions()``
* ``state.ln_fugacity_coefficient_contributions()``

Public contribution family names are ``hard_chain``, ``dispersion``,
``association``, ``ionic``, and ``born``. Inactive terms are retained with zero
values when the native runtime evaluates them as inactive.

Activity coefficients are available through ``state.activity_coefficient(...)``.
Additive activity-coefficient term decomposition is not currently exposed by
the native activity API, so ``state.activity_coefficient_contributions()``
raises ``NotImplementedError`` instead of returning invented terms.

Reactive Regression Benchmarks
------------------------------

Use the benchmark script when changing the reactive batch/context layer:

.. code-block:: powershell

   uv run python scripts/benchmark_reactive_regression.py --warmup 3 --repeat 10 --json build/benchmarks/reactive_regression_main.json
   uv run python scripts/benchmark_reactive_regression.py --case reactive_regression_objective_tiny --warmup 3 --repeat 20 --json build/benchmarks/reactive_regression_objective_main.json
   uv run python scripts/benchmark_reactive_regression.py --case reactive_regression_parameter_perturbation --warmup 3 --repeat 20 --json build/benchmarks/reactive_regression_perturbation_main.json
   uv run python scripts/benchmark_reactive_regression.py --case reactive_regression_pressure_speciation_35_row_surrogate --warmup 0 --repeat 1 --json build/benchmarks/reactive_regression_pressure_speciation_35row_smoke.json

Benchmark JSON excludes failed repeats from timing statistics, records the
number of measured successful repeats, and carries failure messages separately.
The 35-row pressure/speciation surrogate is an opt-in smoke for mixed residual
coverage; it reports ``target_family_counts`` so CI or release notes can prove
that pressure and speciation residual families both ran. It is intentionally
excluded from the default all-case benchmark command because it is a slower
end-to-end mixed residual check.

Native Regression Benchmarks
----------------------------

Use the native benchmark script when changing the fixed-shape native regression
contract, status handling, derivative policy, or public production boundary:

.. code-block:: powershell

   uv run python scripts/benchmark_native_regression.py --warmup 1 --repeat 3 --json build/benchmarks/native_regression_main.json
   uv run python scripts/benchmark_native_regression.py --case native_mea_pressure_speciation_35_row_surrogate --warmup 0 --repeat 1 --json build/benchmarks/native_regression_mea_35row_smoke.json
   uv run python scripts/benchmark_native_ceres_thermo_regression.py --warmup 1 --repeat 3 --json build/benchmarks/native_ceres_thermo_regression.json

The native benchmark includes tiny neutral, generic binary ``k_ij``, and
reactive Born-SSM+DS ``d_born``/``f_solv`` fixtures plus a 35-row public
MEA-style pressure/speciation surrogate. The benchmark payload records target
families, parameter kinds, fixed-shape residual status, derivative backend, and
whether finite differences are allowed in production.

The Ceres thermodynamic benchmark reports ``backend_unavailable`` on builds
without ``EPCSAFT_ENABLE_CERES=ON``. On a Ceres-enabled build it should report
``optimizer_backend="ceres"``, ``native_hot_loop=True``,
``python_objective_used=False``, ``finite_difference_used=False``, and a lower
``final_cost`` than ``initial_cost`` for the supported slice.

The runtime capability entry
``capabilities()["regression"]["reactive_electrolyte_batch_context"]["native_ceres_thermodynamic_regression"]``
tracks the same supported slice. It now advertises:

* ``derivative_backend="cppad_implicit"`` on CppAD-enabled builds
* reactive-speciation ``logK`` rows for ``ideal_mole_fraction``,
  ``concentration``, and supported ``mole_fraction_activity`` standard states
* the supported Born-SSM+DS parameter lane for ``born_radius`` / ``d_born`` and
  ``f_solv`` / ``solvation_factor``

Generic binary interaction parameters such as ``k_ij`` remain blocked on this
thermodynamic regression path and should still use the direct binary-regression
workflow.
