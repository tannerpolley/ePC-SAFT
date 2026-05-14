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
bounded least-squares production path, its supported target families, and the
fact that Python still owns row orchestration and bounded step control while the
thermodynamic calculations remain native.

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

   uv run python scripts/benchmarks/benchmark_reactive_regression.py --warmup 3 --repeat 10 --json build/benchmarks/reactive_regression_main.json
   uv run python scripts/benchmarks/benchmark_reactive_regression.py --case reactive_regression_objective_tiny --warmup 3 --repeat 20 --json build/benchmarks/reactive_regression_objective_main.json
   uv run python scripts/benchmarks/benchmark_reactive_regression.py --case reactive_regression_parameter_perturbation --warmup 3 --repeat 20 --json build/benchmarks/reactive_regression_perturbation_main.json
   uv run python scripts/benchmarks/benchmark_reactive_regression.py --case reactive_regression_pressure_speciation_35_row_surrogate --warmup 0 --repeat 1 --json build/benchmarks/reactive_regression_pressure_speciation_35row_smoke.json

Benchmark JSON excludes failed repeats from timing statistics, records the
number of measured successful repeats, and carries failure messages separately.
The 35-row pressure/speciation surrogate is an opt-in smoke for mixed residual
coverage; it reports ``target_family_counts`` so CI or release notes can prove
that pressure and speciation residual families both ran. It is intentionally
excluded from the default all-case benchmark command because it is a slower
end-to-end mixed residual check.
