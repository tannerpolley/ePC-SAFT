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

Benchmark JSON excludes failed repeats from timing statistics, records the
number of measured successful repeats, and carries failure messages separately.
