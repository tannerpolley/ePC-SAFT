Diagnostics
===========

Runtime Capabilities
--------------------

``epcsaft.capabilities()`` reports installed runtime features, solver backends,
and package workflow status. The important backend labels are:

* ``native``: production native runtime path.
* ``native_ipopt_equilibrium_nlp_required``: public route name is declared, but
  execution is gated until the native Ipopt NLP route builder owns the solve.
* ``batch_residual_evaluation_context``: Python batches rows and formats
  diagnostics for a residual-only context; this is not a production optimizer.

For reactive electrolyte regression, inspect
``capabilities()["regression"]["reactive_electrolyte_batch_context"]["fit_status_contract"]``.
It lists the public fit statuses, the top-level convergence fields, and the
``residual_evaluation_only`` status used while native Ceres derivative coverage is not
routed. The sibling
``mixed_pressure_speciation_residual_context`` capability advertises the
diagnostic residual context, its supported target families, and the fact that
it is not a production optimizer. Thermodynamic calculations remain native
while row batching and diagnostic formatting stay outside the optimizer path.

Reactive Speciation And Bubble Diagnostics
------------------------------------------

Native reactive speciation returns enough diagnostics to prove that the explicit
Ipopt ideal route and exact derivative path were actually used. Activity- and
concentration-coupled reaction constants remain route-gated until their EOS NLP
blocks exist. For accepted ideal routes, check:

* ``solver_language`` is ``c++``.
* ``native_entrypoint`` is ``_solve_chemical_equilibrium_native``.
* ``selected_solver_backend`` is ``native_ipopt``.
* ``problem_class`` is ``homogeneous_ideal_gibbs_speciation``.
* ``reaction_standard_states`` records the public reaction-constant convention.
* ``derivative_backend`` and ``derivative_status`` report ``analytic``.
* ``ipopt_solver_ran`` and ``ipopt_accepted`` describe the native NLP solve.

Reactive electrolyte bubble result fields are retained as the target structured
diagnostics shape, but the public route currently raises ``InputError`` until
the native Ipopt electrolyte bubble route builder owns the solve. After that
route lands, results will contain nested dictionaries:

* ``diagnostics["speciation"]`` is the homogeneous reactive speciation result.
* ``diagnostics["bubble"]`` is the native Ipopt electrolyte bubble-pressure
  result.
* ``partial_pressures`` maps volatile neutral species to pressure contributions.
* ``fugacity_residual_norm`` measures the volatile-neutral fugacity equality
  residual for the bubble solve.

Use these fields together when validating a CO2 + amine + water pressure and
speciation benchmark: reaction, charge, and material residual norms come from
the speciation result; CO2 partial pressure and vapor composition come from the
bubble result.

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
   uv run python scripts/benchmarks/benchmark_reactive_regression.py --case reactive_regression_parameter_shift --warmup 3 --repeat 20 --json build/benchmarks/reactive_regression_parameter_shift_main.json
   uv run python scripts/benchmarks/benchmark_reactive_regression.py --case reactive_regression_pressure_speciation_35_row_surrogate --warmup 0 --repeat 1 --json build/benchmarks/reactive_regression_pressure_speciation_35row_smoke.json

Benchmark JSON excludes failed repeats from timing statistics, records the
number of measured successful repeats, and carries failure messages separately.
The 35-row pressure/speciation surrogate is an opt-in smoke for mixed residual
coverage; it reports ``target_family_counts`` so CI or release notes can prove
that pressure and speciation residual families both ran. It is intentionally
excluded from the default all-case benchmark command because it is a slower
end-to-end mixed residual check.
