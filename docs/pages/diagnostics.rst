Diagnostics
===========

Runtime Capabilities
--------------------

``epcsaft.capabilities()`` reports installed runtime features, solver backends,
and implemented package workflows. The important backend labels are:

* ``native``: production native runtime path.
* ``native_ipopt_equilibrium_nlp``: implemented production equilibrium route
  owned by the native Ipopt NLP layer.
* ``batch_residual_evaluation_context``: Python batches rows and formats
  diagnostics for a residual context; this is not a production optimizer.

The static evidence behind route names, derivative coverage rows, pytest
slices, and ``scripts/dev/validate_project.py`` modes lives in
``epcsaft.capability_evidence``. ``capabilities()`` adapts that registry to the
currently installed native dependencies instead of duplicating route lists in
CLI wrappers. This keeps public command names stable while still tying
capability claims to executable checks.

Use ``uv run python run_pytest.py --native-contracts -q`` for native route
metadata and result-diagnostics checks. That slice is intentionally separate
from the broad native route-builder solver suite, which is slow and should be
targeted only by a single test node or with explicit long-test opt-in.

For reactive electrolyte regression, the
``mixed_pressure_speciation_residual_context`` capability advertises the
diagnostic residual context, its supported target families, and the fact that it
is not a production optimizer. Reactive electrolyte parameter fitting is not
public until native Ceres owns the optimizer route with exact derivatives.

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
* ``derivative_backend`` reports ``analytic``.
* ``ipopt_solver_ran`` and ``ipopt_accepted`` describe the native NLP solve.

Reactive electrolyte bubble result fields are the staged native diagnostics
shape. When the homogeneous speciation stage and the native Ipopt fixed-liquid
electrolyte bubble route are both available, results contain nested
dictionaries:

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
raises ``InputError`` instead of returning invented terms.

Reactive Regression Benchmarks
------------------------------

Use the benchmark script when changing the reactive batch/context layer:

.. code-block:: powershell

   uv run python scripts/benchmarks/benchmark_reactive_regression.py --warmup 3 --repeat 10 --json build/benchmarks/reactive_regression_main.json
   uv run python scripts/benchmarks/benchmark_reactive_regression.py --case reactive_regression_objective_tiny --warmup 3 --repeat 20 --json build/benchmarks/reactive_regression_objective_main.json
   uv run python scripts/benchmarks/benchmark_reactive_regression.py --case reactive_regression_pressure_speciation_35_row_surrogate --warmup 0 --repeat 1 --json build/benchmarks/reactive_regression_pressure_speciation_35row_smoke.json

Benchmark JSON excludes failed repeats from timing statistics, records the
number of measured successful repeats, and carries failure messages separately.
The 35-row pressure/speciation surrogate is an opt-in smoke for mixed residual
coverage; it reports ``target_family_counts`` so CI or release notes can prove
that pressure and speciation residual families both ran. It is intentionally
excluded from the default all-case benchmark command because it is a slower
end-to-end mixed residual check.
