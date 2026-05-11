Parameter Regression
====================

The package-owned regression helpers are record-driven. They accept flat in-memory records, tabular objects, or CSV files and return ``FitResult`` objects. Nothing is written to a parameter folder until ``write_fit_result(...)`` is called.

Supported workflows
-------------------

- ``fit_pure_neutral(...)`` fits nonassociating neutral pure-component ``m``, ``s``, and ``e`` against density and vapor-pressure records with the native least-squares backend.
- ``fit_pure_ion(...)`` fits ion ``s`` and ``e`` by default, and can fit ``d_born`` when requested, with the native least-squares backend and provenance guardrails.
- ``fit_binary_pair(...)`` fits constant binary interaction values from direct VLE x/y records with the native least-squares backend and provenance guardrails.

Ion and binary V1 intentionally do not add dataset manifests or new regression-specific parameter namespaces. The helpers build runtime states from the existing dataset loader and caller-provided records.

Non-native optimizer loops are not an approved production backend for package-owned regression helpers. Python code may prepare records, declare provenance, and call native regression, but coupled electrolyte, reactive, phase-equilibrium, ``d_born``, and ``k_ij`` fitting should use the native backend.

Build prerequisite
------------------

There is no IPOPT prerequisite for the default package build. The supported
developer path is the uv-managed environment plus the direct CMake/pybind11
native build:

.. code-block:: powershell

   uv sync --no-install-project
   uv run python scripts\build_epcsaft.py

For experimental IPOPT work, install the optional Python adapter dependency
with the ``ipopt`` extra or dependency group:

.. code-block:: powershell

   uv sync --extra ipopt

On Windows, ``cyipopt`` builds from source and needs an IPOPT install that
``pkg-config`` can describe. When using conda-forge IPOPT, first make sure the
IPOPT environment also has ``pkg-config``, then let the helper create a corrected
local ``ipopt.pc`` shim and run uv:

.. code-block:: powershell

   conda install -n epcsaft-cyipopt-test -c conda-forge pkg-config
   .\scripts\setup_windows_cyipopt_uv.ps1 -IpoptPrefix C:\ProgramData\miniconda3\envs\epcsaft-cyipopt-test\Library

The IPOPT path uses ``cyipopt`` and remains explicit opt-in through
``solver_backend="ipopt"``. It is not a replacement for native least-squares or
Newton defaults.

Current IPOPT scope
-------------------

The current IPOPT support is an experimental
``bound_constrained_residual_minimization`` backend for selected equilibrium
routes. It does not yet expose a full constrained thermodynamic NLP: material
balances, charge balance, and equilibrium residuals are not formal IPOPT
equality constraints in this phase. Runtime capabilities report
``full_constrained_nlp_available=False`` and ``default_auto_uses_ipopt=False``.

Approximate Hessian diagnostics are intentionally explicit. Gauss-Newton means a
least-squares ``J.T @ J`` callback; L-BFGS means IPOPT limited-memory Hessian
approximation. Both report ``exact_hessian_available=False`` and
``hessian_includes_second_residual_derivatives=False``.

Solver-selection guidance
-------------------------

.. list-table::
   :header-rows: 1

   * - Problem type
     - Preferred default
     - IPOPT role
   * - Scalar bubble/dew variable solve
     - Brent or safeguarded Newton
     - Usually not appropriate
   * - Small smooth residual system
     - Native Newton/trust-region residual solve
     - Explicit fallback or refinement only
   * - Least-squares parameter estimation
     - Gauss-Newton, LM, or trust-region least squares
     - Use only when bounds or constraints dominate
   * - Noisy or nonsmooth black-box workflow
     - Derivative-free or surrogate method
     - Not preferred
   * - Phase equilibrium near active bounds
     - Native safeguarded residual solve first
     - Useful as explicit bounded residual refinement
   * - Large sparse constrained NLP
     - IPOPT or SQP with sparse derivatives
     - Appropriate once constraints are explicit

Create a dataset folder
-----------------------

Start from a user-owned template folder.

.. code-block:: python

   from epcsaft import create_parameter_template

   template_root = create_parameter_template(
       location=r"C:\Users\Tanner\Documents\my_epcsaft_data",
       folder_name="regression_case",
       species=["H2O", "Na+", "Cl-"],
   )

The template gives you user-owned ``pure/`` and ``mixed/binary_interaction/`` CSV files that ``write_fit_result(...)`` can update after a fit.

Neutral records
---------------

``fit_pure_neutral(...)`` accepts flat records with:

- ``T``: temperature in kelvin
- ``P``: vapor pressure in pascals
- ``rho``: liquid molar density in ``mol/m^3``
- or ``rho_kg_m3`` / ``rho_sat_liq_kg_m3``: liquid mass density in ``kg/m^3``
- optional ``phase``: defaults to ``liq``

Example:

.. code-block:: python

   from epcsaft import fit_pure_neutral

   records = [
       {"T": 100.0, "P": 34375.892, "rho_sat_liq_kg_m3": 438.88524},
       {"T": 110.0, "P": 88130.038, "rho_sat_liq_kg_m3": 424.77725},
   ]

   result = fit_pure_neutral(
       records,
       "Methane",
       assoc_scheme="",
       fixed_parameters={
           "MW": 0.0160428,
           "z": 0.0,
           "e_assoc": 0.0,
           "vol_a": 0.0,
           "dielc": 8.0,
           "d_born": 0.0,
           "f_solv": 1.0,
       },
       initial_guess={"m": 1.1, "s": 3.6, "e": 145.0},
   )

Ion records
-----------

``fit_pure_ion(...)`` records require ``T`` and ``P`` plus one composition basis:

- full mole-fraction columns such as ``x_H2O``, ``x_Na+``, and ``x_Cl-``
- or ``molality`` with explicit ``species=[...]`` and ``solvent=...``

Each ion regression problem must include at least one of:

- ``osmotic_coefficient`` or ``osmotic``
- ``mean_ionic_activity``, ``mean_ionic_activity_coefficient``, or ``miac``

Density is optional and is included when ``rho`` or a supported mass-density column is present.
``d_born`` fitting additionally requires electrostatic provenance: dielectric or relative-permittivity data, ion-activity/osmotic data, or an explicit override. Results include ``result.provenance_report`` so downstream workflows can distinguish supported fitted values from provisional diagnostic values.

Example:

.. code-block:: python

   from epcsaft import fit_pure_ion

   records = [
       {"T": 298.15, "P": 101325.0, "molality": 0.1, "osmotic_coefficient": 0.933},
       {"T": 298.15, "P": 101325.0, "molality": 0.2, "mean_ionic_activity": 0.735},
   ]

   result = fit_pure_ion(
       records,
       "Na+",
       dataset="2026_Khudaida",
       species=["H2O", "Na+", "Cl-"],
       solvent="H2O",
       initial_guess={"s": 2.6, "e": 210.0},
       bounds={"s": (2.4, 3.2), "e": (150.0, 300.0)},
   )

Advanced electrolyte options pass through to every runtime state build:

.. code-block:: python

   user_options = {
       "elec_model": {
           "rel_perm": {"rule": "empirical", "differential_mode": "auto"},
           "born_model": {
               "d_Born_mode": 3,
               "solvation_shell_model": True,
               "dielectric_saturation": True,
               "mu_born_model": {"differential_mode": "auto", "comp_dep_delta_d": True},
           },
       }
   }

   result = fit_pure_ion(
       records,
       "Na+",
       dataset="2026_Khudaida",
       species=["H2O", "Na+", "Cl-"],
       solvent="H2O",
       fit_targets=("d_born",),
       initial_guess={"d_born": 3.2},
       user_options=user_options,
   )

Provenance guardrails
---------------------

Use explicit declarations when a fit target needs to document or override its
data basis:

.. code-block:: python

   from epcsaft import BinaryInteraction, FitParameter, validate_regression_provenance

   report = validate_regression_provenance(
       [
           FitParameter("MEAH+", "d_born", source="dielectric_or_ion_activity"),
           BinaryInteraction(("MEA", "H2O"), parameter="k_ij", source="direct_binary_vle"),
       ],
       species=["MEA", "H2O", "MEAH+"],
       charges=[0.0, 0.0, 1.0],
   )

The validator rejects unsupported targets by default:

- same-sign ionic ``k_ij`` targets, because the runtime suppresses same-sign short-range dispersion;
- opposite-sign ion-pair interaction targets without direct electrolyte activity, osmotic, salt-pair, or explicit-override provenance;
- neutral-ion interaction targets without direct neutral-ion/electrolyte provenance;
- ``d_born`` targets backed only by mixed reactive VLE residuals.

``RelativePermittivityResidual`` provides a first-class record/term descriptor for dielectric data:

.. code-block:: python

   from epcsaft import RelativePermittivityResidual

   dielectric_term = RelativePermittivityResidual(
       T=298.15,
       P=101325.0,
       composition={"H2O": 0.8, "MEA": 0.2},
       epsilon_r_exp=65.0,
   ).to_fit_term(species=["H2O", "MEA"])

Reactive electrolyte regression batches
---------------------------------------

Publication-style reactive electrolyte studies can now use the package-owned
batch/context layer instead of packing row loops downstream by hand.

Main entry points:

- ``ReactiveElectrolyteRow``: one reactive speciation or fixed-liquid bubble row
- ``ReactiveElectrolyteBatch``: shared species, balances, reactions, parameter payload, and solver options
- ``ReactiveElectrolyteRegressionContext.from_batch(...)``: compile invariant row/schema metadata once
- ``evaluate_reactive_regression_objective(...)``: evaluate a structured mixed residual objective
- ``fit_reactive_electrolyte_parameters(...)``: production native thermodynamic regression for the supported
  reactive-speciation slice, with an explicit residual-record compatibility boundary for broader legacy cases
- ``summarize_regression_result(...)`` plus ``write_regression_*`` helpers: stable JSON/CSV reporting

Minimal example:

.. code-block:: python

   import epcsaft
   import numpy as np

   params = {
       "m": np.asarray([1.2047, 1.0, 1.0]),
       "s": np.asarray([2.79, 2.82, 2.76]),
       "e": np.asarray([353.95, 230.0, 170.0]),
       "MW": np.asarray([18.01528e-3, 22.989e-3, 35.45e-3]),
       "z": np.asarray([0.0, 1.0, -1.0]),
       "dielc": np.asarray([78.09, 8.0, 8.0]),
       "d_born": np.asarray([0.0, 3.445, 4.1]),
   }

   batch = epcsaft.ReactiveElectrolyteBatch(
       species=["H2O", "Na+", "Cl-"],
       rows=[
           epcsaft.ReactiveElectrolyteRow(
               row_id="row_1",
               T=298.15,
               P_seed=101325.0,
               totals={"water": 0.98, "sodium": 0.01, "chloride": 0.01},
               initial_x=[0.98, 0.01, 0.01],
               balances={
                   "water": {"H2O": 1.0},
                   "sodium": {"Na+": 1.0},
                   "chloride": {"Cl-": 1.0},
               },
               reactions=[],
               vapor_species=["H2O"],
               target_partial_pressures={"H2O": 3166.4},
               target_speciation={"H2O": 0.98},
               source="train",
               split="fit",
           ),
       ],
       balances={
           "water": {"H2O": 1.0},
           "sodium": {"Na+": 1.0},
           "chloride": {"Cl-": 1.0},
       },
       reactions=[],
       vapor_species=["H2O"],
       base_parameters=params,
   )

   context = epcsaft.ReactiveElectrolyteRegressionContext.from_batch(
       species=batch.species,
       rows=batch.rows,
       balances=batch.balances,
       reactions=batch.reactions,
       vapor_species=batch.vapor_species,
       base_parameters=batch.base_parameters,
   )
   result = context.evaluate_objective({"Na+.sigma": 2.85})

Structured outputs follow stable field names:

- per-row: ``row_id``, ``success``, ``message``, ``composition``, ``pressure``,
  ``ln_fugacity``, ``activity_coefficients``, ``density``,
  ``relative_permittivity``, ``residuals``, ``residual_names``,
  ``failure_diagnostics``, ``active_bounds``, ``solver_status``,
  ``elapsed_seconds``, ``cache_stats``, ``warm_start_used``,
  ``warm_start_source``, ``warm_start_failed``, ``fallback_seed_used``,
  ``partial_pressures``, ``y_vap``, ``named_reaction_residuals``, ``source``,
  ``split``, ``metadata``
- batch: ``success_count``, ``failure_count``, ``row_results``, ``residuals``,
  ``residual_names``, ``residual_row_map``, ``diagnostics``, ``cache_stats``,
  ``timing_summary``
- objective: ``objective``, ``metrics``, and the embedded ``batch_result``
- fit: ``success``, ``message``, ``status``, ``termination_reason``,
  ``iterations``, ``objective_initial``, ``objective_final``,
  ``gradient_norm``, ``step_norm``, ``parameter_map``, ``seed_map``,
  ``lower_bounds``, ``upper_bounds``, ``active_bounds``, ``objective_result``,
  ``covariance_available``, ``covariance_matrix``, ``identifiability_status``,
  and ``diagnostics``

Reporting helpers write those schemas without downstream column guessing:

.. code-block:: python

   epcsaft.write_regression_summary(result, "build/regression/summary.json")
   epcsaft.write_regression_row_table(result, "build/regression/rows.csv")
   epcsaft.write_regression_residual_table(result, "build/regression/residuals.csv")

If you are fitting parameters instead of only evaluating an objective, pass the
fit result directly to the same reporting helpers. ``write_regression_parameter_table(...)``
accepts either a raw parameter map plus ``seed_map=...`` or the full
``ReactiveRegressionFitResult`` so seed values, parameter movement, bounds, and
active-bound flags stay aligned with the solved fit payload.

Native regression contract
~~~~~~~~~~~~~~~~~~~~~~~~~~

``epcsaft.native_regression_contract_schema()`` exposes the native regression
result contract that Ceres-backed production fits must satisfy. The contract is
fixed-shape: residual names, row ids, target families, and row diagnostic fields
must remain stable even when a recoverable row failure is represented by penalty
residuals. The canonical top-level statuses are ``converged``,
``max_iterations``, ``line_search_failed``, ``singular_jacobian``,
``all_rows_failed``, ``nonfinite_objective``, ``bounds_inconsistent``,
``invalid_input``, and ``backend_unavailable``. Unsupported native optimizer or
derivative combinations use ``backend_unavailable`` rather than silently falling
back. Production native regression does not allow finite-difference derivatives;
finite-difference comparisons belong behind explicit debug gates such as
``EPCSAFT_ALLOW_FINITE_DIFFERENCE_DEBUG=1``.

The current native Ceres thermodynamic fit slice is intentionally narrow:
``reactive_speciation`` rows, ideal-mole-fraction reaction standard states,
speciation targets, and reaction ``logK`` parameters. Born-SSM+DS ``d_born`` and
``f_solv`` parameters are applied to native mixtures, but they do not yet have
production Ceres sensitivities for the activity/fugacity path. Those parameters
therefore report ``backend_unavailable`` instead of falling back to finite
differences. The missing scalar-templated derivative path is:
``NativeThermoCeresCostFunction::Evaluate`` ->
``evaluate_native_thermo_regression_rows`` ->
``chemical_equilibrium_native`` /
``evaluate_chemical_equilibrium_residual_native`` ->
``activity_coefficients`` ->
``ePCSAFTStateNative::activity_coefficient_native`` ->
``residual_chemical_potential_result_cpp`` ->
``composition_contribution_result_cpp`` ->
``ares_contributions_cpp`` ->
``born_intermediate_state_cpp`` / ``dadx_born_cpp`` for Born-SSM+DS ``d_born``
and ``f_solv``, including any ``solve_density_scoped`` pressure-closure
sensitivities used by activity or concentration standard states.

Binary interaction parameters are a separate regression family. Use
``fit_binary_pair(...)`` for ``k_ij``, ``l_ij``, and ``k_hb_ij`` against direct
binary VLE composition data; do not treat those generic binary fits as reactive
speciation/Born-SSM+DS regression tests.

Reactive-electrolyte bubble-pressure rows are also native-evaluated but not yet
production-differentiated for Ceres. The missing residual system is the coupled
log-pressure and vapor-composition solve used by
``electrolyte_bubble_pressure_native``:

- liquid state fugacity/activity at fixed liquid composition;
- vapor submixture state fugacity at the trial pressure and vapor composition;
- vapor normalization;
- fugacity equality residuals
  ``log(y_i) + log(phi_i^vap) - log(x_i) - log(phi_i^liq)``;
- pressure target residuals through the solved ``P``.

Until those residuals expose analytic, CppAD, or implicit sensitivities with
respect to fitted parameters and continuation variables, the Ceres fit path
returns ``backend_unavailable`` for bubble-pressure rows.

The fit helper defaults to ``backend="native"``. For the currently supported
production slice (reactive speciation rows, ideal-mole-fraction reaction
standard states, linear speciation targets, and reaction ``logK`` parameters),
Python serializes the batch once and calls ``fit_native_thermo_regression(...)``;
the C++ layer owns the thermodynamic row solves, Ceres loop, and derivative
policy. Broader native compatibility cases can still use
``backend="native_residual_records"``, which evaluates the Python objective once
and sends fixed-shape residual records across the native boundary. The legacy
Python Gauss-Newton loop remains available only as ``backend="python_compat"``
for old comparison workflows and is labeled ``production_ready = false`` in
diagnostics.

The native fit helper accepts the compiled batch or context plus the initial
parameter map and optional fit controls:

.. code-block:: python

   fit = epcsaft.fit_reactive_electrolyte_parameters(
       context,
       initial_parameters={"Na+.sigma": 2.85},
       lower_bounds={"Na+.sigma": 2.5},
       upper_bounds={"Na+.sigma": 3.1},
       max_iterations=6,
       tolerance=1e-6,
       derivative_backend="analytic",
   )

For compatibility-only comparisons of the old Python optimizer:

.. code-block:: python

   compat_fit = epcsaft.fit_reactive_electrolyte_parameters(
       context,
       initial_parameters={"Na+.sigma": 2.85},
       lower_bounds={"Na+.sigma": 2.5},
       upper_bounds={"Na+.sigma": 3.1},
       backend="python_compat",
       jacobian_mode="central",
   )

``summarize_regression_result(...)`` returns ``fit_success = null`` for
objective-only results and a boolean for true fit results. Fit results also
carry a production status contract:

.. list-table::
   :header-rows: 1

   * - ``status``
     - Meaning
   * - ``converged``
     - The bounded fit met a stopping tolerance and all objective rows solved.
   * - ``max_iterations``
     - The bounded fit improved or evaluated normally but exhausted the
       configured iteration budget.
   * - ``line_search_failed``
     - The Gauss-Newton step and bounded line search could not improve the
       current objective.
   * - ``failed_rows``
     - The final objective includes failed rows; inspect row diagnostics before
       using the parameters.

The package does not emit placeholder public statuses such as
``bounded_incomplete``. Downstream workflows should branch on ``status``,
``termination_reason``, ``objective_initial``, ``objective_final``,
``gradient_norm``, and ``step_norm`` instead of inventing their own fit-state
labels.

The package-owned micro-benchmark harness for this layer is:

.. code-block:: powershell

   uv run python scripts\benchmark_reactive_regression.py --warmup 3 --repeat 10
   uv run python scripts\benchmark_reactive_regression.py --case reactive_regression_pressure_speciation_35_row_surrogate --warmup 0 --repeat 1
   uv run python scripts\benchmark_native_regression.py --warmup 1 --repeat 3
   uv run python scripts\benchmark_native_regression.py --case native_mea_pressure_speciation_35_row_surrogate --warmup 0 --repeat 1
   uv run python scripts\benchmark_native_ceres_thermo_regression.py --warmup 1 --repeat 3

``benchmark_native_ceres_thermo_regression.py`` reports
``native_hot_loop``, ``python_objective_used``, ``finite_difference_used``,
``initial_cost``, and ``final_cost`` for the supported Ceres thermodynamic
slice. With a default build that does not enable Ceres, the benchmark exits
successfully but reports ``backend_unavailable``. Rebuild with
``EPCSAFT_ENABLE_CERES=ON`` to verify Ceres-owned parameter iteration and
objective decrease.

Binary VLE records
------------------

``fit_binary_pair(...)`` V1 supports VLE x/y records only. Records require:

- ``T``: temperature in kelvin
- ``P``: pressure in pascals
- liquid mole-fraction columns such as ``x_H2O`` and ``x_Ethanol``
- vapor mole-fraction columns such as ``y_H2O`` and ``y_Ethanol``

The V1 constant targets are ``k_ij``, ``l_ij``, and ``k_hb_ij``. Linear temperature models and LLE fitting are future phases and raise ``InputError``. Ion-involving binary targets require explicit provenance and are rejected by default unless they are tied to direct electrolyte/neutral-ion data or an explicit override.

Example:

.. code-block:: python

   from epcsaft import fit_binary_pair

   records = [
       {"T": 330.0, "P": 101325.0, "x_H2O": 0.7, "x_Ethanol": 0.3, "y_H2O": 0.5, "y_Ethanol": 0.5},
       {"T": 340.0, "P": 101325.0, "x_H2O": 0.6, "x_Ethanol": 0.4, "y_H2O": 0.4, "y_Ethanol": 0.6},
   ]

   result = fit_binary_pair(
       records,
       ("H2O", "Ethanol"),
       dataset="2026_Khudaida",
       initial_guess={"k_ij": -0.02},
       bounds={"k_ij": (-0.2, 0.2)},
   )

Inspect and write results
-------------------------

.. code-block:: python

   print(result.success)
   print(result.backend)
   print(result.jacobian_backend)
   print(result.hessian_backend)
   print(result.fitted_values)
   print(result.metrics_by_term)
   print(result.provenance_report)

   from epcsaft import write_fit_result

   written_paths = write_fit_result(result, template_root, overwrite=False)
   print(written_paths)

With ``overwrite=False``, blank template cells can be filled but existing values are protected. Pure ion fits update the target component row in ``pure/``. Binary fits update both symmetric cells in the relevant interaction matrix.

Derivative and Jacobian access
------------------------------

Normal ``FitResult`` payloads report compact derivative metadata:

- ``jacobian_available``
- ``jacobian_backend``
- ``jacobian_fallback_used``
- ``jacobian_fallback_reason``
- ``finite_difference_fallback_count``
- ``hessian_available``
- ``hessian_backend``
- ``hessian_fallback_used``
- ``hessian_fallback_reason``

Large matrices are exposed only through explicit derivative-evaluation helpers. Use ``evaluate_pure_neutral_derivatives(...)`` for the native pure-neutral objective. It returns residuals, gradient, ``jacobian_row_major``, ``jacobian_shape``, and Hessian skeleton fields. Pure-neutral Jacobians use the native autodiff path.

For lower-level generic native records, use ``evaluate_generic_regression_derivatives(..., jacobian_backend="finite_difference")`` when an explicit finite-difference comparison Jacobian is needed. The default ``auto`` mode raises until generic residual state calls have analytic/autodiff coverage.

Reactive electrolyte residual evaluator
---------------------------------------

Use ``evaluate_reactive_electrolyte_bubble_residuals(...)`` when a downstream
project needs a fixed-shape objective for coupled native reactive speciation plus
fixed-liquid electrolyte bubble pressure. The helper is deliberately not an
optimizer and it does not own MEA-specific data, parameter masks, run folders, or
artifact promotion. Downstream code supplies records, targets, species,
balances, reactions, and a ``mixture_factory`` for the current parameters; the
package returns a ``ReactiveElectrolyteRegressionResult`` with residuals,
residual names, per-record diagnostics, and success/failure counts.

The evaluator keeps the residual vector shape stable when a record fails by
inserting bounded penalty residuals. It also forces result-mode reactive bubble
solves internally, so one bad row can be reported without aborting the whole
candidate. Successful row diagnostics include predicted partial pressures,
liquid composition, vapor composition, named reaction residuals, and compact
solver diagnostics so downstream code can write reports without rerunning the
same expensive rows. Keep target magnitudes positive for log-scale pressure and
composition residuals, and use continuation when neighboring rows are ordered by
temperature, loading, or pressure:

.. code-block:: python

   result = epcsaft.evaluate_reactive_electrolyte_bubble_residuals(
       records,
       species=["CO2", "H2O", "MEA", "MEAH+", "HCO3-"],
       mixture_factory=make_mixture_for_candidate,
       balances=balances,
       reactions=reactions,
       vapor_species=["CO2", "H2O"],
       pressure_species=["CO2"],
       speciation_species=["CO2", "MEAH+", "HCO3-"],
       reaction_names=["carbamate", "bicarbonate"],
       continuation="auto",
   )

Pass ``result.residuals`` to a downstream-owned optimizer. Use
``result.record_results`` and ``result.diagnostics`` to decide whether a
candidate failed numerically, hit phase-equilibrium limits, or simply predicts
poorly. Do not treat this helper as a full constrained Gibbs/NLP solve; IPOPT
remains an explicit opt-in refinement elsewhere and is not used automatically.

Derivative availability
-----------------------

.. list-table::
   :header-rows: 1

   * - Method
     - Current Jacobian access
     - Hessian status
   * - Runtime ``dadt()``, ``dadx()``, ``z(return_contribution_terms=True)``, ``mures(return_contribution_terms=True)``
     - Analytical where available, autodiff where implemented; explicit finite difference only where requested
     - Not exposed
   * - Pure-neutral regression
     - Native autodiff Jacobian through ``evaluate_pure_neutral_derivatives(...)``
     - Skeleton metadata only
   * - Generic ion/binary regression
     - Explicit finite-difference Jacobian through ``evaluate_generic_regression_derivatives(..., jacobian_backend="finite_difference")`` until generic autodiff coverage is implemented
     - Skeleton metadata only
   * - Neutral LLE
     - Native Newton solve remains available; derivative callback paths require explicit finite difference until autodiff residual coverage is implemented
     - Skeleton metadata only
   * - Chemical equilibrium / reactive speciation
     - Analytic log-amount Jacobian for ideal-mole-fraction reactions under ``auto``; activity/concentration paths require explicit finite difference until derivative coverage is implemented
     - Opt-in cyipopt accepts Gauss-Newton or L-BFGS approximate Hessian strategies

The Hessian fields are deliberately a contract skeleton for future
IPOPT-compatible optimizer integration. They do not mean exact second-derivative
evaluation is implemented.
