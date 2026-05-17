Parameter Regression
====================

The package-owned regression helpers are record-driven. They accept flat in-memory records, tabular objects, or CSV files and return ``FitResult`` objects. Nothing is written to a parameter folder until ``write_fit_result(...)`` is called.

Supported workflows
-------------------

- ``fit_pure_parameters(...)`` is the easy public wrapper for pure-component
  fits. It accepts issue-style names such as ``species``, ``data_rows``,
  ``parameters_to_fit``, ``fixed_parameters``, ``bounds``, ``weights``,
  ``loss``, ``solver_options``, and ``output_report``. The current
  implementation delegates to ``fit_pure_neutral(...)`` and records the easy
  API metadata on ``result.problem``.
- ``fit_binary_parameters(...)`` is the easy public wrapper for binary
  interaction fits. It accepts a two-component ``species`` pair and delegates
  to ``fit_binary_pair(...)`` for the supported constant binary-interaction
  targets.
- ``fit_liquid_electrolyte_parameters(...)`` fits liquid-electrolyte
  ``s``, ``e``, ``d_born``, ``f_solv``, and ``dielc`` targets through native
  Ceres. Supported row families are density, osmotic coefficient,
  mean-ionic activity coefficient, and relative permittivity. The path is
  liquid-only and reports ``cppad_implicit`` Jacobian metadata; it does not
  imply vapor electrolyte Born support.
- ``fit_pure_neutral(...)`` fits nonassociating neutral pure-component ``m``, ``s``, and ``e`` against density and vapor-pressure records with the native Ceres backend.
- ``fit_pure_ion(...)`` fits ion ``s``, ``e``, and ``d_born`` declarations with native provenance guardrails through the Ceres backend and ``cppad_implicit`` density-root sensitivities for density, osmotic-coefficient, and mean-ionic-activity rows.
- ``fit_binary_pair(...)`` fits supported constant ``k_ij`` binary interaction values from direct VLE x/y records with the native Ceres backend and provenance guardrails. Constant ``l_ij`` and ``k_hb_ij`` declarations are schema-supported, but fitting them raises until native analytic, CppAD, or implicit derivatives are implemented for those targets.

Ion and binary V1 intentionally do not add dataset manifests or new regression-specific parameter namespaces. The helpers build runtime states from the existing dataset loader and caller-provided records.

Non-native optimizer loops are not an approved production backend for package-owned regression helpers. Python code may prepare records, declare provenance, and call native regression, but coupled electrolyte, reactive, phase-equilibrium, ``d_born``, and ``k_ij`` fitting should use a native backend with explicit derivative metadata. Public production helpers route supported generic targets through native Ceres paths with analytic, CppAD, or implicit sensitivities.

The public easy APIs intentionally do not expose non-exact derivative
configuration. Unsupported derivative or optimizer paths raise with diagnostics
tied to the missing analytic, CppAD, or implicit derivative path.

Generic target-row schemas
--------------------------

``TargetRow`` and ``TargetDataset`` provide schema-only, application-neutral
containers for future regression and validation workflows. They validate row
families such as ``pure_density``, ``pure_vapor_pressure``, ``p_rho_t``,
``binary_vle``, ``binary_lle``, ``osmotic_coefficient``,
``mean_ionic_activity``, ``relative_permittivity``, ``activity``,
``fugacity``, ``speciation``, ``vle_partial_pressure``,
``lle_phase_composition``, and ``regularization`` without running optimizer
internals.

These containers deliberately use generic thermodynamic terms. Downstream
projects should keep process-specific metrics, reports, and application names
outside the package-owned public API.

Easy API examples
-----------------

Pure-component fits use the vocabulary a new user would expect while retaining
the native-backed lower-level implementation:

.. code-block:: python

   from epcsaft import fit_pure_parameters

   result = fit_pure_parameters(
       species="Methane",
       data_rows=[
           {"T": 100.0, "P": 34375.892, "rho_sat_liq_kg_m3": 438.88524},
           {"T": 110.0, "P": 88130.038, "rho_sat_liq_kg_m3": 424.77725},
       ],
       parameters_to_fit=("m", "sigma", "epsilon"),
       fixed_parameters={
           "MW": 0.0160428,
           "z": 0.0,
           "e_assoc": 0.0,
           "vol_a": 0.0,
           "dielc": 8.0,
           "d_born": 0.0,
           "f_solv": 1.0,
       },
       bounds={"m": (0.8, 1.4)},
       weights={"density": 2.0},
       loss="linear",
       solver_options={"max_nfev": 20},
   )

Binary interaction fits accept a two-species pair:

.. code-block:: python

   from epcsaft import fit_binary_parameters

   result = fit_binary_parameters(
       species=("H2O", "Ethanol"),
       parameter_set="2026_Khudaida",
       data_rows=[
           {
               "T": 330.0,
               "P": 101325.0,
               "x_H2O": 0.7,
               "x_Ethanol": 0.3,
               "y_H2O": 0.5,
               "y_Ethanol": 0.5,
           },
       ],
       parameters_to_fit=("k_ij",),
       bounds={"k_ij": (-0.2, 0.2)},
   )

Liquid-electrolyte fits use native Ceres for the supported liquid row families:

.. code-block:: python

   from epcsaft import fit_liquid_electrolyte_parameters

   result = fit_liquid_electrolyte_parameters(
       species=("H2O", "Na+", "Cl-"),
       parameter_set="2026_Khudaida",
       data_rows=[
           {
               "T": 298.15,
               "P": 101325.0,
               "x_H2O": 0.98,
               "x_Na+": 0.01,
               "x_Cl-": 0.01,
               "rho": 41.0,
               "osmotic_coefficient": 1.05,
               "mean_ionic_activity": 1.01,
               "epsilon_r_exp": 76.6,
           },
       ],
       parameters_to_fit=("d_born", "f_solv"),
       initial_guess={"d_born": 2.0, "f_solv": 0.5},
       bounds={"d_born": (1.0, 8.0), "f_solv": (0.1, 3.0)},
       weights={"osmotic_coefficient": 1.0},
   )

   assert result.backend == "ceres"
   assert result.jacobian_backend == "cppad_implicit"

Build prerequisite
------------------

There is no IPOPT prerequisite for the default package build. The supported
developer path is the uv-managed environment plus the direct CMake/pybind11
native build:

.. code-block:: powershell

   uv sync --no-install-project
   uv run python scripts\dev\build_epcsaft.py

For native Ipopt development, request a system Ipopt package through the native
build script or PEP 517 environment variables:

.. code-block:: powershell

   uv run python scripts\dev\build_epcsaft.py --profile ipopt --ipopt-dir C:\path\to\Ipopt\lib\cmake\Ipopt

Current IPOPT scope
-------------------

The current package can discover and link a native system Ipopt dependency when
explicitly requested. Production regression routes use native Ceres; production
equilibrium routes are gated to native Ipopt constrained-NLP builders. Runtime
capabilities report the discovered adapter state, public Ipopt route list, and
whether the current build has a constrained-NLP route available.

Solver-selection guidance
-------------------------

.. list-table::
   :header-rows: 1

   * - Problem type
     - Preferred default
     - IPOPT role
   * - Scalar bubble/dew variable solve
     - Native Ipopt constrained NLP route
     - Required for production equilibrium
   * - Small smooth residual system
     - Native Ipopt constrained NLP route for equilibrium; native Ceres for regression residuals
     - Required when the problem is an equilibrium route
   * - Least-squares parameter estimation
     - Native Ceres regression route
     - Not the production regression solver
   * - Noisy or nonsmooth black-box workflow
     - Not a production package route
     - Not a production package route
   * - Phase equilibrium near active bounds
     - Native Ipopt constrained NLP route
     - Required for production equilibrium
   * - Large sparse constrained NLP
     - Native Ipopt constrained NLP route with exact gradients/Jacobians
     - Required for production equilibrium

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

``fit_pure_ion(...)`` records require ``T`` and ``P`` plus one composition basis. The helper validates rows, targets, bounds, and provenance, then optimizes supported ion ``s``, ``e``, and ``d_born`` targets through native Ceres with ``cppad_implicit`` residual Jacobians:

- full mole-fraction columns such as ``x_H2O``, ``x_Na+``, and ``x_Cl-``
- or ``molality`` with explicit ``species=[...]`` and ``solvent=...``

Each ion regression problem must include at least one of:

- ``rho`` or a supported mass-density column
- ``osmotic_coefficient`` or ``osmotic``
- ``mean_ionic_activity``, ``mean_ionic_activity_coefficient``, or ``miac``

``d_born`` declarations additionally require electrostatic provenance: dielectric or relative-permittivity data, ion-activity/osmotic data, or an explicit override. Results include ``result.provenance_report`` so downstream workflows can distinguish supported fitted values from provisional diagnostic values.

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
               P=101325.0,
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
  ``failure_diagnostics``, ``active_bounds``, ``elapsed_seconds``,
  ``partial_pressures``, ``y_vap``,
  ``named_reaction_residuals``, ``source``, ``split``, ``metadata``
- batch: ``success_count``, ``failure_count``, ``row_results``, ``residuals``,
  ``residual_names``, ``residual_row_map``, ``diagnostics``, ``cache_stats``,
  ``timing_summary``
- objective: ``objective``, ``metrics``, and the embedded ``batch_result``

Reporting helpers write those schemas without downstream column guessing:

.. code-block:: python

   epcsaft.write_regression_summary(result, "build/regression/summary.json")
   epcsaft.write_regression_row_table(result, "build/regression/rows.csv")
   epcsaft.write_regression_residual_table(result, "build/regression/residuals.csv")
   epcsaft.write_regression_parameter_table(
       {"Na+.sigma": 2.85},
       "build/regression/parameters.csv",
       seed_map={"Na+.sigma": 2.80},
   )

Reactive-electrolyte parameter fitting is intentionally not a public API until
native Ceres owns that optimizer with exact residual derivatives. Downstream
workflows should use ``evaluate_reactive_regression_objective(...)`` for
residual diagnostics and ``write_regression_parameter_table(...)`` for seed or
candidate parameter maps.

The package-owned micro-benchmark harness for this layer is:

.. code-block:: powershell

   uv run python scripts\benchmarks\benchmark_reactive_regression.py --warmup 3 --repeat 10
   uv run python scripts\benchmarks\benchmark_reactive_regression.py --case reactive_regression_pressure_speciation_35_row_surrogate --warmup 0 --repeat 1

Binary VLE records
------------------

``fit_binary_pair(...)`` V1 supports VLE x/y records only. Records require:

- ``T``: temperature in kelvin
- ``P``: pressure in pascals
- liquid mole-fraction columns such as ``x_H2O`` and ``x_Ethanol``
- vapor mole-fraction columns such as ``y_H2O`` and ``y_Ethanol``

The V1 native optimizer target is constant ``k_ij`` through Ceres with ``cppad_implicit`` Jacobians, including neutral associating binaries where the constant-pressure response combines CppAD explicit EOS terms with association site-fraction implicit density sensitivities. Constant ``l_ij`` and ``k_hb_ij`` remain schema-supported targets, but fitting them raises until a native analytic, CppAD, or implicit derivative path is registered. Linear temperature models and LLE fitting are future phases and raise ``InputError``. Ion-involving binary targets require explicit provenance and are rejected by default unless they are tied to direct electrolyte/neutral-ion data or an explicit override.

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
   print(result.optimizer_backend)
   print(result.derivative_backend)
   print(result.jacobian_backend)
   print(result.fitted_values)
   print(result.initial_parameters)
   print(result.final_parameters)
   print(result.parameter_movement)
   print(result.metrics_by_term)
   print(result.residual_block_norms)
   print(result.target_family_summaries)
   print(result.source_summaries)
   print(result.provenance_report)

   from epcsaft import write_fit_result

   written_paths = write_fit_result(result, template_root, overwrite=False)
   print(written_paths)

With ``overwrite=False``, blank template cells can be filled but existing values are protected. Pure ion fits update the target component row in ``pure/``. Binary fits update both symmetric cells in the relevant interaction matrix.

Derivative and Jacobian access
------------------------------

Normal ``FitResult`` payloads report compact derivative metadata:

- ``optimizer_backend``
- ``derivative_backend``
- ``objective_initial``
- ``objective_final``
- ``initial_parameters``
- ``final_parameters``
- ``parameter_movement``
- ``parameter_map``
- ``active_bounds``
- ``row_diagnostics``
- ``source_summaries``
- ``target_family_summaries``
- ``residual_block_norms``
- ``jacobian_available``
- ``jacobian_backend``

Large matrices are exposed only through explicit derivative-evaluation helpers. Use ``evaluate_pure_neutral_derivatives(...)`` for the native pure-neutral objective. It returns residuals, gradient, ``jacobian_row_major``, and ``jacobian_shape``. Pure-neutral Jacobians use CppAD derivatives plus implicit density sensitivities.

Generic native-record derivative helpers are exposed only after the target family returns exact analytical, CppAD, or implicit derivative matrices from the native Ceres route. Public generic fitting does not route through non-native optimizer loops.

Reactive electrolyte diagnostic objective
-----------------------------------------

Use ``ReactiveElectrolyteRegressionContext.from_batch(...)`` and
``evaluate_reactive_regression_objective(...)`` when a downstream project needs a
fixed-shape diagnostic objective for coupled reactive-electrolyte rows. The
package exposes this as fixed-parameter objective evaluation only; public
reactive-electrolyte parameter fitting stays absent until native Ceres owns the
optimizer and exact derivative path. Downstream code supplies records, targets,
species, balances, reactions, and a ``mixture_factory`` for the current
parameters; the package returns a ``ReactiveRegressionObjectiveResult`` with
residuals, residual names, per-record diagnostics, and success/failure counts.

The diagnostic context keeps the residual vector shape stable when a record
fails by inserting bounded penalty residuals. Successful row diagnostics include
predicted partial pressures, liquid composition, vapor composition, named
reaction residuals, and compact solver diagnostics so downstream code can write
reports without rerunning the same expensive rows. Keep target magnitudes
positive for log-scale pressure and composition residuals:

.. code-block:: python

   batch = epcsaft.ReactiveElectrolyteBatch(
       species=["CO2", "H2O", "MEA", "MEAH+", "HCO3-"],
       rows=rows,
       balances=balances,
       reactions=reactions,
       vapor_species=["CO2", "H2O"],
       mixture_factory=make_mixture_for_candidate,
   )
   result = epcsaft.evaluate_reactive_regression_objective(
       batch,
       objective=epcsaft.ReactiveRegressionObjective(
           residual_weights={"partial_pressure": 1.0, "speciation": 1.0},
           failure_penalty=8.0,
       ),
   )

Use ``result.residuals``, ``result.record_results``, and
``result.diagnostics`` as diagnostic evidence for a fixed candidate parameter
map. Do not treat this helper as parameter fitting or as a full constrained
Gibbs/NLP solve; Ipopt equilibrium routes are explicit native constrained-NLP
routes and are not used automatically.

Derivative availability
-----------------------

.. list-table::
   :header-rows: 1

   * - Method
     - Current Jacobian access
     - Second-derivative exposure
   * - Runtime ``dadt()``, ``dadx()``, ``z(return_contribution_terms=True)``, ``mures(return_contribution_terms=True)``
     - Analytical where available, CppAD where implemented; unsupported derivative paths raise clearly
     - Not exposed
   * - Pure-neutral regression
     - Native CppAD/implicit Jacobian through ``evaluate_pure_neutral_derivatives(...)``
     - Solver-internal only; exact Hessian callbacks are not exposed
   * - Generic ion/binary regression
     - Binary ``k_ij`` fitting uses native Ceres ``cppad_implicit`` Jacobians; other generic residual families raise until analytic or CppAD coverage is implemented
     - Solver-internal only; exact Hessian callbacks are not exposed
   * - Neutral LLE
     - Native Ipopt constrained-NLP route when compiled; exact gradient/Jacobian callbacks are route-internal
     - Ipopt limited-memory Hessian handling is solver-internal
   * - Chemical equilibrium / reactive speciation
     - Explicit native Ipopt ideal-mole-fraction route uses analytic derivatives when Ipopt is compiled; activity/concentration paths raise until EOS derivative NLP blocks exist
     - Ipopt limited-memory Hessian handling is solver-internal

Regression result payloads report optimizer and first-derivative metadata only.
Solver-internal Hessian handling is not exposed as a package derivative backend.
