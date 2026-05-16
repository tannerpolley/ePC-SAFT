Equilibrium Cookbook
====================

Use this page to choose and validate an equilibrium, speciation, or
repeated-property workflow. The public API is Python, but package-owned
thermodynamic solves run through native C++ kernels unless this page explicitly
says otherwise.

Fast validation
---------------

Before changing code that depends on equilibrium behavior, run:

.. code-block:: powershell

   uv run python run_pytest.py --equilibrium-api -q

Use ``uv run python scripts/dev/validate_project.py quick`` before publishing a
change. Use ``--equilibrium-confidence`` only for the bounded electrolyte
confidence checks; full scientific reports remain explicit opt-ins.

Capability and solver selection
-------------------------------

Use ``epcsaft.capabilities()`` to detect what the current install can do. Treat
these fields as routing hints, not as proof that a physical case is valid.

.. list-table::
   :header-rows: 1

   * - Capability
     - Status
     - Use
   * - Neutral TP flash, LLE, stability
     - Production native
     - Default phase-equilibrium workflows for neutral systems.
   * - Electrolyte LLE
     - Production native
     - Fixed-species charge-neutral LLE, preferably after stability checks.
   * - Reactive speciation
     - Production native
     - Homogeneous chemical equilibrium for a single phase, including activity-
       and concentration-coupled standard states.
   * - Electrolyte bubble pressure
     - Production native, scoped
     - Fixed-liquid electrolyte bubble pressure with neutral vapor species only.
   * - Reactive electrolyte bubble
     - Native substeps with fixed-liquid handoff, scoped
     - Native speciation followed by native fixed-liquid electrolyte bubble
       pressure for neutral vapor species.
   * - IPOPT
     - Optional opt-in bridge
     - Bound-constrained residual-minimization refinement only; not full Gibbs/NLP.

.. list-table::
   :header-rows: 1

   * - Solver option
     - Choose when
     - Do not use when
   * - ``solver_backend="auto"``
     - You want the supported native default.
     - You expect IPOPT to run automatically.
   * - Native Newton/default
     - Ordinary equilibrium/speciation solves and continuation.
     - You need active-bound NLP refinement.
   * - Native least squares
     - Package regression helpers.
     - You need a Python optimizer loop.
   * - ``jacobian_backend="auto"``
     - You want the native chemical-equilibrium default: analytic, CppAD, or implicit sensitivities where available, clear failures otherwise.
     - You need strict failure when a specific derivative backend is unavailable.
   * - ``jacobian_backend="cppad"``
     - You need a CppAD residual derivative path and want unsupported routes to fail loudly.
     - You expect a fallback derivative substitute.
   * - ``differential_mode="autodiff"``
     - You need implemented autodiff derivative paths.
     - You need an automatic analytical substitute.
   * - ``solver_backend="ipopt"``
     - You are testing the future native Ipopt constrained-NLP route once it is implemented.
     - You need the current native equilibrium route.

Neutral VLE, LLE, and stability
-------------------------------

Use explicit mixture methods in new code. Use string-dispatched
``mixture.equilibrium(kind=...)`` only as a compatibility shim for older
scripts.

.. code-block:: python

   import numpy as np
   import epcsaft

   mixture = epcsaft.ePCSAFTMixture.from_params(
       {
           "MW": np.asarray([32.042e-3, 84.147e-3]),
           "m": np.asarray([1.5255, 2.5303]),
           "s": np.asarray([3.23, 3.8499]),
           "e": np.asarray([188.9, 278.11]),
           "e_assoc": np.asarray([2899.5, 0.0]),
           "vol_a": np.asarray([0.035176, 0.0]),
           "assoc_scheme": ["2B", None],
           "k_ij": np.asarray([[0.0, 0.051], [0.051, 0.0]]),
       },
       species=["Methanol", "Cyclohexane"],
   )

   z = np.asarray([0.55, 0.45])
   stability = mixture.stability_tp(T=298.15, P=101325.0, z=z)
   lle = mixture.lle_tp(T=298.15, P=101325.0, z=z)
   assert lle.split_detected
   print(lle.phase_labels, lle.diagnostics)

For neutral TP flash diagnostics, check these fields before concluding that a speedup changed solver semantics:

- ``diagnostics["neutral_fast_path"]``
- ``diagnostics["neutral_fallback_used"]``
- ``diagnostics["neutral_fallback_reason"]``

Neutral bubble/dew public methods now raise until native Ipopt route builders
replace the removed Python scalar route.

The fast-path flag means the current native or local-first neutral route handled
the solve directly. The ``neutral_fallback_*`` keys are legacy diagnostics:
they indicate the optimized route used the more conservative neutral path while
preserving the same result contract.

Derivative policy
-----------------

Equilibrium solvers do not execute perturbation-approximated Jacobians. Solver and
result diagnostics report the derivative status explicitly:

- ``thermodynamic_backend``
- ``solver_backend``
- ``derivative_backend``
- ``derivative_status``
- ``not_available_reason``
- ``solved_internal_states``
- ``derivative_backend_by_block``
- ``implicit_sensitivity_blocks``
- ``residual_norm_by_block``
- ``best_state_available``
- ``association_solver_status``

Supported derivative labels are ``analytic``, ``cppad``,
``analytic_implicit`` and ``cppad_implicit``.
Unsupported combinations report ``not_available``. ``auto`` never falls
back to unsupported derivative approximations.

Sequential Reactive Workflow Boundary
-------------------------------------

Reactive equilibrium examples in this package should use fixed or literature
reaction constants first. Create those inputs with
``ReactionDefinition(log_equilibrium_constant=...)`` or the more explicit
``ReactionDefinition.from_literature_constant(...)`` helper, then run chemical
speciation before the phase route.

.. code-block:: python

   reaction = epcsaft.ReactionDefinition.from_literature_constant(
       {"A": -1.0, "B": 1.0},
       log_equilibrium_constant=log_k_literature,
       name="literature_a_to_b",
       standard_state="ideal_mole_fraction",
       source="Smith et al. table 2",
   )

   result = epcsaft.solve_reactive_staged_equilibrium(
       species=["A", "B"],
       mixture_factory=mixture_factory,
       T=298.15,
       P=1.0e5,
       balances={"total": {"A": 1.0, "B": 1.0}},
       totals={"total": 1.0},
       reactions=[reaction],
       initial_x=[0.5, 0.5],
       phase_kind="tp_flash",
   )

The resulting diagnostics identify the sequential coupling level, keep
reaction-constant fitting as ``secondary_optional``, and state that full
simultaneous reactive NLP is not the active route. Fit pure, binary, or
electrolyte ePC-SAFT parameters after this fixed-constant speciation and phase
handoff. Do not make reaction-constant regression a prerequisite for those
parameter fits.

Electrolyte LLE
---------------

Use stability first. If the native LLE solve collapses, provide explicit
charge-neutral ``initial_phases`` and inspect ``diagnostics["seed_attempts"]``.

.. code-block:: python

   import numpy as np
   import epcsaft

   species = ["H2O", "TBP", "[emim][tcb]", "Li+", "Cl-"]
   feed = np.asarray([0.9549141, 0.0290154, 0.00603255, 0.00501896, 0.00501896])
   mixture = epcsaft.ePCSAFTMixture.from_dataset("2024_Hubach", species, feed, 294.15)

   initial_phases = {
       "aq": np.asarray([0.9762254, 0.0147531, 0.00108794, 0.00396628, 0.00396628]),
       "org": np.asarray([0.55, 0.30, 0.10, 0.025, 0.025]),
       "phase_fraction": 0.05,
   }
   result = mixture.electrolyte_lle_tp(
       T=294.15,
       P=101325.0,
       z=feed,
       initial_phases=initial_phases,
       options=epcsaft.EquilibriumOptions(
           max_iterations=180,
           timeout_seconds=15.0,
           max_total_objective_evaluations=5000,
       ),
   )
   print(result.diagnostics)

For diagnostic matrices, use the budget diagnostics rather than killing a
worker process:

.. code-block:: python

   options = epcsaft.EquilibriumOptions(
       max_iterations=24,
       timeout_seconds=8.0,
       max_seed_attempts=4,
       max_total_objective_evaluations=8000,
   )

Keep strict exceptions for single scientific claims. For broad downstream
sweeps where the caller needs to record every row, add
``return_best_effort=True`` and check ``result.split_detected`` plus
``result.diagnostics["acceptance_gate"]`` before using any phases:

.. code-block:: python

   result = mixture.electrolyte_lle_tp(
       T=294.15,
       P=101325.0,
       z=feed,
       options=epcsaft.EquilibriumOptions(
           timeout_seconds=8.0,
           max_seed_attempts=4,
           return_best_effort=True,
       ),
   )
   if not result.split_detected:
       print(result.diagnostics["acceptance_gate"], result.diagnostics["seed_attempts"])

Reactive speciation
-------------------

The caller owns species, balances, totals, reactions, standard states, and
initial composition. Use ``error_mode="result"`` only for diagnostic sweeps.

.. code-block:: python

   import epcsaft

   species = ["H2O", "NaCl", "Na+", "Cl-"]
   mixture = epcsaft.ePCSAFTMixture.from_dataset("2026_Khudaida", species, [0.998, 0.001, 0.0005, 0.0005], 298.15)
   result = epcsaft.solve_reactive_speciation(
       species=species,
       mixture_factory=lambda x, T, P: mixture,
       T=298.15,
       P=101325.0,
       balances={
           "water": {"H2O": 1.0},
           "sodium": {"NaCl": 1.0, "Na+": 1.0},
           "chloride": {"NaCl": 1.0, "Cl-": 1.0},
       },
       totals={"water": 0.998, "sodium": 0.0015, "chloride": 0.0015},
       reactions=[
           epcsaft.ReactionDefinition(
               {"NaCl": -1.0, "Na+": 1.0, "Cl-": 1.0},
               log_equilibrium_constant=-4.0,
               name="salt_dissociation",
           )
       ],
       initial_x=[0.998, 0.001, 0.0005, 0.0005],
       options=epcsaft.ReactiveSpeciationOptions(),
   )
   print(result.x, result.named_reaction_residuals)
   print(result.diagnostics["jacobian_backend"])

With ``jacobian_backend="auto"``, the native chemical-equilibrium residual uses
the analytic log-amount Jacobian for ideal, activity-coupled, and
concentration-coupled standard states. Activity and fugacity terms are
evaluated inside the native residual, not from a cached external activity
vector. Check these fields before treating a result as a production
activity-coupled solve:

- ``diagnostics["solver_language"] == "c++"``
- ``diagnostics["activity_model"]``
- ``diagnostics["reaction_standard_states"]``
- ``diagnostics["activity_or_fugacity_terms_in_residual"]``
- ``diagnostics["derivative_backend_by_block"]["reactive_speciation_variables"]``
- ``diagnostics["implicit_sensitivity_blocks"]``
- ``diagnostics["implicit_solve_results"]["reactive_speciation_variables"]``

Request a specific derivative backend only when unsupported routes should fail
loudly. The required solved composition sensitivity is reported as an implicit
solved-state block rather than by differentiating through solver iterations.

Electrolyte bubble and reactive bubble
--------------------------------------

Use electrolyte bubble pressure only for a fixed liquid composition and neutral
vapor species. Ions remain liquid-only.

.. code-block:: python

   bubble = mixture.electrolyte_bubble_p(
       T=298.15,
       x_liq=[0.998, 0.001, 0.0005, 0.0005],
       vapor_species=["H2O"],
       options=epcsaft.ElectrolyteBubbleOptions(initial_pressure=101325.0),
   )

For reactive electrolyte bubbles, speciation runs first and the equilibrated
liquid goes into the fixed-liquid bubble solve. Keep strict default behavior
for single points; use ``error_mode="result"`` for sweeps that must continue.
This route is suitable for native-backed volatile-neutral partial-pressure
checks from an equilibrated liquid. It is not a simultaneous reactive NLP.

.. code-block:: python

   reactive_bubble = epcsaft.solve_reactive_electrolyte_bubble(
       species=species,
       mixture_factory=lambda x, T, P: mixture,
       T=298.15,
       P_seed=101325.0,
       balances={
           "water": {"H2O": 1.0},
           "sodium": {"NaCl": 1.0, "Na+": 1.0},
           "chloride": {"NaCl": 1.0, "Cl-": 1.0},
       },
       totals={"water": 0.998, "sodium": 0.0015, "chloride": 0.0015},
       reactions=[],
       initial_x=[0.998, 0.001, 0.0005, 0.0005],
       vapor_species=["H2O"],
       options=epcsaft.ReactiveElectrolyteBubbleOptions(error_mode="result"),
   )
   print(reactive_bubble.success, reactive_bubble.diagnostics)

For a pressure/speciation proof, check both nested diagnostic blocks. The
speciation block should report native activity coupling and implicit sensitivity;
the bubble block should report ``native_entrypoint`` as
``_solve_electrolyte_bubble_native`` plus ``partial_pressures`` and
``fugacity_residual_norm`` for the volatile neutral species. Charged species
must be absent from ``y_vap`` unless the caller has explicitly defined a charged
vapor model.

Repeated fugacity/property loops
--------------------------------

Use ``P`` when pressure closure matters. Pass ``rho_guess`` or ``rho_seed`` only
as a density-solve hint. Use direct ``rho`` only when the supplied density is
the physical closure and pressure mismatch is acceptable or separately audited.

.. code-block:: python

   rows = [
       {"T": 298.15, "P": 101325.0, "x": [1.0], "phase": "liq"},
       {"T": 300.15, "P": 101325.0, "x": [1.0], "phase": "liq"},
   ]
   values = epcsaft.evaluate_fugacity_coefficients_batch(mixture, rows=rows)
   next_seed = values[-1]["density"]
   next_value = epcsaft.evaluate_fugacity_coefficients(
       mixture,
       T=301.15,
       P=101325.0,
       x=[1.0],
       rho_seed=next_seed,
   )

IPOPT route status
------------------

``solver_backend="ipopt"`` is reserved for the native Ipopt constrained-NLP
adapter. The old Python adapter has been removed. Until the native adapter is
wired to public equilibrium routes, an explicit IPOPT request raises
``InputError``.

.. code-block:: python

   try:
       result = mixture.electrolyte_lle_tp(
           T=294.15,
           P=101325.0,
           z=feed,
           options=epcsaft.EquilibriumOptions(solver_backend="ipopt"),
       )
   except epcsaft.InputError as exc:
       print("IPOPT was requested before the native adapter was routed:", exc)

