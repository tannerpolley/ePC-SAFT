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

Use ``uv run python scripts/validate_project.py quick`` before publishing a
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
     - Homogeneous chemical equilibrium for a single phase.
   * - Electrolyte bubble pressure
     - Production native, scoped
     - Fixed-liquid electrolyte bubble pressure with neutral vapor species only.
   * - Reactive electrolyte bubble
     - Staged production native, scoped
     - Native speciation followed by fixed-liquid electrolyte bubble pressure.
   * - IPOPT
     - Experimental opt-in
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
     - You want the native chemical-equilibrium default: analytic, CppAD, or implicit sensitivities where available, clear ``backend_unavailable`` diagnostics otherwise.
     - You need strict failure when a specific derivative backend is unavailable.
   * - ``jacobian_backend="autodiff"``
     - You need an existing legacy Eigen forward-mode path and want unsupported routes to fail loudly.
     - You need a fallback to an analytical formula.
   * - ``jacobian_backend="cppad"``
     - You need a CppAD residual derivative path and want unsupported routes to return ``backend_unavailable``.
     - You need a finite-difference fallback.
   * - ``differential_mode="autodiff"``
     - You need implemented autodiff derivative paths.
     - You need a fallback to analytical derivatives.
   * - ``solver_backend="ipopt"``
     - You explicitly installed ``cyipopt`` and want residual-minimization refinement.
     - You need full constrained Gibbs minimization.

Neutral VLE, LLE, and stability
-------------------------------

Use explicit mixture methods in new code. Use string-dispatched
``mixture.equilibrium(kind=...)`` only for compatibility with older scripts.

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

For neutral TP flash and neutral bubble/dew diagnostics, check these fields before concluding that a speedup changed solver semantics:

- ``diagnostics["neutral_fast_path"]``
- ``diagnostics["neutral_fallback_used"]``
- ``diagnostics["neutral_fallback_reason"]``

The fast-path flag means the current native or local-first neutral route handled the solve directly. A fallback flag means the optimized route dropped to the more conservative path while preserving the same result contract.

Derivative policy
-----------------

Equilibrium solvers do not execute finite-difference Jacobians. Solver and
result diagnostics report the derivative status explicitly:

- ``thermodynamic_backend``
- ``solver_backend``
- ``derivative_backend``
- ``derivative_status``
- ``backend_unavailable_reason``
- ``solved_internal_states``
- ``derivative_backend_by_block``
- ``implicit_sensitivity_blocks``
- ``residual_norm_by_block``
- ``best_state_available``
- ``association_solver_status``

Supported derivative labels are ``analytic``, ``cppad``,
``analytic_implicit``, ``cppad_implicit``, and ``legacy_eigen_forward``.
Unsupported combinations report ``backend_unavailable``. ``auto`` never falls
back to finite differences.

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

With ``jacobian_backend="auto"``, ideal mole-fraction standard states use the
analytic native Jacobian. Activity- or concentration-coupled standard states
raise ``backend_unavailable`` until analytic or autodiff residual derivatives
are implemented. Request ``jacobian_backend="autodiff"`` only when unsupported
derivative paths should fail loudly.

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

IPOPT opt-in
------------

``solver_backend="ipopt"`` requires the optional ``cyipopt`` dependency. If it
is not importable, the package raises ``InputError`` and does not silently fall
back to Newton.

.. code-block:: python

   try:
       result = mixture.electrolyte_lle_tp(
           T=294.15,
           P=101325.0,
           z=feed,
           options=epcsaft.EquilibriumOptions(solver_backend="ipopt"),
       )
   except epcsaft.InputError as exc:
       print("IPOPT was requested but is unavailable:", exc)

