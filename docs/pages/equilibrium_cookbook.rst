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
   * - Neutral TP flash and LLE
     - Native Ipopt route when compiled
     - Uses native constrained-NLP route builders; no Python solve loop exists.
   * - Neutral bubble/dew pressure and temperature
     - Native Ipopt route when compiled
     - Uses native constrained-NLP route builders; no Python solve loop exists.
   * - Neutral stability
     - Native Ipopt route required
     - Fails loudly until native Ipopt stability route builders own production use.
   * - Electrolyte LLE
     - Native Ipopt route when compiled
     - Requires an Ipopt-enabled build; no alternate public optimizer path exists.
   * - Reactive speciation
     - Explicit Ipopt ideal route
     - Homogeneous ``ideal_mole_fraction`` chemical equilibrium when Ipopt is
       compiled; activity and concentration standard states are route-gated.
   * - Electrolyte bubble pressure
     - Native Ipopt route when compiled
     - Fixed liquid composition with neutral vapor species; ions remain liquid-only.
   * - Reactive electrolyte bubble
     - Staged native route when compiled
     - Uses native speciation followed by the native fixed-liquid electrolyte
       bubble route for scoped supported inputs.
   * - IPOPT
     - Optional native constrained-NLP backend
     - Owns implemented equilibrium routes when the extension is compiled.

.. list-table::
   :header-rows: 1

   * - Solver option
     - Choose when
     - Do not use when
   * - ``solver_backend="auto"``
     - You want the supported native default.
     - You expect a route without native Ipopt ownership to run.
   * - ``jacobian_backend="auto"``
     - You want the native chemical-equilibrium default: analytic, CppAD, or implicit sensitivities where available, clear failures otherwise.
     - You need strict failure when a specific derivative backend cannot run.
   * - ``jacobian_backend="cppad"``
     - You need a CppAD residual derivative path and want unsupported routes to fail loudly.
     - You expect a substitute derivative backend.
   * - ``solver_backend="ipopt"``
     - You are testing an implemented native Ipopt constrained-NLP route.
     - You need the current native equilibrium route.

Neutral VLE, LLE, and stability
-------------------------------

Use typed problem objects or explicit mixture methods in new code.
``mixture.equilibrium(kind=...)`` remains available for string-dispatched
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
   caps = epcsaft.capabilities()
   if caps["equilibrium"]["neutral_lle_flash"]["available"]:
       lle = mixture.lle_tp(T=298.15, P=101325.0, z=z)
       assert lle.split_detected
       print(lle.phase_labels, lle.diagnostics)
   else:
       print("neutral LLE needs an Ipopt-enabled native build")

Neutral TP flash, neutral LLE, fixed-temperature bubble/dew pressure, and
fixed-pressure bubble/dew temperature route through native Ipopt when that
extension is compiled. Stability remains native-Ipopt-gated.
Downstream sweeps should inspect ``epcsaft.capabilities()`` or catch typed
route errors instead of manufacturing alternate phase-equilibrium results.

Derivative policy
-----------------

Equilibrium solvers do not execute perturbation-approximated Jacobians. Solver and
result diagnostics report the derivative backend explicitly:

- ``thermodynamic_backend``
- ``solver_backend``
- ``derivative_backend``
- ``derivative_backend_by_block``
- ``implicit_sensitivity_blocks``
- ``residual_norm_by_block``

Supported derivative labels are ``analytic``, ``cppad``,
``analytic_implicit`` and ``cppad_implicit``.
Unsupported combinations raise at the route boundary. ``auto`` never switches
to substitute derivative approximations.

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

Electrolyte stability remains route-gated until a native Ipopt stability NLP
builder owns that calculation. Call the native electrolyte LLE route with the
feed specification only; the native Ipopt route builder owns the canonical
initial point and user-provided phase seeds are rejected.

.. code-block:: python

   import numpy as np
   import epcsaft

   species = ["H2O", "TBP", "[emim][tcb]", "Li+", "Cl-"]
   feed = np.asarray([0.9549141, 0.0290154, 0.00603255, 0.00501896, 0.00501896])
   mixture = epcsaft.ePCSAFTMixture.from_dataset("2024_Hubach", species, feed, 294.15)

   result = mixture.electrolyte_lle_tp(
       T=294.15,
       P=101325.0,
       z=feed,
       options=epcsaft.EquilibriumOptions(
           max_iterations=180,
           timeout_seconds=15.0,
       ),
   )
   print(result.diagnostics)

For diagnostic matrices, keep strict package diagnostics visible and bound the
native route with ``timeout_seconds`` rather than reintroducing Python-side
search controls:

.. code-block:: python

   options = epcsaft.EquilibriumOptions(
       max_iterations=24,
       timeout_seconds=8.0,
   )

Keep strict exceptions for single scientific claims. For broad downstream
sweeps where the caller needs to record every row, catch ``SolutionError`` and
store ``exc.args[1]`` diagnostics before moving to the next row:

.. code-block:: python

   try:
       result = mixture.electrolyte_lle_tp(
           T=294.15,
           P=101325.0,
           z=feed,
           options=epcsaft.EquilibriumOptions(
               timeout_seconds=8.0,
           ),
       )
   except epcsaft.SolutionError as exc:
       diagnostics = exc.args[1]
       print(diagnostics["acceptance_gate"])

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
the analytic derivative path for the native-Ipopt ``ideal_mole_fraction`` route.
Activity- and concentration-coupled routes raise until their EOS derivative NLP
blocks exist. Check these fields before treating a result as a production
reactive-speciation solve:

- ``diagnostics["solver_language"] == "c++"``
- ``diagnostics["activity_model"]``
- ``diagnostics["reaction_standard_states"]``
- ``diagnostics["derivative_backend"]``
- ``diagnostics["ipopt_solver_ran"]``
- ``diagnostics["ipopt_accepted"]``

Request a specific derivative backend only when unsupported routes should fail
loudly. Ipopt limited-memory Hessian behavior is solver-internal and is not a
package derivative backend.

Electrolyte bubble and reactive bubble
--------------------------------------

Electrolyte bubble pressure uses the native Ipopt fixed-liquid bubble-pressure
route when Ipopt is compiled. Reactive electrolyte bubble pressure uses staged
native speciation followed by the native fixed-liquid bubble route for scoped
supported inputs.

.. code-block:: python

   caps = epcsaft.capabilities()
   assert caps["equilibrium"]["electrolyte_bubble_pressure"]["available"] == caps["optimizers"]["ipopt"]["available"]

Repeated fugacity/property loops
--------------------------------

Use ``P`` when pressure closure matters. Pass ``rho_guess`` only as a
density-solve hint. Use direct ``rho`` only when the supplied density is the
physical closure and pressure mismatch is acceptable or separately audited.

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
       rho_guess=next_seed,
   )

IPOPT route status
------------------

``solver_backend="ipopt"`` is reserved for the native Ipopt constrained-NLP
adapter. The old Python adapter has been removed. Implemented public routes use
the native adapter when Ipopt is compiled; methods without native Ipopt
ownership raise ``InputError``.

.. code-block:: python

   try:
       result = mixture.electrolyte_lle_tp(
           T=294.15,
           P=101325.0,
           z=feed,
           options=epcsaft.EquilibriumOptions(solver_backend="ipopt"),
       )
   except epcsaft.InputError as exc:
       print("IPOPT route is dependency-gated in this build or method:", exc)

