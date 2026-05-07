Electrolyte VLE And Reactive Speciation
=======================================

Native-only solver policy
-------------------------

Package-owned regression and equilibrium solves must run through the native
C++ backend exposed by ``epcsaft._core``. Python package code may normalize
inputs, build request payloads, validate units and compositions, and convert
native payloads into structured result objects, but it must not contain
production nonlinear equilibrium or regression optimizers.

Consequently, electrolyte bubble-pressure and composed reactive electrolyte
bubble-pressure entry points route through native C++ kernels. The Python layer
only coordinates native speciation, native phase-equilibrium calls, fixed-liquid
bubble-pressure handoffs, structured results, and diagnostics.

Homogeneous reactive speciation
-------------------------------

``solve_reactive_speciation(...)`` solves one homogeneous reactive phase using
the C++ chemical-equilibrium kernel. The caller owns the chemistry: species
labels, material balances, totals, reactions, equilibrium constants, and the
``mixture_factory`` used to create the native ePC-SAFT mixture.

The native backend solves material balances, charge balance, and reaction
residuals in log mole amounts. For ion-containing mixtures it evaluates
ePC-SAFT component activity coefficients. For nonionic mixtures it uses neutral
ePC-SAFT fugacity-reference activities.

.. code-block:: python

   result = epcsaft.solve_reactive_speciation(
       species=["H2O", "NaCl", "Na+", "Cl-"],
       mixture_factory=lambda x, T, P: mixture,
       T=298.15,
       P=1.0e5,
       balances={
           "water_total": {"H2O": 1.0},
           "sodium_total": {"NaCl": 1.0, "Na+": 1.0},
           "chloride_total": {"NaCl": 1.0, "Cl-": 1.0},
       },
       totals={"water_total": 0.998, "sodium_total": 0.0015, "chloride_total": 0.0015},
       reactions=[
           epcsaft.ReactionDefinition(
               stoichiometry={"NaCl": -1.0, "Na+": 1.0, "Cl-": 1.0},
               log_equilibrium_constant=-4.0,
               name="salt_dissociation",
           )
       ],
       initial_x=[0.998, 0.001, 0.0005, 0.0005],
       options=epcsaft.ReactiveSpeciationOptions(),
   )

   print(result.x)
   print(result.mass_balance_residuals)
   print(result.named_reaction_residuals)

``ReactiveSpeciationOptions`` keeps strict failure behavior by default. Failed
native solves raise ``SolutionError``. For diagnostic grid work, opt into
``return_best_effort=True`` to receive a
``ReactiveSpeciationResult(success=False, ...)`` from the best finite native
payload instead. Diagnostics include ``best_x``,
``best_activity_coefficients``, family residual norms, and named reaction
residuals. ``mass_tolerance``, ``charge_tolerance``, and
``reaction_tolerance`` may be set separately; each defaults to ``tolerance``
when omitted.

Reactive stability handoff
--------------------------

``kind="reactive_stability"`` first solves homogeneous chemical equilibrium
with the native chemical solver, then sends the equilibrated composition to the
existing native stability route.

.. code-block:: python

   result = mixture.equilibrium(
       kind="reactive_stability",
       T=298.15,
       P=1.013e5,
       z=feed_guess,
       balances=balances,
       totals=totals,
       reactions=reactions,
       options=epcsaft.ReactiveSpeciationOptions(),
   )

The returned ``StabilityResult`` diagnostics include
``reactive_chemical_equilibrium`` and ``reactive_feed_composition`` so callers
can inspect the chemical-equilibrated feed used for tangent-plane-distance
analysis. This route is a native stability/handoff coordinator, not a full
rigorous reactive flash solver.

Electrolyte bubble pressure
---------------------------

``kind="electrolyte_bubble_pressure"`` and
``solve_reactive_electrolyte_bubble(...)`` are native-backed fixed-liquid
bubble-pressure workflows. Ions remain liquid-only. Vapor fugacity is evaluated
with a neutral vapor submixture built from the declared vapor species.
These paths do not fall back to Python pressure, speciation, or regression
loops.

For reactive electrolyte bubbles, the package first solves homogeneous
chemical equilibrium with the native chemical solver and then hands the
equilibrated liquid composition to the native electrolyte bubble solver.
``ReactiveElectrolyteBubbleResult`` reports both strict speciation success and
phase-handoff success:

- ``speciation_strict_success`` reflects the standalone native chemical
  equilibrium tolerances requested in ``ReactiveSpeciationOptions``.
- ``speciation_phase_handoff_success`` reflects whether finite mass, charge,
  and reaction residuals are accurate enough for phase-equilibrium handoff.
- ``bubble_success`` reflects the native electrolyte bubble-pressure solve.

The handoff tolerances default to ``1e-8`` for mass, ``1e-8`` for charge, and
``1e-5`` for reaction residuals. They can be adjusted with
``ReactiveElectrolyteBubbleOptions.phase_handoff_mass_tolerance``,
``phase_handoff_charge_tolerance``, and
``phase_handoff_reaction_tolerance``. This lets downstream workflows keep
strict standalone chemical-equilibrium tolerances while accepting a looser,
documented phase-equilibrium handoff envelope when the bubble solve itself is
well converged.
