Electrolyte VLE And Reactive Speciation
=======================================

The first package-level electrolyte VLE route is a bubble-pressure calculation
for an ion-containing liquid phase with neutral vapor species. Ions remain
liquid-only in this V1 API, and vapor fugacity is evaluated with a neutral
vapor-only submixture built from the same pure and binary parameters.

Electrolyte bubble pressure
---------------------------

Use ``kind="electrolyte_bubble_pressure"`` when the liquid composition is a
fixed true-species electrolyte composition and the vapor contains only declared
neutral volatile species.

.. code-block:: python

   import numpy as np
   import epcsaft

   result = mixture.equilibrium(
       kind="electrolyte_bubble_pressure",
       T=313.15,
       x_liq=np.asarray([0.02, 0.979, 0.0005, 0.0005]),
       volatile_species=["CO2", "H2O"],
       vapor_species=["CO2", "H2O"],
       nonvolatile_species=["Na+", "Cl-"],
       options=epcsaft.ElectrolyteBubbleOptions(initial_pressure=1.0e5),
   )

   assert result.success
   print(result.P)
   print(result.y_vap)
   print(result.partial_pressures)

``z`` is accepted as an alias only when ``x_liq`` is omitted. The liquid
composition must be charge neutral, and ``vapor_species`` must not include ions.
Nonconverged solves raise ``SolutionError`` with JSON-like diagnostics containing
the pressure bracket, iteration history, and state failure count.

Reactive speciation
-------------------

``solve_reactive_speciation(...)`` solves one homogeneous reactive phase using
log mole amounts. The caller owns the chemistry: species labels, material
balances, totals, reactions, equilibrium constants, and the ``mixture_factory``
used to evaluate activity coefficients.

The package chemical-equilibrium implementation is the C++ homogeneous
chemical-equilibrium kernel. The native backend solves material balances,
charge balance, and reaction residuals in log mole amounts, evaluates ePC-SAFT
component activity coefficients for ion-containing mixtures, and uses neutral
ePC-SAFT fugacity-reference activities for nonionic mixtures. This is the fast
inner primitive intended for MEA/CO2/H2O speciation packages and for
Ascani-style reactive phase-equilibrium workflows.

For the native backend, the final activity-coupled residual solve is warm
started by an ideal convex Gibbs minimization in reaction-extent space. The
ideal soft start is only used when it produces a finite composition that
improves the full ePC-SAFT residual; otherwise the solver falls back to the
caller-provided initial composition. The runtime ladder is therefore:
ideal convex Gibbs soft start, full ePC-SAFT activity residual Newton solve,
then the phase-equilibrium handoff diagnostics.

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
           )
       ],
       initial_x=[0.998, 0.001, 0.0005, 0.0005],
       options=epcsaft.ReactiveSpeciationOptions(),
   )

   print(result.x)
   print(result.mass_balance_residuals)
   print(result.reaction_residuals)

The native solver does not require SciPy at runtime. Its diagnostics include a
``phase_equilibrium_handoff`` payload with the equilibrated composition and
activity basis. That handoff is deliberately shaped for the existing native
stability, LLE, and VLE routes: chemical-equilibrate the feed, run phase
stability or a non-reactive flash, chemical-equilibrate any phase seeds, then
use those seeds in a rigorous reactive phase calculation as that layer is
added.

The same route is also available from the mixture equilibrium facade:

.. code-block:: python

   result = mixture.equilibrium(
       kind="chemical_equilibrium",
       T=298.15,
       P=1.0e5,
       z=[0.998, 0.001, 0.0005, 0.0005],
       balances={...},
       totals={...},
       reactions=[...],
       options=epcsaft.ReactiveSpeciationOptions(),
   )

Reactive stability handoff
--------------------------

``kind="reactive_stability"`` runs the first two steps of the Ascani-style
workflow from the public facade: it first solves the homogeneous chemical
equilibrium with the native chemical solver, then sends the equilibrated
composition to the existing native stability route.

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
analysis. This route is intentionally a stability/handoff coordinator, not yet
a full rigorous reactive flash solver.
