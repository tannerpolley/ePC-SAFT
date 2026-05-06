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
   )

   print(result.x)
   print(result.mass_balance_residuals)
   print(result.reaction_residuals)

The solver uses a small damped finite-difference least-squares loop and does
not require SciPy. It is intended as a reusable building block for project-owned
reaction packages such as MEA absorbers rather than a hard-coded chemistry
model.
