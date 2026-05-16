Electrolyte VLE And Reactive Speciation
=======================================

Native-default solver policy
----------------------------

Package-owned regression and equilibrium solves must run through the native
C++ backend exposed by ``epcsaft._core``. Python package code may normalize
inputs, build request payloads, validate units and compositions, and convert
native payloads into structured result objects. New production equilibrium
ownership is moving to native Ipopt constrained NLPs; existing native residual
routes remain transitional until each Ipopt route builder lands.

Consequently, electrolyte bubble-pressure and composed reactive electrolyte
bubble-pressure entry points now fail loudly until native Ipopt route builders
own those solves. The Python layer may coordinate native speciation, native
phase-equilibrium calls, structured results, and diagnostics only after those
native route builders exist.

``ReactiveSpeciationOptions`` accepts ``solver_backend="auto" | "ipopt"``.
``auto`` keeps the current native chemical-equilibrium solve while the full
migration is underway. Explicit ``ipopt`` routes homogeneous
``ideal_mole_fraction`` speciation through the native Ipopt constrained-NLP
adapter when the extension is built with Ipopt; activity and concentration
standard states still require the EOS derivative NLP blocks.

The target IPOPT formulation is constrained Gibbs minimization with formal
material, charge, and reaction constraints plus exact gradients/Jacobians from
analytic or CppAD-backed callbacks. Limited-memory Hessian handling is allowed
only as Ipopt solver-internal behavior and is not reported as a package
derivative backend.

Solver-selection policy is intentionally conservative. ``auto`` keeps the native
chemical-equilibrium solver until the corresponding native Ipopt route is
implemented and validated.

Homogeneous reactive speciation
-------------------------------

``solve_reactive_speciation(...)`` solves one homogeneous reactive phase using
the C++ chemical-equilibrium kernel. The caller owns the chemistry: species
labels, material balances, totals, reactions, equilibrium constants, and the
``mixture_factory`` used to create the native ePC-SAFT mixture.

The near-term reactive workflow is sequential and fixed-constant first:

1. pass fixed or literature reaction constants through
   ``ReactionDefinition.log_equilibrium_constant`` or
   ``ReactionDefinition.from_literature_constant(...)``;
2. evaluate activity-coupled reactive speciation with ePC-SAFT states;
3. hand the equilibrated composition to phase or electrolyte-equilibrium
   routes; and
4. regress ePC-SAFT pure, binary, or electrolyte parameters against that fixed
   chemistry boundary.

Reaction-constant fitting is therefore an optional later refinement, not the
default workflow and not a blocking dependency for pure, binary, or electrolyte
parameter regression.

The native backend solves material balances, charge balance, and reaction
residuals in log mole amounts. Activity coefficients are evaluated only when
the selected reaction standard state needs them or when
``ReactiveSpeciationOptions.activity_output="always"`` is requested. This keeps
``standard_state="concentration"`` and ``"ideal_mole_fraction"`` workflows from
paying for unused activity calls under the default ``activity_output="auto"``.
The derivative default is native-owned and diagnostic-friendly:
``jacobian_backend="auto"`` uses the analytic log-amount Jacobian for
``standard_state="ideal_mole_fraction"``. Activity- or concentration-coupled
standard states raise until analytic, CppAD, or implicit residual derivatives
are implemented. Diagnostics report
``requested_jacobian_backend``, ``derivative_backend``, and ``derivative_status``
for strict route-boundary failures. Approximate Jacobian substitutes are not
supported.

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
``error_mode="result"`` or ``return_best_effort=True`` to receive a
``ReactiveSpeciationResult(success=False, ...)`` from the best native
payload instead. Diagnostics include ``best_x``,
``best_activity_coefficients``, family residual norms, and named reaction
residuals. ``mass_tolerance``, ``charge_tolerance``, and
``reaction_tolerance`` may be set separately; each defaults to ``tolerance``
when omitted.

For repeated nearby states, use ``solve_reactive_speciation_sweep(...)`` with
``continuation="auto"``. Each successful result carries a
``continuation_state`` containing the composition and lightweight counters for
the next point. Failed points return fixed-shape diagnostic result objects when
``error_mode="result"`` is selected, so downstream sweeps can continue.
Speciation diagnostics report iteration, residual, Jacobian, density, state,
and activity-evaluation counts.

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
``solve_reactive_electrolyte_bubble(...)`` are reserved public contract names.
Their previous fixed-liquid pressure-search implementation is disabled
because production equilibrium routes must be native Ipopt NLPs, not
package-owned scalar pressure searches.

The target reactive electrolyte bubble route will first solve homogeneous
chemical equilibrium with the native chemical solver and then hand the
equilibrated liquid composition to the native Ipopt electrolyte bubble solver.
``ReactiveElectrolyteBubbleResult`` is retained as the structured result shape
for that route:

- ``speciation_strict_success`` reflects the standalone native chemical
  equilibrium tolerances requested in ``ReactiveSpeciationOptions``.
- ``speciation_phase_handoff_success`` reflects whether finite mass, charge,
  and reaction residuals are accurate enough for phase-equilibrium handoff.
- ``bubble_success`` will reflect the native Ipopt electrolyte bubble-pressure solve.

The handoff tolerances default to ``1e-8`` for mass, ``1e-8`` for charge, and
``1e-5`` for reaction residuals. They can be adjusted with
``ReactiveElectrolyteBubbleOptions.phase_handoff_mass_tolerance``,
``phase_handoff_charge_tolerance``, and
``phase_handoff_reaction_tolerance``. This lets downstream workflows keep
strict standalone chemical-equilibrium tolerances while accepting a looser,
documented phase-equilibrium handoff envelope when the bubble solve itself is
well converged.

Strict default behavior currently raises ``InputError`` for the route-pending
bubble stage. Diagnostic result-mode sweeps become meaningful again after the
native Ipopt bubble route builder lands.

Solved-state derivative boundary
--------------------------------

Sequential reactive workflows contain nested solved states: association site
fractions, reactive speciation variables, density roots, bubble-pressure roots,
and phase-equilibrium variables. These blocks must report derivatives through
analytic residual derivatives, CppAD residual partials, implicit sensitivities,
or strict route-boundary failures when coverage is incomplete. Approximate
derivative substitutes are not supported.

For active association, the preferred sequential derivative boundary is to solve
the association site fractions, compute residual partials, solve the implicit
sensitivity system, and propagate those sensitivities into activity, fugacity,
or chemical-potential derivatives. Direct CppAD through the association
fixed-point iteration is not required for the sequential fixed-constant workflow.

Future simultaneous root systems, Ceres residual systems, or Ipopt constrained
NLP formulations may be useful for harder coupled systems, but they are future
architecture options. They are not the default fixed-literature-constant
workflow and are not required for reactive electrolyte VLE handoffs.
