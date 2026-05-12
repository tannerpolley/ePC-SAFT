Electrolyte VLE And Reactive Speciation
=======================================

Native-default solver policy
----------------------------

Package-owned regression and equilibrium solves must run through the native
C++ backend exposed by ``epcsaft._core``. Python package code may normalize
inputs, build request payloads, validate units and compositions, and convert
native payloads into structured result objects. The default production paths
remain native Newton and native least-squares.

Consequently, electrolyte bubble-pressure and composed reactive electrolyte
bubble-pressure entry points route through native C++ kernels. The Python layer
only coordinates native speciation, native phase-equilibrium calls, fixed-liquid
bubble-pressure handoffs, structured results, and diagnostics.

``ReactiveSpeciationOptions`` accepts ``solver_backend="auto" | "newton" |
"ipopt"`` and ``hessian_strategy="gauss_newton" | "lbfgs"``. ``auto`` keeps
the current native chemical-equilibrium solve. ``ipopt`` is explicit opt-in,
uses the optional ``cyipopt`` package, and raises ``InputError`` when requested
without ``cyipopt`` instead of falling back to Newton. The adapter solves a
bounded log-amount residual-minimization NLP through native residual/Jacobian
callbacks and preserves mass, charge, and reaction residual diagnostics.

The current IPOPT formulation is experimental residual minimization, not full
constrained Gibbs minimization. Diagnostics report ``formulation``,
``ipopt_success``, residual/physical gate status, and approximate Hessian
metadata separately. ``hessian_strategy="gauss_newton"`` supplies a
least-squares ``J.T @ J`` callback; ``hessian_strategy="lbfgs"`` delegates to
IPOPT limited-memory Hessian approximation. Neither route includes exact second
residual derivatives.

Solver-selection policy is intentionally conservative. ``auto`` keeps the native
chemical-equilibrium solver. IPOPT is appropriate only as an explicit
bound-constrained residual-minimization experiment until material, charge, and
reaction constraints are exposed as formal NLP constraints.

Homogeneous reactive speciation
-------------------------------

``solve_reactive_speciation(...)`` solves one homogeneous reactive phase using
the C++ chemical-equilibrium kernel. The caller owns the chemistry: species
labels, material balances, totals, reactions, equilibrium constants, and the
``mixture_factory`` used to create the native ePC-SAFT mixture.

The native backend solves material balances, charge balance, and reaction
residuals in log mole amounts. Activity coefficients are evaluated only when
the selected reaction standard state needs them or when
``ReactiveSpeciationOptions.activity_output="always"`` is requested. This keeps
``standard_state="concentration"`` and ``"ideal_mole_fraction"`` workflows from
paying for unused activity calls under the default ``activity_output="auto"``.
The derivative default is native-owned and diagnostic-friendly:
``jacobian_backend="auto"`` uses the analytic log-amount Jacobian for
``standard_state="ideal_mole_fraction"``. Activity- or concentration-coupled
standard states currently raise ``backend_unavailable`` until analytic or
autodiff residual derivatives are implemented. Diagnostics report
``requested_jacobian_backend``, ``derivative_backend``, ``derivative_status``,
and ``backend_unavailable_reason`` so users can see the active support boundary.
Explicit ``jacobian_backend="autodiff"`` remains strict and raises when the
requested derivative path is unavailable.

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

Strict default behavior still raises on bubble-stage failure. For diagnostic
sweeps, set ``ReactiveElectrolyteBubbleOptions(error_mode="result")`` to return
a structured ``ReactiveElectrolyteBubbleResult(success=False, ...)`` containing
the successful speciation payload plus the failed bubble diagnostics. Complete
examples are in :doc:`equilibrium_cookbook`.
