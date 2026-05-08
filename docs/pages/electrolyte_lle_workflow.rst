Electrolyte LLE Workflow
========================

Hard electrolyte LLE cases should be handled as fixed-species phase-equilibrium
problems before chemical or speciation equilibrium is introduced.

Recommended workflow
--------------------

1. Run ``kind="electrolyte_stability"`` to confirm whether the feed is unstable.
2. If ``min_tpd`` is negative and ``kind="electrolyte_lle"`` collapses, provide
   explicit charge-neutral ``initial_phases``.
3. Check ``diagnostics["seed_attempts"]`` to see which seed families were tried.
4. For curves, solve one point with strong initial phases and use continuation
   through ``equilibrium_curve``.

Bounded diagnostic runs
-----------------------

Long downstream sweeps should bound hard fixed-species LLE attempts with
``EquilibriumOptions`` rather than relying on an outer process timeout. The
default values are ``None`` and keep the existing robust solve behavior. When a
budget is provided and exhausted, strict mode raises ``SolutionError`` with
structured diagnostics instead of hanging indefinitely.

Useful budget controls are:

- ``timeout_seconds`` for native wall-clock budget checks inside the LLE route.
- ``max_seed_attempts`` to stop after a bounded number of seed families.
- ``max_density_failures`` to stop density-heavy diagnostics when repeated
  phase-state construction fails.
- ``max_total_objective_evaluations`` to bound transformed Gibbs objective
  work in exploratory matrices.

On budget stop, diagnostics include ``acceptance_gate =
"predictive_budget_exhausted"``, ``budget_exceeded``, ``budget_trigger``,
``elapsed_seconds``, ``objective_evaluation_count``, and the requested budget
values. Use these fields to distinguish "the case is physically hard or
unaccepted" from "the calling script killed the process before the package
could report diagnostics."

Strict versus best-effort results
---------------------------------

Strict behavior remains the default. If fixed-species electrolyte LLE does not
pass the predictive acceptance gates, the Python binding raises
``SolutionError`` with JSON-safe diagnostics.

For downstream sweeps, set ``return_best_effort=True`` to receive a structured
``EquilibriumResult`` instead of an exception when the native route reaches a
finite but unaccepted diagnostic state:

.. code-block:: python

   result = mix.equilibrium(
       kind="electrolyte_lle",
       T=294.15,
       P=1.013e5,
       z=feed,
       options=epcsaft.EquilibriumOptions(
           timeout_seconds=8.0,
           max_seed_attempts=4,
           return_best_effort=True,
       ),
   )

   if not result.split_detected:
       gate = result.diagnostics["acceptance_gate"]
       seed_attempts = result.diagnostics["seed_attempts"]

Best-effort mode does not weaken acceptance. A returned result with
``split_detected=False`` is a diagnostic payload, not a solved phase split.
When the native route has a noncollapsed best candidate, diagnostics report
``best_effort_phases_returned=True`` and include the candidate phases for
inspection or continuation experiments; otherwise the result carries the same
failure diagnostics without phases.

Optional IPOPT surface
----------------------

``EquilibriumOptions`` accepts ``solver_backend="auto" | "newton" | "ipopt"``
and ``hessian_strategy="gauss_newton" | "lbfgs"``. The default ``auto`` keeps
the current native Newton path. ``solver_backend="ipopt"`` is explicit opt-in
and requires the optional ``cyipopt`` dependency; if ``cyipopt`` is unavailable,
the package raises ``InputError`` and does not silently fall back to Newton.

The cyipopt adapter solves a bounded transformed-variable min-residual NLP for
electrolyte LLE by calling native residual/Jacobian callbacks. Acceptance still
uses the package engineering gates: residual norm, material balance, charge
balance, phase distance, stability diagnostics, and no collapsed split.

This is not yet a full constrained thermodynamic NLP. Material and charge
balances are checked by acceptance gates and residual penalties rather than
being exposed to IPOPT as equality constraints. Diagnostics separate
``ipopt_success`` from ``residual_gate_success``, ``physical_gate_success``, and
``accepted``. Approximate Hessian choices are reported as ``hessian_strategy``
and ``hessian_kind``; ``exact_hessian_available`` remains ``False``.

When ``solver_backend="ipopt"`` is requested without explicit
``initial_phases``, the adapter first asks the native transformed Newton route
for a seed and then refines that split with IPOPT. If that seed cannot be
generated, diagnostics report ``ipopt_seed_failure`` and the residual evaluator
falls back to its native default seed.

Solver-selection policy
-----------------------

``solver_backend="auto"`` remains conservative and does not switch to IPOPT just
because ``cyipopt`` is installed. Use IPOPT explicitly when active bounds or a
near-bound residual-minimization refinement matters. Prefer the native Newton
route for ordinary fixed-species electrolyte LLE solves and continuation unless
the diagnostics show a bound or split-acceptance failure worth refining.

Hubach-style example
--------------------

.. code-block:: python

   import numpy as np
   import epcsaft
   from epcsaft import ePCSAFTMixture

   species = ["H2O", "TBP", "[emim][tcb]", "Li+", "Cl-"]
   feed = np.asarray([
       0.9549140976171719,
       0.029015432040304778,
       0.006032546048615457,
       0.005018962146953909,
       0.005018962146953909,
   ])
   initial_phases = {
       "aq": np.asarray([0.9762253659128125, 0.014753086358215556, 0.0010879432090689022, 0.003966275944161999, 0.003966275944161999]),
       "org": np.asarray([0.55, 0.30, 0.10, 0.025, 0.025]),
       "phase_fraction": 0.05,
   }

   mix = ePCSAFTMixture.from_dataset("2024_Hubach", species, feed, 294.15)
   result = mix.equilibrium(
       kind="electrolyte_lle",
       T=294.15,
       P=1.013e5,
       z=feed,
       initial_phases=initial_phases,
       options=epcsaft.EquilibriumOptions(
           max_iterations=180,
           tolerance=1.0e-8,
           timeout_seconds=15.0,
           max_total_objective_evaluations=5000,
       ),
   )

   assert result.split_detected
   assert result.diagnostics["phase_distance"] > 1.0e-4
