Electrolyte LLE Workflow
========================

Hard electrolyte LLE cases should be handled as fixed-species phase-equilibrium
problems before chemical or speciation equilibrium is introduced.

Recommended workflow
--------------------

1. Run ``kind="electrolyte_stability"`` to confirm whether the feed is unstable.
2. If ``min_tpd`` is negative, provide explicit charge-neutral
   ``initial_phases`` so the future native Ipopt route has a well-defined
   request payload.
3. For curves, solve one point with accepted phase data and use continuation
   only after the native route builder owns that production path.

Bounded diagnostic runs
-----------------------

Long downstream sweeps should keep strict package failures visible. The public
facade no longer accepts old electrolyte LLE seed-family or density-budget
controls. ``timeout_seconds`` remains available as the wall-clock option for
the native Ipopt route, while unsupported route requests raise typed package
errors instead of running Python-side search logic.

Strict failure results
----------------------

If fixed-species electrolyte LLE does not pass the predictive acceptance gates,
the Python binding raises ``SolutionError`` with JSON-safe diagnostics. Broad
downstream sweeps should catch that exception and record the diagnostics rather
than asking the package to return an unaccepted phase result.

.. code-block:: python

   try:
       result = mix.equilibrium(
           kind="electrolyte_lle",
           T=294.15,
           P=1.013e5,
           z=feed,
           options=epcsaft.EquilibriumOptions(
               timeout_seconds=8.0,
           ),
       )
   except epcsaft.SolutionError as exc:
       diagnostics = exc.args[1]
       gate = diagnostics["acceptance_gate"]

Native IPOPT plan
-----------------

``EquilibriumOptions`` accepts ``solver_backend="auto" | "ipopt"``. The
default ``auto`` validates the public request shape and then raises while the
native Ipopt electrolyte LLE route builder is pending. ``solver_backend="ipopt"``
is reserved for that native constrained-NLP adapter and also raises until the
route is implemented.

The planned native adapter will expose material balance, charge balance, phase
amounts, thermodynamic objective terms, and derivative callbacks as formal NLP
blocks. The old Python IPOPT bridge has been removed.

Solver-selection policy
-----------------------

``solver_backend="auto"`` remains conservative and does not switch to any
package-owned solve loop while the native adapter is pending. Treat electrolyte
LLE calls as route-gated until the native Ipopt implementation lands.

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
       ),
   )

   assert result.split_detected
   assert result.diagnostics["phase_distance"] > 1.0e-4
