Electrolyte LLE Workflow
========================

Hard electrolyte LLE cases should be handled as fixed-species phase-equilibrium
problems before chemical or speciation equilibrium is introduced.

Recommended workflow
--------------------

1. Run ``kind="electrolyte_stability"`` to confirm whether the feed is unstable.
2. If ``min_tpd`` is negative, call ``kind="electrolyte_lle"`` with the feed
   specification only. The native route builder owns the canonical initial
   point for the production solve.
3. For curves, pass independent feed/specification points. Ordered curve calls
   use each route builder's canonical initial point instead of caller-provided
   phase seeds.

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

Native Ipopt route
------------------

``EquilibriumOptions`` accepts ``solver_backend="auto" | "ipopt"``. The
default ``auto`` validates the public request shape and uses the native Ipopt
electrolyte LLE route when the extension was built with Ipopt. Local builds
without Ipopt raise a typed dependency error before any package-owned solve
loop can run.

The native adapter exposes material balance, charge balance, phase amounts,
thermodynamic objective terms, and derivative callbacks as formal NLP blocks.
The old Python IPOPT bridge has been removed.

Solver-selection policy
-----------------------

``solver_backend="auto"`` remains conservative and does not switch to any
package-owned solve loop. Treat electrolyte LLE calls as native Ipopt routes:
if the native route cannot run, the package raises instead of accepting a
caller-provided phase seed or alternate Python solve path.

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
   mix = ePCSAFTMixture.from_dataset("2024_Hubach", species, feed, 294.15)
   result = mix.equilibrium(
       kind="electrolyte_lle",
       T=294.15,
       P=1.013e5,
       z=feed,
       options=epcsaft.EquilibriumOptions(
           max_iterations=180,
           tolerance=1.0e-8,
           timeout_seconds=15.0,
       ),
   )

   assert result.split_detected
   assert result.diagnostics["phase_distance"] > 1.0e-4
