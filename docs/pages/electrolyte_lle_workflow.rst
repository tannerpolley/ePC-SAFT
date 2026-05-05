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
       options=epcsaft.EquilibriumOptions(max_iterations=180, tolerance=1.0e-8),
   )

   assert result.split_detected
   assert result.diagnostics["phase_distance"] > 1.0e-4
