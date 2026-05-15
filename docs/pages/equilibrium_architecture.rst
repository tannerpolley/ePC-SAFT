Equilibrium Architecture
========================

The recommended direct API remains method based:

* ``mixture.flash_tp(...)``
* ``mixture.stability_tp(...)``
* ``mixture.bubble_p(...)`` and ``mixture.dew_p(...)``
* ``mixture.lle_tp(...)``
* ``mixture.electrolyte_lle_tp(...)``
* ``mixture.electrolyte_bubble_p(...)``
* ``mixture.chemical_equilibrium(...)``
* ``mixture.reactive_electrolyte_bubble_p(...)``

For agents or workflow builders that need a serializable problem object, use
``mixture.solve_equilibrium(problem)`` with one of:

* ``TPFlash``
* ``StabilityAnalysis``
* ``BubblePoint``
* ``DewPoint``
* ``LLEProblem``
* ``ElectrolyteLLEProblem``
* ``ElectrolyteBubblePoint``
* ``ReactiveSpeciationProblem``
* ``ReactiveElectrolyteBubbleProblem``

Example
-------

.. code-block:: python

   import epcsaft

   result = mixture.solve_equilibrium(
       epcsaft.TPFlash(T=220.0, P=1.0e5, z=[0.1, 0.3, 0.6])
   )

Solver Selection
----------------

``solver_backend="auto"`` uses the production native defaults. IPOPT is an
explicit opt-in residual-minimization bridge; it is not selected automatically
and is not a full constrained Gibbs/NLP formulation. Use
``epcsaft.capabilities()`` to check which optional solver paths are available
in the current install.

Repeated State Work
-------------------

For many property calls, keep the loop downstream-owned. Feed each successful
state density into the next pressure-closed state as ``rho_guess`` or
``rho_seed``. Use direct ``rho=...`` only when density is the closure variable,
not when exact pressure closure is required.
