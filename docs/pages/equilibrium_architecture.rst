Equilibrium Architecture
========================

The recommended direct API remains method based:

* ``mixture.flash_tp(...)``
* ``mixture.stability_tp(...)``
* ``mixture.lle_tp(...)``
* ``mixture.electrolyte_lle_tp(...)``
* ``mixture.chemical_equilibrium(...)``

Neutral ``mixture.bubble_p(...)`` and ``mixture.dew_p(...)`` remain public API
names, but they currently fail loudly until the native Ipopt route builders own
those solves.
``mixture.electrolyte_bubble_p(...)`` and
``mixture.reactive_electrolyte_bubble_p(...)`` follow the same route-pending
policy until native Ipopt electrolyte bubble builders land.

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

``solver_backend="auto"`` is conservative: it may use a validated native
default only where one exists, and otherwise raises at the route boundary.
``ipopt`` is the explicit native constrained-NLP backend for production
equilibrium routes as they land. The current public Ipopt route is homogeneous
ideal reactive speciation; broader multiphase, electrolyte, and EOS-coupled
equilibrium routes remain route-builder work. Use ``epcsaft.capabilities()`` to
check which optional solver paths are available in the current install.

The convex Gibbs formulation is limited to homogeneous ideal reaction or
speciation subkernels and validation tests. Full ePC-SAFT multiphase,
electrolyte, density-coupled, or association-coupled equilibrium should be
treated as a thermodynamic constrained NLP, not as a globally convex problem.
Production equilibrium routes require exact analytic or CppAD Jacobians. Native
Ceres owns package regression solves, while CppAD and implicit sensitivities
provide derivative payloads where the route is validated.

Repeated State Work
-------------------

For many property calls, keep the loop downstream-owned. Feed each successful
state density into the next pressure-closed state as ``rho_guess`` or
``rho_seed``. Use direct ``rho=...`` only when density is the closure variable,
not when exact pressure closure is required.
