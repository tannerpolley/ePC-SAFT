Equilibrium Architecture
========================

The recommended direct API remains method based:

* ``mixture.flash_tp(...)``
* ``mixture.stability_tp(...)``
* ``mixture.lle_tp(...)``
* ``mixture.electrolyte_lle_tp(...)``
* ``mixture.chemical_equilibrium(...)``

Neutral ``mixture.bubble_p(...)`` and ``mixture.dew_p(...)`` remain public API
names and are implemented native Ipopt pressure routes when Ipopt is compiled.
``mixture.electrolyte_bubble_p(...)`` follows the same native Ipopt dependency
boundary. Public routes that still lack a production native route builder fail
loudly at the route boundary and are not advertised as implemented
``capabilities()["equilibrium"]`` entries.

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
* ``ReactivePhaseEquilibriumProblem``
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
equilibrium routes as they land. Current implemented native Ipopt capability
entries cover homogeneous ideal reactive speciation, neutral TP/LLE, neutral
bubble/dew pressure, electrolyte LLE, and fixed-liquid electrolyte bubble
pressure when Ipopt is compiled. Stability, temperature bubble/dew, and
reactive phase-equilibrium solves remain route-builder work. Use
``epcsaft.capabilities()`` to check implemented solver paths in the current
install.

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
state density into the next pressure-closed state as ``rho_guess``. Use direct ``rho=...`` only when density is the closure variable,
not when exact pressure closure is required.
