Equilibrium Architecture
========================

The recommended direct API remains method based:

* ``mixture.flash_tp(...)``
* ``mixture.stability_tp(...)``
* ``mixture.lle_tp(...)``
* ``mixture.electrolyte_lle_tp(...)``
* ``mixture.chemical_equilibrium(...)``

Neutral ``mixture.bubble_p(...)``, ``mixture.bubble_t(...)``,
``mixture.dew_p(...)``, and ``mixture.dew_t(...)`` remain public API names and
are implemented native Ipopt bubble/dew routes when Ipopt is compiled.
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

Request Normalization
---------------------

``mixture.equilibrium(kind=...)`` is a compatibility facade over the same typed
problem objects. The facade normalizes non-reactive string requests into one
typed problem plus route diagnostics before solving, so ``kind="auto"``,
explicit neutral routes, electrolyte LLE, electrolyte stability, bubble/dew,
and fixed-liquid electrolyte bubble pressure share the same result-stamping
contract as ``mixture.solve_equilibrium(problem)``.

Reactive convenience routes remain explicit. ``chemical_equilibrium``,
``reactive_staged_equilibrium``, ``reactive_lle``,
``reactive_electrolyte_lle``, ``reactive_stability``, and
``reactive_electrolyte_bubble_pressure`` are selected before non-reactive
normalization so their specialized option checks and native route boundaries do
not become implicit fallback behavior.

Solver Selection
----------------

``solver_backend="auto"`` is conservative: it may use a validated native
default only where one exists, and otherwise raises at the route boundary.
``ipopt`` is the explicit native constrained-NLP backend for production
equilibrium routes as they land. Current implemented native Ipopt capability
entries cover homogeneous ideal, activity, and concentration reactive
speciation, neutral TP flash, neutral LLE, neutral stability, fixed-temperature
and fixed-pressure neutral bubble/dew routes, electrolyte LLE, and fixed-liquid
electrolyte bubble pressure when Ipopt is compiled. Reactive phase-equilibrium
solves remain route-builder work. Use
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
