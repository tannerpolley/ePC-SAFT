API Reference
=============

Package boundary guidance lives in :doc:`package_architecture`. The public API
continues to support top-level imports while subsystem modules provide clearer
navigation for new equilibrium and regression code.

Mixture model
-------------

.. autoclass:: epcsaft.ePCSAFTMixture
   :members:
   :undoc-members:
   :no-index:

State model
-----------

.. autoclass:: epcsaft.ePCSAFTState
   :members:
   :undoc-members:
   :no-index:

The ``epcsaft.ePCSAFTState.composition_derivative_residual_helmholtz()`` method returns the structured composition-derivative breakdown used by ``residual_chemical_potential()``.
The ``epcsaft.ePCSAFTState.activity_coefficient(...)`` method returns component activity coefficients by default, and mean-ionic values when called with ``mean_ionic_form=True``.
The ``epcsaft.ePCSAFTState.density()`` method returns molar density by default and also accepts ``units="mass"`` for ``kg/m^3`` when molecular weights are available. Use ``molar_density()`` and ``mass_density()`` when you want the unit choice to be explicit in code.
``ePCSAFTMixture.state(...)`` accepts exactly one of ``P`` or ``rho``. With ``P``, optional ``rho_guess`` seeds the native pressure-density solve without changing the pressure-closure contract. With ``rho``, the state is evaluated directly at the supplied molar density.
``rho_seed`` is accepted as an alias for ``rho_guess`` in repeated-state workflows.
``ePCSAFTMixture.check_density(...)`` returns pressure-residual diagnostics for an externally supplied density and target pressure.
The ``epcsaft.ePCSAFTState.method_aliases()`` method returns the canonical state-method abbreviation map, and those aliases can be called directly on the state object. Example: ``ares()`` maps to ``residual_helmholtz()``.
``fugacity_coefficient()`` now defaults to the natural-log form. Pass ``natural_log=False`` to get the coefficient form instead, and use ``return_contribution_terms=True`` when you need the structured contribution payload.
The primitive contribution methods ``compressibility_factor()``, ``residual_helmholtz()``, ``temperature_derivative_residual_helmholtz()``, ``residual_chemical_potential()``, and ``fugacity_coefficient()`` accept ``return_contribution_terms=True`` to return a structured payload with ``total`` and per-contribution ``terms``. ``composition_derivative_residual_helmholtz()`` already returns its structured payload directly.
For ``fugacity_coefficient(..., return_contribution_terms=True)``, the returned ``terms`` stay in natural-log form and the structured payload also includes ``terms_total_natural_log``.

Activity results
----------------

.. autoclass:: epcsaft.ActivityCoefficientResult
   :members:
   :undoc-members:
   :no-index:

Runtime metadata
----------------

.. autofunction:: epcsaft.runtime_build_info

.. autofunction:: epcsaft.capabilities

Equilibrium and speciation
--------------------------

For runnable examples and solver-selection guidance, start with
:doc:`equilibrium_cookbook`.

Explicit mixture methods are preferred for new code: ``flash_tp(...)``,
``lle_tp(...)``, ``stability_tp(...)``, ``electrolyte_lle_tp(...)``,
``electrolyte_stability_tp(...)``, ``electrolyte_bubble_p(...)``,
``chemical_equilibrium(...)``, and ``reactive_staged_equilibrium(...)``.
The string-dispatched ``ePCSAFTMixture.equilibrium(kind=...)`` API remains
supported and routes through the explicit methods.

Neutral ``bubble_p(...)``, ``bubble_t(...)``, ``dew_p(...)``, and ``dew_t(...)``
remain declared API names, but their previous Python scalar route has been
removed. They raise ``InputError`` until native Ipopt route builders replace
that path.

``ePCSAFTMixture.equilibrium(kind="electrolyte_bubble_pressure", ...)`` uses the
native backend for fixed-liquid electrolyte bubble-pressure solves. Its current
scope keeps ions liquid-only and permits neutral vapor species.
``solve_reactive_electrolyte_bubble(...)`` and the matching sweep helper first
run native chemical speciation, then call the same native fixed-liquid
electrolyte bubble-pressure workflow.
Use ``ReactiveElectrolyteBubbleOptions(error_mode="result")`` for diagnostic
sweeps that should return structured nonconverged bubble-stage failures instead
of raising immediately.

``evaluate_fugacity_coefficients_batch(...)`` is the intended lightweight helper
for downstream-owned repeated property loops. It reuses the previous row's
resolved density as ``rho_seed`` when pressure closure is requested and
``continuation="auto"`` is selected.

``EquilibriumOptions`` also exposes optional diagnostic work budgets for
electrolyte LLE: ``timeout_seconds``, ``max_seed_attempts``,
``max_density_failures``, and ``max_total_objective_evaluations``. Defaults are
``None`` and preserve the robust solver path. When a configured budget is
exhausted, strict mode raises ``SolutionError`` with
``diagnostics["acceptance_gate"] == "predictive_budget_exhausted"`` and budget
fields so downstream sweeps can record a structured nonconvergence instead of
using only outer process timeouts.

.. autoclass:: epcsaft.EquilibriumOptions
   :members:
   :undoc-members:
   :no-index:

.. autoclass:: epcsaft.ElectrolyteBubbleOptions
   :members:
   :undoc-members:
   :no-index:

.. autoclass:: epcsaft.ElectrolyteBubbleResult
   :members:
   :undoc-members:
   :no-index:

.. autoclass:: epcsaft.ReactionDefinition
   :members:
   :undoc-members:
   :no-index:

.. autoclass:: epcsaft.ReactiveSpeciationOptions
   :members:
   :undoc-members:
   :no-index:

.. autoclass:: epcsaft.ReactiveSpeciationResult
   :members:
   :undoc-members:
   :no-index:

.. autofunction:: epcsaft.solve_reactive_speciation

.. autofunction:: epcsaft.solve_reactive_speciation_sweep

.. autofunction:: epcsaft.evaluate_fugacity_coefficients

.. autofunction:: epcsaft.evaluate_fugacity_coefficients_batch

.. autofunction:: epcsaft.validate_dataset_bundle

.. autoclass:: epcsaft.ReactiveElectrolyteBubbleOptions
   :members:
   :undoc-members:
   :no-index:

.. autoclass:: epcsaft.ReactiveElectrolyteBubbleResult
   :members:
   :undoc-members:
   :no-index:

.. autofunction:: epcsaft.solve_reactive_electrolyte_bubble

.. autofunction:: epcsaft.solve_reactive_electrolyte_bubble_sweep

Regression helpers
------------------

The regression API remains Python-facing. The easy public entry points
``fit_pure_parameters(...)``, ``fit_binary_parameters(...)``, and
``fit_liquid_electrolyte_parameters(...)`` provide stable problem contracts for
pure-component, binary, and liquid-electrolyte parameter fitting. The lower-level
native-backed helpers remain available for code that already targets them.

.. autoclass:: epcsaft.FitBounds
   :members:
   :undoc-members:
   :no-index:

.. autoclass:: epcsaft.FitParameter
   :members:
   :undoc-members:
   :no-index:

.. autoclass:: epcsaft.BinaryInteraction
   :members:
   :undoc-members:
   :no-index:

.. autoclass:: epcsaft.RelativePermittivityResidual
   :members:
   :undoc-members:
   :no-index:

.. autoclass:: epcsaft.FitTerm
   :members:
   :undoc-members:
   :no-index:

.. autoclass:: epcsaft.FitProblem
   :members:
   :undoc-members:
   :no-index:

.. autoclass:: epcsaft.FitResult
   :members:
   :undoc-members:
   :no-index:

.. autoclass:: epcsaft.ReactiveElectrolyteRegressionResult
   :members:
   :undoc-members:
   :no-index:

.. autofunction:: epcsaft.load_regression_records

.. autofunction:: epcsaft.validate_regression_provenance

.. autofunction:: epcsaft.fit_pure_parameters

.. autofunction:: epcsaft.fit_binary_parameters

.. autofunction:: epcsaft.fit_liquid_electrolyte_parameters

.. autofunction:: epcsaft.fit_pure_neutral

.. autofunction:: epcsaft.fit_pure_ion

.. autofunction:: epcsaft.fit_binary_pair

.. autofunction:: epcsaft.evaluate_reactive_electrolyte_bubble_residuals

.. autofunction:: epcsaft.write_fit_result


