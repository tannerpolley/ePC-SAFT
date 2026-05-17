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

Neutral ``bubble_p(...)`` and ``dew_p(...)`` use native Ipopt fixed-temperature
route builders when Ipopt is compiled. ``bubble_t(...)`` and ``dew_t(...)`` use
native Ipopt fixed-pressure temperature route builders on the same dependency
boundary.

``ePCSAFTMixture.equilibrium(kind="electrolyte_bubble_pressure", ...)`` and
``solve_reactive_electrolyte_bubble(...)`` use the native Ipopt fixed-liquid
electrolyte bubble route when Ipopt is compiled and the scoped staged inputs are
supported. They do not expose public pressure or vapor-composition seed knobs;
the route builder owns the canonical initial point.

``evaluate_fugacity_coefficients_batch(...)`` is the intended lightweight helper
for downstream-owned repeated property loops. It reuses the previous row's
resolved density as ``rho_guess`` when pressure closure is requested and
``continuation="auto"`` is selected.

``EquilibriumOptions`` keeps ``timeout_seconds`` as the public wall-clock option
for the native Ipopt equilibrium routes. Old electrolyte LLE seed and density
budget controls are no longer accepted by the public facade; callers should
catch the typed ``InputError``/``SolutionError`` diagnostics emitted by the
single native route attempt.

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

.. autofunction:: epcsaft.load_regression_records

.. autofunction:: epcsaft.validate_regression_provenance

.. autofunction:: epcsaft.fit_pure_parameters

.. autofunction:: epcsaft.fit_binary_parameters

.. autofunction:: epcsaft.fit_liquid_electrolyte_parameters

.. autofunction:: epcsaft.fit_pure_neutral

.. autofunction:: epcsaft.fit_pure_ion

.. autofunction:: epcsaft.fit_binary_pair

.. autoclass:: epcsaft.ReactiveElectrolyteBatch
   :members:
   :undoc-members:
   :no-index:

.. autoclass:: epcsaft.ReactiveElectrolyteBatchOptions
   :members:
   :undoc-members:
   :no-index:

.. autoclass:: epcsaft.ReactiveElectrolyteRow
   :members:
   :undoc-members:
   :no-index:

.. autoclass:: epcsaft.ReactiveElectrolyteRegressionContext
   :members:
   :undoc-members:
   :no-index:

.. autoclass:: epcsaft.ReactiveRegressionObjective
   :members:
   :undoc-members:
   :no-index:

.. autoclass:: epcsaft.ReactiveRegressionObjectiveResult
   :members:
   :undoc-members:
   :no-index:

.. autofunction:: epcsaft.evaluate_reactive_regression_objective

.. autofunction:: epcsaft.write_fit_result


