API Reference
=============

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

Regression helpers
------------------

The regression API remains Python-facing, but ``fit_pure_neutral(...)`` now delegates to the native least-squares regression engine linked into ``epcsaft.epcsaft``. The shipped v1 scope is a nonassociating neutral pure-component fit of \(m\), \(s\), and \(e\) against liquid-density and pure-VLE records.

.. autoclass:: epcsaft.FitBounds
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

.. autofunction:: epcsaft.fit_pure_neutral

.. autofunction:: epcsaft.write_fit_result


