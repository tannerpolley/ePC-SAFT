Parameter Regression
====================

Phase 1 of the regression workflow is intentionally narrow:

- supported workflow: ``fit_pure_neutral(...)``
- supported targets: ``m``, ``s``, and ``e``
- supported property basis: liquid density and vapor pressure for one neutral component
- implementation owner: native least-squares workflow inside the main ``epcsaft.epcsaft`` extension
- current model restriction: nonassociating neutral pure-component fits only

Ion and binary regression are deferred until a later phase.

The public Python call shape stays the same, but the optimization loop itself is no longer Python/SciPy-owned. The current public workflow uses deterministic transformed-space starts and solves the supported pure-neutral regression surface with the native least-squares backend only.
Build prerequisite
------------------

There is no IPOPT prerequisite for the current package build. The supported developer path is the active Conda environment with the normal package build dependencies only.

The public Python call shape stays the same, but the optimization loop itself is no longer Python/SciPy-owned. The current public workflow is fully native and least-squares-only.

Current data basis
------------------

The current neutral-reference basis lives under ``data/pure_component/`` in this source checkout.

- ``hydrocarbon_basis_workbook_reference.csv`` stores workbook reference targets for ``Methane``, ``Ethane``, and ``Propane``
- the ``*_nist_saturation.csv`` files store saturation pressure and saturated liquid density values used as a clean starting dataset basis

Those files are useful for reproducing the workbook targets and for iterating on the neutral-only regression implementation.

Create a dataset folder
-----------------------

Start from a user-owned template folder.

.. code-block:: python

   from epcsaft import create_parameter_template

   template_root = create_parameter_template(
       location=r"C:\Users\Tanner\Documents\my_epcsaft_data",
       folder_name="hydrocarbon_case",
       species=["Methane"],
   )

The template gives you a user-owned ``pure/`` CSV that ``write_fit_result(...)`` can update after a fit.

Expected record schema
----------------------

``fit_pure_neutral(...)`` accepts flat records. The practical phase-1 schema is:

- ``T``: temperature in kelvin
- ``P``: vapor pressure in pascals
- ``rho``: liquid molar density at the same saturation state, in ``mol/m^3``
- or ``rho_kg_m3`` / ``rho_sat_liq_kg_m3``: liquid mass density in ``kg/m^3``; the regression helper converts this to molar density using ``MW``
- optional ``phase``: defaults to ``liq``

The records can come from your own CSV files, spreadsheet exports, or in-memory data structures. The regression helper does not require any package-owned data path.

Example:

.. code-block:: python

   records = [
       {"T": 100.0, "P": 34375.892, "rho_sat_liq_kg_m3": 438.88524, "phase": "liq"},
       {"T": 110.0, "P": 88130.038, "rho_sat_liq_kg_m3": 424.77725, "phase": "liq"},
   ]

Fit ``m``, ``s``, and ``e``
---------------------------

.. code-block:: python

   from epcsaft import fit_pure_neutral

   result = fit_pure_neutral(
       records,
       "Methane",
       assoc_scheme="",
       fixed_parameters={
           "MW": 0.0160428,
           "z": 0.0,
           "e_assoc": 0.0,
           "vol_a": 0.0,
           "dielc": 8.0,
           "d_born": 0.0,
           "f_solv": 1.0,
       },
       initial_guess={"m": 1.1, "s": 3.6, "e": 145.0},
   )

Phase-1 rules:

- only ``m``, ``s``, and ``e`` are fitted
- other pure-component fields are treated as fixed metadata
- the objective and exact first derivatives are evaluated natively
- the fit is read-only until you explicitly call ``write_fit_result(...)``

Inspect the result
------------------

.. code-block:: python

   print(result.success)
   print(result.backend)
   print(result.fitted_values)
   print(result.metrics_by_term)

Write the fitted values back
----------------------------

``write_fit_result(...)`` updates only the target component row in the user-owned dataset folder.

.. code-block:: python

   from epcsaft import write_fit_result

   written_paths = write_fit_result(result, template_root, overwrite=False)
   print(written_paths)

With ``overwrite=False``, blank template cells can be filled but existing values are protected.
