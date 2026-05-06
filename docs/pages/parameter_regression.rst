Parameter Regression
====================

The package-owned regression helpers are record-driven. They accept flat in-memory records, tabular objects, or CSV files and return ``FitResult`` objects. Nothing is written to a parameter folder until ``write_fit_result(...)`` is called.

Supported workflows
-------------------

- ``fit_pure_neutral(...)`` fits nonassociating neutral pure-component ``m``, ``s``, and ``e`` against density and vapor-pressure records with the native least-squares backend.
- ``fit_pure_ion(...)`` fits ion ``s`` and ``e`` by default, and can fit ``d_born`` when requested, with the native least-squares backend.
- ``fit_binary_pair(...)`` fits constant binary ``k_ij``, ``l_ij``, and ``k_hb_ij`` values from VLE x/y records with the native least-squares backend.

Ion and binary V1 intentionally do not add dataset manifests or new regression-specific parameter namespaces. The helpers build runtime states from the existing dataset loader and caller-provided records.

Build prerequisite
------------------

There is no IPOPT prerequisite for the current package build. The supported developer path is the uv-managed environment plus the direct CMake/pybind11 native build:

.. code-block:: powershell

   uv sync --no-install-project
   uv run python scripts\build_epcsaft.py

Create a dataset folder
-----------------------

Start from a user-owned template folder.

.. code-block:: python

   from epcsaft import create_parameter_template

   template_root = create_parameter_template(
       location=r"C:\Users\Tanner\Documents\my_epcsaft_data",
       folder_name="regression_case",
       species=["H2O", "Na+", "Cl-"],
   )

The template gives you user-owned ``pure/`` and ``mixed/binary_interaction/`` CSV files that ``write_fit_result(...)`` can update after a fit.

Neutral records
---------------

``fit_pure_neutral(...)`` accepts flat records with:

- ``T``: temperature in kelvin
- ``P``: vapor pressure in pascals
- ``rho``: liquid molar density in ``mol/m^3``
- or ``rho_kg_m3`` / ``rho_sat_liq_kg_m3``: liquid mass density in ``kg/m^3``
- optional ``phase``: defaults to ``liq``

Example:

.. code-block:: python

   from epcsaft import fit_pure_neutral

   records = [
       {"T": 100.0, "P": 34375.892, "rho_sat_liq_kg_m3": 438.88524},
       {"T": 110.0, "P": 88130.038, "rho_sat_liq_kg_m3": 424.77725},
   ]

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

Ion records
-----------

``fit_pure_ion(...)`` records require ``T`` and ``P`` plus one composition basis:

- full mole-fraction columns such as ``x_H2O``, ``x_Na+``, and ``x_Cl-``
- or ``molality`` with explicit ``species=[...]`` and ``solvent=...``

Each ion regression problem must include at least one of:

- ``osmotic_coefficient`` or ``osmotic``
- ``mean_ionic_activity``, ``mean_ionic_activity_coefficient``, or ``miac``

Density is optional and is included when ``rho`` or a supported mass-density column is present.

Example:

.. code-block:: python

   from epcsaft import fit_pure_ion

   records = [
       {"T": 298.15, "P": 101325.0, "molality": 0.1, "osmotic_coefficient": 0.933},
       {"T": 298.15, "P": 101325.0, "molality": 0.2, "mean_ionic_activity": 0.735},
   ]

   result = fit_pure_ion(
       records,
       "Na+",
       dataset="2026_Khudaida",
       species=["H2O", "Na+", "Cl-"],
       solvent="H2O",
       initial_guess={"s": 2.6, "e": 210.0},
       bounds={"s": (2.4, 3.2), "e": (150.0, 300.0)},
   )

Advanced electrolyte options pass through to every runtime state build:

.. code-block:: python

   user_options = {
       "elec_model": {
           "rel_perm": {"rule": "empirical", "differential_mode": "auto"},
           "born_model": {
               "d_Born_mode": 3,
               "solvation_shell_model": True,
               "dielectric_saturation": True,
               "mu_born_model": {"differential_mode": "auto", "comp_dep_delta_d": True},
           },
       }
   }

   result = fit_pure_ion(
       records,
       "Na+",
       dataset="2026_Khudaida",
       species=["H2O", "Na+", "Cl-"],
       solvent="H2O",
       fit_targets=("d_born",),
       initial_guess={"d_born": 3.2},
       user_options=user_options,
   )

Binary VLE records
------------------

``fit_binary_pair(...)`` V1 supports VLE x/y records only. Records require:

- ``T``: temperature in kelvin
- ``P``: pressure in pascals
- liquid mole-fraction columns such as ``x_H2O`` and ``x_Ethanol``
- vapor mole-fraction columns such as ``y_H2O`` and ``y_Ethanol``

The V1 targets are constant ``k_ij``, ``l_ij``, and ``k_hb_ij``. Linear temperature models and LLE fitting are future phases and raise ``InputError``.

Example:

.. code-block:: python

   from epcsaft import fit_binary_pair

   records = [
       {"T": 330.0, "P": 101325.0, "x_H2O": 0.7, "x_Ethanol": 0.3, "y_H2O": 0.5, "y_Ethanol": 0.5},
       {"T": 340.0, "P": 101325.0, "x_H2O": 0.6, "x_Ethanol": 0.4, "y_H2O": 0.4, "y_Ethanol": 0.6},
   ]

   result = fit_binary_pair(
       records,
       ("H2O", "Ethanol"),
       dataset="2026_Khudaida",
       initial_guess={"k_ij": -0.02},
       bounds={"k_ij": (-0.2, 0.2)},
   )

Inspect and write results
-------------------------

.. code-block:: python

   print(result.success)
   print(result.backend)
   print(result.jacobian_backend)
   print(result.hessian_backend)
   print(result.fitted_values)
   print(result.metrics_by_term)

   from epcsaft import write_fit_result

   written_paths = write_fit_result(result, template_root, overwrite=False)
   print(written_paths)

With ``overwrite=False``, blank template cells can be filled but existing values are protected. Pure ion fits update the target component row in ``pure/``. Binary fits update both symmetric cells in the relevant interaction matrix.

Derivative and Jacobian access
------------------------------

Normal ``FitResult`` payloads report compact derivative metadata:

- ``jacobian_available``
- ``jacobian_backend``
- ``jacobian_fallback_used``
- ``jacobian_fallback_reason``
- ``finite_difference_fallback_count``
- ``hessian_available``
- ``hessian_backend``
- ``hessian_fallback_used``
- ``hessian_fallback_reason``

Large matrices are exposed only through explicit derivative-evaluation helpers. Use ``evaluate_pure_neutral_derivatives(...)`` for the native pure-neutral objective. It returns residuals, gradient, ``jacobian_row_major``, ``jacobian_shape``, and Hessian skeleton fields. Pure-neutral Jacobians use the native autodiff path.

For lower-level generic native records, use ``evaluate_generic_regression_derivatives(...)``. It returns residuals plus ``jacobian_row_major`` and ``jacobian_shape``. The current generic-regression Jacobian backend is reported honestly as ``finite_difference`` until the residual state calls are scalar-templated for autodiff.

Derivative availability
-----------------------

.. list-table::
   :header-rows: 1

   * - Method
     - Current Jacobian access
     - Hessian status
   * - Runtime ``dadt()``, ``dadx()``, ``z(return_contribution_terms=True)``, ``mures(return_contribution_terms=True)``
     - Analytical where available, autodiff where implemented, finite-difference fallback with metadata
     - Not exposed
   * - Pure-neutral regression
     - Native autodiff Jacobian through ``evaluate_pure_neutral_derivatives(...)``
     - Skeleton metadata only
   * - Generic ion/binary regression
     - Explicit finite-difference Jacobian through ``evaluate_generic_regression_derivatives(...)`` until generic autodiff coverage is implemented
     - Skeleton metadata only
   * - Neutral LLE
     - Native finite-difference Newton Jacobian with fallback diagnostics; autodiff residual boundary is planned
     - Skeleton metadata only
   * - Chemical equilibrium / reactive speciation
     - Native finite-difference Newton Jacobian with fallback diagnostics; autodiff residual boundary is planned
     - Skeleton metadata only

The Hessian fields are deliberately a contract skeleton for future IPOPT-compatible optimizer integration. They do not mean IPOPT support or second-derivative evaluation is implemented.
