Package Guide
=============

This guide groups the public package methods by task.

Constructing a model
--------------------

Use `ePCSAFTMixture` to load parameters and create states.

- `create_parameter_template(...)` makes a new user-owned dataset folder with blank files you can fill in.
- `fit_pure_neutral(...)` is the current phase-1 regression entrypoint for neutral-component `m`, `s`, and `e`.
- `write_fit_result(...)` writes a neutral regression result back into a user-owned dataset folder when you explicitly want to persist it.
- `ePCSAFTMixture.from_dataset(...)` loads a folder path you created yourself. If you are working from a checkout of this repository, it can also load the example folders under `data/reference/epcsaft_parameters/`.
- `ePCSAFTMixture.from_params(...)` builds a mixture from a resolved parameter dict.
- `ePCSAFTMixture.state(...)` creates a new state snapshot from temperature, composition, and either pressure or density.
  Pressure-based states solve and cache density during construction; density-based states still compute pressure on demand.

For the current regression basis, see the neutral reference files under `data/reference/pure_component/` in this source checkout.

Evaluating state properties
---------------------------

Use `ePCSAFTState` for point-property calculations.

- State methods use the `T`, `x`, phase, and cached density/pressure already bound to the state object. You do not pass `T`, `P`, or `rho` again when calling `residual_helmholtz()`, `compressibility_factor()`, or similar methods.
- `density()` returns molar density in `mol/m^3` by default.
- `density(units="mass")` returns mass density in `kg/m^3` when molecular weights are available.
- `molar_density()` and `mass_density()` are explicit aliases for the same two density bases.
- `pressure()` returns pressure in pascals.
- `compressibility_factor()` returns the compressibility factor.
- `residual_helmholtz()`, `temperature_derivative_residual_helmholtz()`, `residual_enthalpy()`, `residual_entropy()`, and `residual_gibbs()` return residual properties.
- `composition_derivative_residual_helmholtz()` returns the structured per-contribution composition-derivative breakdown used by `residual_chemical_potential()`.
- `fugacity_coefficient()` now defaults to the natural-log form. Pass `natural_log=False` if you need the coefficient form. Use `return_contribution_terms=True` on the same method when you want the structured contribution breakdown.
- `residual_chemical_potential()` returns residual chemical potentials.
- `method_aliases()` returns the canonical short-name map for state methods, and the aliases can be called directly. For example, `ares()` maps to `residual_helmholtz()` and `mures()` maps to `residual_chemical_potential()`.
- `return_contribution_terms=True` is available on the primitive contribution methods `compressibility_factor()`, `residual_helmholtz()`, `temperature_derivative_residual_helmholtz()`, `residual_chemical_potential()`, and `fugacity_coefficient()`. `composition_derivative_residual_helmholtz()` always returns its structured payload directly.
- For `fugacity_coefficient(..., return_contribution_terms=True)`, the contribution payload remains in natural-log form under `terms`, while `total` follows the chosen `natural_log` basis. The structured payload also includes `terms_total_natural_log`.

Activity-coefficient data
-------------------------

Use `activity_coefficient(...)` when you need electrolyte activity outputs.

- `activity_coefficient(...)` returns component activity coefficients by default.
- `activity_coefficient(..., mean_ionic_form=True, basis="mole")` and `activity_coefficient(..., mean_ionic_form=True, basis="molality")` return mean-ionic values on the requested basis.
- `solvation_free_energy()` returns ion solvation data.
- `osmotic_coefficient()` returns the osmotic coefficient.
- `relative_permittivity()` returns the dielectric diagnostics used by electrolyte terms.

Result objects
--------------

The package also exposes a lower-level activity result type if you want to inspect the values returned by the activity-coefficient calculations:

- `ActivityCoefficientResult`

Example
-------

.. code-block:: python

   import numpy as np
   from epcsaft import create_parameter_template, ePCSAFTMixture

   template_root = create_parameter_template(
       location=r"C:\Users\Tanner\Documents\my_epcsaft_data",
       folder_name="water_salt_case",
       species=["H2O", "Na+", "Cl-"],
   )
   mixture = ePCSAFTMixture.from_dataset(template_root, ["Na+", "Cl-", "H2O"], np.asarray([1e-4, 1e-4, 0.9998]), 298.15)
   state = mixture.state(T=298.15, x=np.asarray([1e-4, 1e-4, 0.9998]), P=1.0e5)
   # Pressure-based states solve and cache density during construction.

   rho_molar = state.density()
   rho_mass = state.density(units="mass")
   activity = state.activity_coefficient(species=["Na+", "Cl-", "H2O"])
   mean_ionic = state.activity_coefficient(
       species=["Na+", "Cl-", "H2O"],
       mean_ionic_form=True,
       basis="mole",
   )
   gamma_na_cl = mean_ionic["Na+Cl-"]


