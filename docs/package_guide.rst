Package Guide
=============

This guide groups the public package methods by task.

Constructing a model
--------------------

Use `PCSAFTMixture` to load parameters and create states.

- `PCSAFTMixture.from_dataset(...)` loads packaged dataset parameters.
- `PCSAFTMixture.from_params(...)` builds a mixture from a resolved parameter dict.
- `PCSAFTMixture.state(...)` creates a new state snapshot from temperature, composition, and either pressure or density.

Evaluating state properties
---------------------------

Use `PCSAFTState` for point-property calculations.

- `density()` and `pressure()` return the bound intensive variable.
- `Z()` returns the compressibility factor.
- `a_res()`, `dadt()`, `h_res()`, `s_res()`, and `g_res()` return residual properties.
- `lnfugcoef()` and `lnfugcoef_terms()` expose fugacity information.
- `mu_res()` and `gamma()` return residual chemical-potential style outputs.
- `cp(...)` estimates the heat capacity from an Aly-Lee fit.

Activity-coefficient data
-------------------------

Use `actcoeff(...)` when you need all electrolyte activity outputs together.

- `actcoeff(...)` bundles component values, ion values, mean-ionic values, and osmotic information.
- `miac(...)` and `miac_m(...)` return mean-ionic activities on the mole-fraction and molality bases.
- `gsolv(...)` returns ion solvation data.
- `osmoticC()` returns the osmotic coefficient.
- `dielectric_eval()` returns the dielectric diagnostics used by electrolyte terms.

Flash and vaporization
----------------------

- `flashTQ(...)` solves a temperature-quality flash.
- `flashPQ(...)` solves a pressure-quality flash.
- `Hvap(...)` returns the vaporization pressure result for the current state.

Result objects
--------------

The package returns structured result objects for solver-style methods:

- `ActivityCoeffResult`
- `FlashResult`
- `VaporizationResult`
- `PhaseResult`

Example
-------

.. code-block:: python

   import numpy as np
   from pcsaft import PCSAFTMixture

   mixture = PCSAFTMixture.from_dataset("2012_Held", ["Na+", "Cl-", "H2O"], np.asarray([1e-4, 1e-4, 0.9998]), 298.15)
   state = mixture.state(T=298.15, x=np.asarray([1e-4, 1e-4, 0.9998]), P=1.0e5)

   rho = state.density()
   act = state.actcoeff(species=["Na+", "Cl-", "H2O"])
   gamma_na_cl = act.mean_ionic_x()["Na+Cl-"]
