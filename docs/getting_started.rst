Getting Started
===============

Install and build
-----------------

Most users install the published package from PyPI:

.. code-block:: bash

   pip install epcsaft

``epcsaft`` includes a compiled Cython/C++ extension. If a wheel is available for your Python version and platform, pip installs it automatically. If a wheel is not available, pip falls back to a source build, which requires a working native build toolchain.

To install the package from a source checkout of this repository:

.. code-block:: bash

   pip install .

For editable development from this source tree:

.. code-block:: bash

   python scripts/build_epcsaft.py

On Windows, a compiled ``.pyd`` file may appear during local builds. That file is a build artifact and is normally produced by the installer or build backend, not copied into a project manually.

If you want to call pip directly, use:

.. code-block:: bash

   pip install -e . --config-settings editable_mode=compat

If you change the Cython or C++ sources, rerun ``python scripts/build_epcsaft.py`` to refresh the editable install.

Create your own parameter folder
--------------------------------

If you want to build your own parameter folder, use ``create_parameter_template(...)``.

.. code-block:: python

   from epcsaft import create_parameter_template

   template_root = create_parameter_template(
       location=r"C:\Users\Tanner\Documents\my_epcsaft_data",
       folder_name="water_salt_case",
       species=["H2O", "Na+", "Cl-"],
   )

This creates the expected ``pure/``, ``mixed/``, and ``user_options.json`` files for you to fill in. After you add your parameters, pass ``template_root`` into ``ePCSAFTMixture.from_dataset(...)``.

If you want to fit missing entries instead of typing them by hand, the regression helpers can work directly with the same folder structure. See ``parameter_regression`` for the current phase-1 neutral ``m/s/e`` workflow.

If you are working from a checkout of this repository, you can also inspect the example parameter folders in ``data/epcsaft_parameters/``. They are there for comparison and reference while you build your own folders.

For a complete list of ``user_options.json`` settings, see ``user_options``.

Minimal example
---------------

.. code-block:: python

   import numpy as np
   from epcsaft import create_parameter_template, ePCSAFTMixture

   template_root = create_parameter_template(
       location=r"C:\Users\Tanner\Documents\my_epcsaft_data",
       folder_name="water_salt_case",
       species=["H2O", "Na+", "Cl-"],
   )
   mixture = ePCSAFTMixture.from_dataset(
       template_root,
       ["H2O", "Na+", "Cl-"],
       np.asarray([0.9998, 1e-4, 1e-4]),
       298.15,
   )
   state = mixture.state(T=298.15, x=np.asarray([0.9998, 1e-4, 1e-4]), P=1.0e5)
   # Pressure-based states solve and cache density during construction.

   rho_molar = state.density()
   rho_mass = state.density(units="mass")
   z = state.compressibility_factor()
   activity = state.activity_coefficient(species=["H2O", "Na+", "Cl-"])
   mean_ionic = state.activity_coefficient(
       species=["H2O", "Na+", "Cl-"],
       mean_ionic_form=True,
       basis="molality",
   )

Key objects
-----------

- Use `ePCSAFTMixture` to load parameters and create states.
- Use `ePCSAFTState` for pressure, density, residual properties, fugacity, and activity-coefficient calls.
- Use `activity_coefficient(...)` with `mean_ionic_form=True` when you want mean-ionic activity data on a chosen basis.


