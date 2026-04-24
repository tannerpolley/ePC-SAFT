Getting Started
===============

Install and build
-----------------

Most users install the published package from PyPI:

.. code-block:: bash

   pip install epcsaft

``epcsaft`` includes a compiled C++ extension exposed through pybind11. If a wheel is available for your Python version and platform, pip installs it automatically. If a wheel is not available, pip falls back to a source build, which requires a working native build toolchain.

For development from this source tree, use ``uv`` and the direct CMake build loop:

.. code-block:: bash

   uv sync --no-install-project
   uv run python scripts/build_epcsaft.py --clean
   uv run python scripts/codex_doctor.py

For tests, use:

.. code-block:: bash

   uv run python run_pytest.py tests/test_runtime.py -q

For package artifacts, use:

.. code-block:: bash

   uv build

On Windows, a compiled ``.pyd`` file may appear during local builds. That file is a build artifact and is normally produced by the installer or build backend, not copied into a project manually.

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


