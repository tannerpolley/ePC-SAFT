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
   uv run python scripts/build_epcsaft.py
   uv run python scripts/codex_doctor.py
   uv run python run_pytest.py --confidence -q

This is the default new-agent validation sequence: sync, normal native build,
doctor, then the confidence suite.

For tests, use:

.. code-block:: bash

   uv run python run_pytest.py tests/test_runtime.py -q

Direct pytest also works:

.. code-block:: bash

   uv run python -m pytest tests/test_runtime.py -q

For Codex and Windows work, prefer ``uv run python run_pytest.py ...`` because the wrapper manages pytest temporary directories more predictably. Set ``EPCSAFT_PYTEST_TEMP_ROOT`` when you want the wrapper to use an opt-in external pytest temp root instead of its default repo-local generated temp area.

``--generic`` is the faster validation slice. ``--confidence`` is the default runtime-confidence check because it also runs the native runtime contract tests.

For native/equation work, use the native/equation debugging guide after the normal build and confidence sequence.

Use ``uv run python scripts/build_epcsaft.py --clean`` only as a repair step for stale CMake state or stale/locked ``_core`` artifacts. If a ``_core*.pyd`` is locked, stop the importing Python/test/IDE process before running the clean repair.

``CMakePresets.json`` is optional Windows MinGW convenience for IDEs and manual CMake use. The canonical local native build remains ``uv run python scripts/build_epcsaft.py``.

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


