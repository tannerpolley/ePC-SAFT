Overview
========

``epcsaft`` is a Python package for electrolyte PC-SAFT thermodynamic
calculations. The public interface is Python, while the equation-of-state
runtime and package-owned equilibrium/regression kernels are implemented in
native C++ through ``pybind11``.

Current release: ``1.5.2``

What the package does
---------------------

Use ``epcsaft`` to build PC-SAFT/ePC-SAFT mixtures, evaluate thermodynamic
states, compute fugacity and activity coefficients, run supported
phase-equilibrium workflows, and fit supported parameter sets against tabular
data.

The main user objects are:

- ``ePCSAFTMixture``: stores species parameters and creates states.
- ``ePCSAFTState``: evaluates density, pressure, residual properties, fugacity
  coefficients, activity coefficients, and diagnostics.
- ``create_parameter_template(...)``: creates a user-owned parameter folder.
- ``fit_pure_neutral(...)`` and related helpers: fit supported parameter sets.
- ``capabilities()``: reports available runtime and solver paths.

Install
-------

The standard install command is:

.. code-block:: powershell

   python -m pip install epcsaft

With ``uv``:

.. code-block:: powershell

   uv add epcsaft

The current public release is also available from GitHub.

Install from the current GitHub release:

``https://github.com/tannerpolley/ePC-SAFT/releases/tag/v1.5.2``

If a wheel matching your platform is attached to the release, install it
directly:

.. code-block:: powershell

   python -m pip install C:\path\to\epcsaft-1.5.2-*.whl

To install from the tagged source:

.. code-block:: powershell

   python -m pip install "epcsaft @ git+https://github.com/tannerpolley/ePC-SAFT.git@v1.5.2"

With ``uv``:

.. code-block:: powershell

   uv add "epcsaft @ git+https://github.com/tannerpolley/ePC-SAFT.git@v1.5.2"

Source builds require Python ``>=3.9``, a C++ compiler, CMake, and Ninja or
another CMake generator. Python 3.13 is the current project smoke-test
baseline.

For a normal local source install:

.. code-block:: powershell

   git clone https://github.com/tannerpolley/ePC-SAFT.git
   cd ePC-SAFT
   python -m pip install .

For an editable install while changing Python files:

.. code-block:: powershell

   python -m pip install -e .

Editable installs use the same native build backend as wheel installs. Python
source changes are picked up from the checkout. If you change C++ sources,
pybind bindings, CMake files, or build metadata, rerun
``python -m pip install -e .``.

Verify the install
------------------

.. code-block:: python

   import epcsaft

   print(epcsaft.__version__)
   print(epcsaft.runtime_build_info())
   print(epcsaft.capabilities())

Quick example
-------------

.. code-block:: python

   import numpy as np
   from epcsaft import ePCSAFTMixture

   mixture = ePCSAFTMixture.from_params(
       {
           "m": np.asarray([2.8149]),
           "s": np.asarray([3.7169]),
           "e": np.asarray([285.69]),
       },
       species=["Toluene"],
   )

   state = mixture.state(T=320.0, x=np.asarray([1.0]), P=101325.0)

   print(state.density())
   print(state.pressure())
   print(state.compressibility_factor())
   print(state.fugacity_coefficient())

Pressure, density, and seeds
----------------------------

State construction uses exactly one closure variable:

- ``state(..., P=...)`` solves the EOS pressure-density closure.
- ``state(..., rho=...)`` evaluates properties at the supplied molar density.
- ``state(..., P=..., rho_guess=...)`` still solves exact pressure closure, but
  seeds the density solve with a previous good density.

.. code-block:: python

   base = mixture.state(T=320.0, x=np.asarray([1.0]), P=101325.0)
   next_state = mixture.state(
       T=321.0,
       x=np.asarray([1.0]),
       P=101325.0,
       rho_guess=base.density(),
   )

   audit = mixture.check_density(
       T=320.0,
       x=np.asarray([1.0]),
       P=101325.0,
       rho=base.density(),
   )
   print(audit["within_tolerance"], audit["pressure_residual"])

Parameter data
--------------

Most users should create and own their parameter folders:

.. code-block:: python

   from epcsaft import create_parameter_template

   template_root = create_parameter_template(
       location=r"C:\path\to\my_epcsaft_data",
       folder_name="water_salt_case",
       species=["H2O", "Na+", "Cl-"],
   )

After filling in the generated files, load them with
``ePCSAFTMixture.from_dataset(...)``. The source checkout contains
reference/example datasets under ``data/reference/epcsaft_parameters/`` for
comparison and validation, but users should not assume every installed wheel
contains those source-checkout reference folders.

Equilibrium and speciation
--------------------------

Use ``capabilities()`` and the cookbook before wiring a high-level equilibrium
workflow. The package includes native-backed paths for neutral phase
equilibrium, electrolyte LLE, fixed-liquid electrolyte bubble pressure,
reactive speciation, reactive staged equilibrium, and scoped reactive
electrolyte bubble pressure.

Important boundaries:

- Electrolyte bubble pressure is for fixed liquid composition with neutral
  vapor species; ions remain liquid-only.
- Reactive electrolyte bubble pressure is staged: native speciation first,
  then fixed-liquid electrolyte bubble pressure.
- IPOPT is optional, experimental, and explicit opt-in.
- Downstream case-study models should own their own data, balances, run
  matrices, and acceptance criteria.

Where to go next
----------------

- :doc:`release_installation` for release downloads and source-build notes.
- :doc:`getting_started` for the first local calculation.
- :doc:`user_parameter_templates` to build your own parameter folder.
- :doc:`user_options` for supported ``user_options.json`` settings.
- :doc:`equilibrium_cookbook` for equilibrium and speciation examples.
- :doc:`package_guide` for task-based API guidance.
- :doc:`api_reference` for the full API reference.

License
-------

GNU General Public License v3.0.
