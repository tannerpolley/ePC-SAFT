Getting Started
===============

Install
-------

For a released package, use the standard package-manager command once PyPI
publishing is enabled:

.. code-block:: powershell

   python -m pip install epcsaft

Until then, install the current release from GitHub:

.. code-block:: powershell

   python -m pip install "epcsaft @ git+https://github.com/tannerpolley/ePC-SAFT.git@v1.5.1"

For a local editable checkout:

.. code-block:: powershell

   git clone https://github.com/tannerpolley/ePC-SAFT.git
   cd ePC-SAFT
   python -m pip install -e .

Source and editable installs build a native C++ extension. They require Python
``>=3.9``, a C++ compiler, CMake, and Ninja or another CMake generator. See
:doc:`release_installation` for the full install matrix.

Verify the install
------------------

.. code-block:: python

   import epcsaft

   print(epcsaft.__version__)
   print(epcsaft.runtime_build_info())

Create a mixture
----------------

For a one-component example, pass a parameter dictionary directly:

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
   print(state.compressibility_factor())
   print(state.fugacity_coefficient())

Use pressure or density closure
-------------------------------

Every state uses exactly one closure variable:

- ``P`` asks the EOS to solve for density.
- ``rho`` evaluates the model directly at a supplied molar density.
- ``rho_guess`` is allowed only with ``P`` and seeds the pressure-density solve.

.. code-block:: python

   base = mixture.state(T=320.0, x=np.asarray([1.0]), P=101325.0)
   nearby = mixture.state(
       T=321.0,
       x=np.asarray([1.0]),
       P=101325.0,
       rho_guess=base.density(),
   )

   density_state = mixture.state(T=320.0, x=np.asarray([1.0]), rho=base.density())
   audit = mixture.check_density(
       T=320.0,
       x=np.asarray([1.0]),
       P=101325.0,
       rho=density_state.density(),
   )
   print(audit["within_tolerance"], audit["pressure_residual"])

Create a parameter folder
-------------------------

For real systems, keep your parameter data in a folder you control:

.. code-block:: python

   from epcsaft import create_parameter_template

   template_root = create_parameter_template(
       location=r"C:\path\to\my_epcsaft_data",
       folder_name="water_salt_case",
       species=["H2O", "Na+", "Cl-"],
   )

Fill in the generated files, then load the folder:

.. code-block:: python

   import numpy as np
   from epcsaft import ePCSAFTMixture

   species = ["H2O", "Na+", "Cl-"]
   x = np.asarray([0.9998, 1e-4, 1e-4])
   mixture = ePCSAFTMixture.from_dataset(template_root, species, x, 298.15)

Next steps
----------

- :doc:`user_parameter_templates` explains the parameter-folder layout.
- :doc:`user_options` lists supported model options.
- :doc:`package_guide` gives task-based examples.
- :doc:`equilibrium_cookbook` shows phase-equilibrium and speciation workflows.
- :doc:`api_reference` lists the public API.
