Getting Started
===============

Install and build
-----------------

Install the package with pip:

.. code-block:: bash

   pip install .

For editable development:

.. code-block:: bash

   pip install -e . --no-build-isolation

Minimal example
---------------

.. code-block:: python

   import numpy as np
   from pcsaft import PCSAFTMixture

   mixture = PCSAFTMixture.from_dataset(
       "2012_Held",
       ["Na+", "Cl-", "H2O"],
       np.asarray([1e-4, 1e-4, 0.9998]),
       298.15,
   )
   state = mixture.state(T=298.15, x=np.asarray([1e-4, 1e-4, 0.9998]), P=1.0e5)

   rho = state.density()
   z = state.Z()
   act = state.actcoeff(species=["Na+", "Cl-", "H2O"])

Key objects
-----------

- Use `PCSAFTMixture` to load parameters and create states.
- Use `PCSAFTState` for pressure, density, residual properties, fugacity, and activity-coefficient calls.
- Use `ActivityCoeffResult` to inspect ion, solvent, and mean-ionic data in structured form.
