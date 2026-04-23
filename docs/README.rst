Overview
========

.. image:: https://badge.fury.io/py/epcsaft.svg
    :target: https://badge.fury.io/py/epcsaft
.. image:: https://img.shields.io/badge/License-GPLv3-blue.svg
    :target: /LICENSE
.. image:: https://img.shields.io/pypi/dm/epcsaft
    :alt: PyPI - Downloads
    :target: https://pypistats.org/packages/epcsaft
.. image:: https://readthedocs.org/projects/epcsaft/badge/?version=latest
   :target: http://epcsaft.readthedocs.io/?badge=latest

``epcsaft`` is a Cython/C++ implementation of the ePC-SAFT equation of state with dipole, association, and electrolyte terms.

Where to start
--------------

Most users should install from PyPI:

.. code-block:: bash

   pip install epcsaft

``epcsaft`` includes a compiled Cython/C++ extension. If a wheel is available for your Python version and platform, pip installs it automatically. If a wheel is not available, pip falls back to a source build, which requires a working native build toolchain.

For a source checkout of this repository:

.. code-block:: bash

   pip install .

For editable development from this source tree:

.. code-block:: bash

   python scripts/install_dev.py
   python scripts/build_epcsaft.py

``install_dev.py`` creates or repairs the editable install. ``build_epcsaft.py`` is the fast native/Cython iteration command and rebuilds the in-place extension only when tracked build inputs are stale. Use ``python run_pytest.py`` for tests and ``python scripts/build_dist.py`` for a source distribution, wheel, and wheel smoke check. If you want to call pip directly for the editable install, use ``pip install -e . --no-build-isolation --config-settings editable_mode=compat``.

Then choose the path that matches what you want to do:

- :doc:`getting_started` for the quickest first run
- :doc:`user_parameter_templates` to build your own parameter folder
- :doc:`parameter_regression` for the phase-1 neutral ``m/s/e`` regression workflow
- :doc:`user_options` to see every supported ``user_options.json`` setting
- :doc:`package_guide` for a task-based guide to the public API
- :doc:`api_reference` for the full method list

Simple example
--------------

.. code-block:: python

  import numpy as np
  from epcsaft import ePCSAFTMixture

  mixture = ePCSAFTMixture.from_params(
      {"m": np.asarray([2.8149]), "s": np.asarray([3.7169]), "e": np.asarray([285.69])},
      species=["Toluene"],
  )
  state = mixture.state(T=320.0, x=np.asarray([1.0]), P=101325.0)
  # Pressure-based states solve and cache density during construction.
  print(state.density())               # mol/m^3 by default
  print(state.density(units="mass"))   # kg/m^3 when MW is available
  print(state.molar_density())         # explicit molar-density alias
  print(state.ares())                  # short alias for residual_helmholtz()
  print(state.ares(return_contribution_terms=True)["terms"]["hc"])

Project Layout
--------------

- ``src/epcsaft/``: installable runtime package and Cython/C++ sources
- ``data/epcsaft_parameters/``: source-checkout example parameter datasets for inspection, comparison, and tests
- ``data/``: other datasets and figures that are not required by the package
- ``docs/``: user documentation and reference material

Public API
----------

The main entry points are ``create_parameter_template``, ``ePCSAFTMixture``, ``ePCSAFTState``, ``fit_pure_neutral(...)``, and the structured result objects returned by solver-style methods.

If you want to build your own parameter folder, see :doc:`user_parameter_templates`.

Author
------

Tanner Polley

License
-------

GNU General Public License v3.0


