=======
PC-SAFT
=======

.. image:: https://badge.fury.io/py/pcsaft.svg
    :target: https://badge.fury.io/py/pcsaft
.. image:: https://img.shields.io/badge/License-GPLv3-blue.svg
    :target: /LICENSE
.. image:: https://img.shields.io/pypi/dm/pcsaft
    :alt: PyPI - Downloads
    :target: https://pypistats.org/packages/pcsaft
.. image:: https://readthedocs.org/projects/pcsaft/badge/?version=latest
   :target: http://pcsaft.readthedocs.io/?badge=latest

Introduction
------------
``pcsaft`` is a Cython/C++ implementation of the PC-SAFT equation of state with dipole, association, and electrolyte terms.

Install
-------

.. code-block:: bash

   pip install .

For editable development:

.. code-block:: bash

   pip install -e . --no-build-isolation

Example
-------

.. code-block:: python

  import numpy as np
  from pcsaft import PCSAFTMixture

  mixture = PCSAFTMixture.from_params(
      {"m": np.asarray([2.8149]), "s": np.asarray([3.7169]), "e": np.asarray([285.69])},
      species=["Toluene"],
  )
  state = mixture.state(T=320.0, x=np.asarray([1.0]), P=101325.0)
  print(state.density())

Package Layout
--------------

- ``src/pcsaft/``: installable runtime package, Cython/C++ sources, and packaged parameter datasets
- ``data/``: datasets and figures that are not required at runtime
- ``docs/``: user documentation and reference material

Public API
----------

The main entry points are ``PCSAFTMixture``, ``PCSAFTState``, and the structured result objects returned by solver-style methods.

Author
------

Tanner Polley

License
-------

GNU General Public License v3.0
