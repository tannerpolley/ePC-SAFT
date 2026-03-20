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

Repository roles
----------------
This repository currently serves three purposes:

- the installable ``pcsaft`` package
- the package-owned runtime parameter datasets used by ``pcsaft.parameters``
- the in-repo validation and paper-reproduction workspace under ``scripts/``

Introduction
------------
This package implements the PC-SAFT equation of state. In addition to the hard chain and dispersion terms, these functions also include dipole, association and ion terms for use with these types of compounds. When the ion term is included it is also called electrolyte PC-SAFT (ePC-SAFT).

Documentation
-------------
Documentation for the package is available on `Read the Docs`_.

Example
-------
.. code-block:: python

  import numpy as np
  from pcsaft import pcsaft_den

  # Toluene
  x = np.asarray([1.])
  m = np.asarray([2.8149])
  s = np.asarray([3.7169])
  e = np.asarray([285.69])
  pyargs = {'m':m, 's':s, 'e':e}

  t = 320 # K
  p = 101325 # Pa
  den = pcsaft_den(t, p, x, pyargs, phase='liq')
  print('Density of toluene at {} K: {} mol m^-3'.format(t, den))

  # Water using default 2B association scheme
  x = np.asarray([1.])
  m = np.asarray([1.2047])
  e = np.asarray([353.95])
  volAB = np.asarray([0.0451])
  eAB = np.asarray([2425.67])

  t = 274
  p = 101325
  s = np.asarray([2.7927 + 10.11*np.exp(-0.01775*t) - 1.417*np.exp(-0.01146*t)]) # temperature dependent sigma is used for better accuracy
  pyargs = {'m':m, 's':s, 'e':e, 'e_assoc':eAB, 'vol_a':volAB}
  den = pcsaft_den(t, p, x, pyargs, phase='liq')
  print('Density of water at {} K: {} mol m^-3'.format(t, den))
  
  # Water using 4C association scheme
  x = np.asarray([1.])
  m = np.asarray([1.2047])
  e = np.asarray([353.95])
  volAB = np.asarray([0.0451])
  eAB = np.asarray([2425.67])
  assoc_schemes = ['4c']

  t = 274
  p = 101325
  s = np.asarray([2.7927 + 10.11*np.exp(-0.01775*t) - 1.417*np.exp(-0.01146*t)]) # temperature dependent sigma is used for better accuracy
  pyargs = {'m':m, 's':s, 'e':e, 'e_assoc':eAB, 'vol_a':volAB, 'assoc_scheme':assoc_schemes}
  den = pcsaft_den(t, p, x, pyargs, phase='liq')
  print('Density of water at {} K: {} mol m^-3'.format(t, den))

Dependencies
------------

The runtime package depends on NumPy and SciPy. Building from source also requires Cython_ and a C++ toolchain. The repository vendors the Eigen_ headers used by the extension build.

Python package
--------------

The supported install paths are now standard pip installs:

::

  pip install .

or, for local development inside a source checkout:

::

  pip install -e . --no-build-isolation

This repository builds wheels in CI and supports source builds from the packaged Cython/C++ sources plus the vendored Eigen headers.

Source layout
-------------

The installable runtime code now lives under ``src/pcsaft/``. Package-owned runtime parameter datasets are bundled with the package under ``src/pcsaft/data/pcsaft_parameters`` and loaded through ``pcsaft.parameters``. The broader analysis datasets remain under the top-level ``data/`` workspace.

For local rebuilds in this repo, use:

::

  python scripts/build_pcsaft.py

That command performs an editable install using the active environment rather than relying on the legacy ``setup.py build_ext --inplace`` flow. If the editable build is already current, it skips the reinstall. Use ``python scripts/build_pcsaft.py --force`` to force a fresh editable reinstall.

Analysis and validation scripts are expected to run from the active ``PC-SAFT`` environment with ``pcsaft`` installed editable. A source checkout by itself is not a supported package-import path.

Author
------

- **Zach Baird** - zmeri_

License
-------

This project is licensed under the GNU General Public License v3.0

Acknowledgments
---------------

When developing these functions the code from two other groups was used as references

- Code from Joachim Gross (https://www.th.bci.tu-dortmund.de/cms/de/Forschung/PC-SAFT/Download/index.html)
- The MATLAB/Octave program written by Angel Martin and others (http://hpp.uva.es/open-source-software-eos/)

.. _`Clapeyron.jl`: https://github.com/ClapeyronThermo/Clapeyron.jl
.. _`teqp`: https://github.com/usnistgov/teqp
.. _`Read the Docs`: https://pcsaft.readthedocs.io/en/latest/
.. _Cython: http://cython.org/
.. _Eigen: https://github.com/eigenteam/eigen-git-mirror
.. _pcsaft: https://pypi.org/project/pcsaft/
.. _`Cython documentation`: http://docs.cython.org/en/latest/src/quickstart/build.html
.. _`earlier version`: https://github.com/zmeri/PC-SAFT/tree/b43bf568c4dc1907316422d5c3f7b809e9725848
.. _zmeri: https://github.com/zmeri
