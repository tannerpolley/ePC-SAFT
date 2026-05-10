Release Installation
====================

Current release: ``1.5.0``

The current official download route is the GitHub release:

``https://github.com/tannerpolley/ePC-SAFT/releases/tag/v1.5.0``

Install from a wheel
--------------------

If the release includes a wheel matching your Python version and platform,
download it and install it directly:

.. code-block:: powershell

   python -m pip install C:\path\to\epcsaft-1.5.0-*.whl

This is the simplest install path because the native extension is already
built for your platform.

Install from tagged source
--------------------------

If no compatible wheel is available, install from the tagged Git source:

.. code-block:: powershell

   python -m pip install "epcsaft @ git+https://github.com/tannerpolley/ePC-SAFT.git@v1.5.0"

With ``uv``:

.. code-block:: powershell

   uv add "epcsaft @ git+https://github.com/tannerpolley/ePC-SAFT.git@v1.5.0"

Source builds require:

- Python ``>=3.9``.
- A C++ compiler for your platform.
- CMake.
- Ninja or another working CMake generator.
- Network access to download build requirements declared in ``pyproject.toml``.

Python 3.13 is the current project smoke-test baseline.

Install from a local source archive
-----------------------------------

Download the release source archive, extract it, then run:

.. code-block:: powershell

   cd C:\path\to\ePC-SAFT-1.5.0
   python -m pip install .

For downstream projects that use ``uv`` and a local path dependency, prefer a
normal reinstall after package changes:

.. code-block:: powershell

   uv sync --reinstall-package epcsaft
   uv run --no-sync python -m your_downstream_smoke

PyPI status
-----------

The package metadata is prepared for distribution under the name ``epcsaft``.
Once PyPI publishing is enabled, the install command will be:

.. code-block:: powershell

   python -m pip install epcsaft

Until then, use the GitHub release or tagged Git source install path.

Optional IPOPT support
----------------------

IPOPT support is optional and experimental. Install it only when the platform
already has the required IPOPT and ``cyipopt`` prerequisites:

.. code-block:: powershell

   python -m pip install "epcsaft[ipopt] @ git+https://github.com/tannerpolley/ePC-SAFT.git@v1.5.0"

IPOPT is explicit opt-in only. It is not selected by ``solver_backend="auto"``
and it is not a full constrained Gibbs/NLP replacement.

Verify the install
------------------

.. code-block:: python

   import epcsaft

   print(epcsaft.__version__)
   print(epcsaft.runtime_build_info())
   print(epcsaft.capabilities())

If import fails after a source build on Windows, make sure no Python process is
holding an old ``epcsaft._core`` extension open and reinstall from a clean
environment.
