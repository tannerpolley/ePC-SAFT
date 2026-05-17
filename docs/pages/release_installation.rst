Installation
============

Current release: ``1.5.2``

Install from PyPI
-----------------

The standard install command is:

.. code-block:: powershell

   python -m pip install epcsaft

With ``uv``:

.. code-block:: powershell

   uv add epcsaft

The current public release is also available from GitHub:

``https://github.com/tannerpolley/ePC-SAFT/releases/tag/v1.5.2``

Install from a wheel
--------------------

If the release includes a wheel matching your Python version and platform,
download it and install it directly:

.. code-block:: powershell

   python -m pip install C:\path\to\epcsaft-1.5.2-*.whl

This is the simplest install path because the native extension is already
built for your platform.

Install from tagged source
--------------------------

If no compatible wheel is available, install from the tagged Git source:

.. code-block:: powershell

   python -m pip install "epcsaft @ git+https://github.com/tannerpolley/ePC-SAFT.git@v1.5.2"

With ``uv``:

.. code-block:: powershell

   uv add "epcsaft @ git+https://github.com/tannerpolley/ePC-SAFT.git@v1.5.2"

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

   cd C:\path\to\ePC-SAFT-1.5.2
   python -m pip install .

Editable source install
-----------------------

Use an editable install when you are changing Python files and want imports to
come directly from the checkout:

.. code-block:: powershell

   git clone https://github.com/tannerpolley/ePC-SAFT.git
   cd ePC-SAFT
   python -m pip install -e .

With ``uv``:

.. code-block:: powershell

   uv pip install -e .

Editable installs use the same native build backend as wheel installs. Python
source changes are picked up from the checkout. If you change C++ sources,
pybind bindings, CMake files, or build metadata, rerun the editable install
command so the native extension is rebuilt.

Local path dependency
---------------------

For a project that depends on a local checkout, use a path dependency:

.. code-block:: toml

   dependencies = [
       "epcsaft @ file:///C:/path/to/ePC-SAFT",
   ]

After package changes, refresh the installed dependency:

.. code-block:: powershell

   uv sync --reinstall-package epcsaft

Use ``uv run --no-sync ...`` for follow-up downstream commands when you do not
want an implicit sync to rebuild the package again.

Optional IPOPT support
----------------------

IPOPT support is a native build dependency, not a Python extra. Request it only
when the platform already has a native Ipopt package that CMake can discover
through an install root or package config directory:

.. code-block:: powershell

   $env:EPCSAFT_PEP517_ENABLE_IPOPT = "1"
   $env:EPCSAFT_PEP517_IPOPT_ROOT = "C:\path\to\Ipopt"
   python -m pip install "epcsaft @ git+https://github.com/tannerpolley/ePC-SAFT.git@v1.5.2"

Use ``EPCSAFT_PEP517_IPOPT_DIR`` instead when the install provides an
``IpoptConfig.cmake`` directory.
Runtime processes that execute Ipopt on Windows must expose the Ipopt ``bin``
directory through both ``PATH`` and ``EPCSAFT_RUNTIME_DLL_DIRS``.

IPOPT remains opt-in at build time. When Ipopt is compiled, the public
``solver_backend="auto"`` selector uses the native ideal reactive-speciation
route where its assumptions hold. Broader public equilibrium routes still
require native route-builder work before they can use Ipopt as a production
solver.

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
