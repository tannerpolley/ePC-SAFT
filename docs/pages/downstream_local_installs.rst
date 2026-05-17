Local Source Installs
=====================

This page is for projects that use a local ePC-SAFT checkout before or instead
of installing a published wheel.

Editable install
----------------

Use an editable install when you are changing Python files in the ePC-SAFT
checkout:

.. code-block:: powershell

   cd C:\path\to\ePC-SAFT
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

Use a path dependency when another project should install a local ePC-SAFT
checkout:

.. code-block:: toml

   dependencies = [
       "epcsaft @ file:///C:/Users/Tanner/Documents/git/ePC-SAFT",
   ]

Recommended local dependency loop
---------------------------------

Install or refresh the local package once, then prove the installed package and
the downstream repo-local integration wrapper without implicit sync:

.. code-block:: powershell

   $env:UV_CACHE_DIR = "$PWD\.uv-cache"
   uv sync --reinstall-package epcsaft
   uv run --no-sync python -m epcsaft
   uv run --no-sync python scripts/check_epcsaft_integration.py --mode dev

The repo's current CI path-install smoke uses Python 3.13. The package
metadata remains broader, but Python 3.13 is the baseline to match when
checking current local-install behavior.

Use ``uv run --no-sync`` after a known-good install because ordinary ``uv run`` is allowed to sync the environment. For a local native path dependency, that sync can rebuild ePC-SAFT when the downstream goal is only to run a smoke test.

Do not start multiple downstream ``uv run`` commands in parallel until the reinstall has completed. If parallel checks are needed, run the reinstall once first, then use ``uv run --no-sync`` for every parallel check.

The current downstream repos ``MEA-Thermodynamics``, ``Lithium_Extraction``,
and ``MEA-Absorption-Column`` all expose the repo-local install check
``scripts/check_epcsaft_integration.py --mode dev``. Run that first in the
downstream repo after reinstalling the local package, then run one real
workflow command from that repo separately. The real workflow runs are not
replaced by package-side smokes; issue #119 tracks them as a later release-gate
phase.

Build directory behavior
------------------------

PEP 517 wheel builds use an isolated temporary native build directory by default. This avoids repeated downstream path installs writing into the shared source checkout ``build/`` tree, which is the common source of Windows ``_core*.pyd`` lock races.

If you intentionally want a persistent build directory for a downstream reinstall, set:

.. code-block:: powershell

   $env:EPCSAFT_PEP517_BUILD_DIR = "$PWD\.uv-cache\epcsaft-build"
   uv sync --reinstall-package epcsaft

This keeps the package-install build tree around so a repeated full install can
reuse CMake/Ninja dependency state instead of rebuilding Ceres from scratch.

Reusable Ceres package
----------------------

Full package installs keep Ceres enabled by default. When repeated local or
downstream path installs spend most of their time compiling Ceres, build Ceres
once and point package installs at the prebuilt CMake package:

.. code-block:: powershell

   cd C:\path\to\ePC-SAFT
   uv run python scripts\dev\build_system_ceres.py --parallel 4

The helper prints the exact environment variables to reuse the result. The
manual form is:

.. code-block:: powershell

   $env:EPCSAFT_PEP517_CERES_DIR = "C:\path\to\ePC-SAFT\build\system-ceres\2.2.0\install\lib\cmake\Ceres"
   $env:EPCSAFT_PEP517_USE_SYSTEM_CERES = "1"
   $env:EPCSAFT_PEP517_BUILD_DIR = "$PWD\.uv-cache\epcsaft-build"
   uv sync --reinstall-package epcsaft

``EPCSAFT_PEP517_CERES_DIR`` must point at the directory containing
``CeresConfig.cmake``. The build backend then passes
``EPCSAFT_USE_SYSTEM_CERES=ON`` and ``Ceres_DIR=...`` to CMake. If the
environment variables are not set, package installs still use the package
default Ceres ``FetchContent`` path.

For normal ePC-SAFT source development, keep using the explicit in-place dev build:

.. code-block:: powershell

   uv run python scripts\dev\build_epcsaft.py
   uv run python scripts\dev\build_epcsaft.py --build-only --parallel 10

The default dev-script build is the required native dependency profile: Ceres ON, CppAD ON, and Ipopt OFF. Editable, wheel, and downstream path installs follow the same Ceres/CppAD requirement. New dev build configurations prefer Ninja when available. Existing ``build/dev`` trees keep their configured generator until you run the coordinated repair command ``uv run python scripts\dev\build_epcsaft.py --clean --generator ninja``.

Windows ``_core`` lock failures
-------------------------------

If a build reports ``Permission denied`` while writing ``_core*.pyd``, a Python process is usually still importing the extension. Stop downstream tests, Python REPLs, IDE runs, and parallel workers that imported ``epcsaft._core``. Then run exactly one reinstall/build command before starting checks again.

Runtime metadata
----------------

Downstream projects can confirm which package source they are using:

.. code-block:: python

   import epcsaft

   print(epcsaft.__version__)
   print(epcsaft.__git_commit__)
   print(epcsaft.runtime_build_info())

``runtime_build_info()`` reports package version, source path/commit when discoverable, native extension path, Python version, and platform information.

Capability discovery
--------------------

Use ``capabilities()`` before wiring high-level downstream workflows:

.. code-block:: python

   import epcsaft

   caps = epcsaft.capabilities()
   assert caps["equilibrium"]["neutral_tp_flash"]["status"] in {"available", "ipopt_dependency_required"}
   assert caps["equilibrium"]["electrolyte_bubble_pressure"]["status"] in {"available", "ipopt_dependency_required"}
   assert "reactive_electrolyte_bubble" not in caps["equilibrium"]
   assert "ReactiveElectrolyteBubbleProblem" in caps["equilibrium"]["problem_objects"]["classes"]

Native EOS/property calls and native regression helpers are available.
Implemented equilibrium capability entries are native Ipopt routes gated only
by whether the current install was compiled with Ipopt. Public route names that
still lack a production native route builder remain public contracts, but they
are not advertised as implemented capabilities.

For routing examples and the production/opt-in solver table, see
:doc:`equilibrium_cookbook`.

Capability status summary
~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1

   * - Area
     - Status
     - Notes
   * - Neutral TP flash and LLE
     - Native Ipopt route when compiled
     - Requires an Ipopt-enabled build; Python does not provide an alternate solve loop.
   * - Neutral bubble/dew pressure and temperature
     - Native Ipopt route when compiled
     - Requires an Ipopt-enabled build; Python does not provide an alternate solve loop.
   * - Neutral stability
     - Native Ipopt route when compiled
     - Uses native TPD route builders; Python does not provide an alternate solve loop.
   * - Electrolyte LLE
     - Native Ipopt route when compiled
     - Requires an Ipopt-enabled build.
   * - Reactive speciation
     - Explicit Ipopt ideal route
     - ``solver_backend="ipopt"`` supports homogeneous ``ideal_mole_fraction`` when Ipopt is compiled; activity and concentration routes require their native EOS derivative NLP blocks.
   * - Electrolyte bubble pressure
     - Native Ipopt route when compiled
     - Fixed liquid composition with neutral vapor species; ions remain liquid-only.
   * - Reactive electrolyte bubble
     - Staged native route when compiled
     - Requires native speciation followed by the native Ipopt fixed-liquid electrolyte bubble route.
   * - IPOPT
     - Optional native NLP backend
     - Owns production equilibrium solves as route builders land; Python does not provide an alternate optimizer path.

Package-side generic contract smoke coverage
--------------------------------------------

The upstream package tests cover three downstream-shaped generic API contracts
without adding downstream-specific public APIs:

* Reactive speciation with generic target rows for speciation, volatile partial
  pressure, and activity observations.
* Electrolyte LLE with generic solvent-feed, salt-molality, phase-composition,
  mean-ionic-activity, and regularization rows.
* Reactive electrolyte bubble pressure with generic speciation, fugacity, and
  partial-pressure rows.

Downstream projects should build their own project-specific metrics outside
``epcsaft``. The package boundary is the generic problem, result, capability,
and regression-target schema. Public API names are intentionally not tied to
MEA, lithium extraction, absorption columns, distribution coefficients,
selectivity, or other application labels.

These package-side tests do not count as the real downstream workflow proof
required by issue #119. They verify that a local install exposes the generic
package contracts that downstream repos consume. The required downstream proof
is still one real recorded workflow run each in ``MEA-Thermodynamics``,
``Lithium_Extraction``, and ``MEA-Absorption-Column`` after the local install
check has passed.

The smoke tests also assert that the public derivative contract keeps finite
difference route. Use ``analytic``, ``cppad``, ``analytic_implicit``, or
``cppad_implicit`` derivative routes where derivatives are required.
