Downstream Local Installs
=========================

This page is for downstream projects that depend on a local ePC-SAFT checkout through a uv path dependency, for example:

.. code-block:: toml

   dependencies = [
       "epcsaft @ file:///C:/Users/Tanner/Documents/git/ePC-SAFT",
   ]

Recommended downstream loop
---------------------------

Install or refresh the local package once, then run downstream checks without implicit sync:

.. code-block:: powershell

   $env:UV_CACHE_DIR = "$PWD\.uv-cache"
   uv sync --reinstall-package epcsaft
   uv run --no-sync python -m unittest tests.test_epcsaft_ionic -v
   uv run --no-sync python -m MEA.epcsaft_ionic.smoke

The repo's current CI path-install smoke uses Python 3.13. The package metadata remains broader, but Python 3.13 is the baseline to match when checking current downstream local-install behavior.

Use ``uv run --no-sync`` after a known-good install because ordinary ``uv run`` is allowed to sync the environment. For a local native path dependency, that sync can rebuild ePC-SAFT when the downstream goal is only to run a smoke test.

Do not start multiple downstream ``uv run`` commands in parallel until the reinstall has completed. If parallel checks are needed, run the reinstall once first, then use ``uv run --no-sync`` for every parallel check.

Build directory behavior
------------------------

PEP 517 wheel builds use an isolated temporary native build directory by default. This avoids repeated downstream path installs writing into the shared source checkout ``build/`` tree, which is the common source of Windows ``_core*.pyd`` lock races.

If you intentionally want a persistent build directory for a downstream reinstall, set:

.. code-block:: powershell

   $env:EPCSAFT_PEP517_BUILD_DIR = "$PWD\.uv-cache\epcsaft-build"
   uv sync --reinstall-package epcsaft

For normal ePC-SAFT source development, keep using the explicit in-place dev build:

.. code-block:: powershell

   uv run python scripts\build_epcsaft.py
   uv run python scripts\build_epcsaft.py --build-only --parallel 10

New dev build configurations prefer Ninja when available. Existing ``build/dev`` trees keep their configured generator until you run the coordinated repair command ``uv run python scripts\build_epcsaft.py --clean --generator ninja``.

Editable installs
-----------------

The current custom build backend is intended for wheel/path installs and does not expose a PEP 660 editable-install workflow as the normal downstream path. For downstream co-development, prefer reinstalling the local path dependency after package changes:

.. code-block:: powershell

   $env:UV_CACHE_DIR = "$PWD\.uv-cache"
   uv sync --reinstall-package epcsaft
   uv run --no-sync python -m your_downstream_smoke

When working inside the ePC-SAFT checkout itself, use ``uv run python scripts\build_epcsaft.py --build-only --parallel 10`` after the initial build rather than relying on ``pip install -e``.

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
   assert caps["equilibrium"]["neutral_tp_flash"]["available"]
   assert caps["equilibrium"]["electrolyte_bubble_pressure"]["available"]
   assert caps["equilibrium"]["reactive_electrolyte_bubble"]["available"]

Native neutral TP flash, neutral LLE, electrolyte LLE, electrolyte bubble pressure, reactive speciation, reactive electrolyte bubble pressure, and native regression helpers are available. The electrolyte bubble-pressure workflow is scoped to fixed liquid composition with neutral vapor species; ions remain liquid-only. The reactive electrolyte bubble workflow performs native chemical speciation before the native fixed-liquid electrolyte bubble-pressure solve.

For agent-facing routing examples and the production/experimental solver table,
see :doc:`equilibrium_cookbook`.

Capability status summary
~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1

   * - Area
     - Status
     - Notes
   * - Neutral TP flash, LLE, stability
     - Production native
     - Use the explicit mixture methods in new code.
   * - Electrolyte LLE
     - Production native
     - Fixed-species, charge-neutral LLE; use stability and explicit seeds for hard cases.
   * - Reactive speciation
     - Production native
     - Homogeneous one-phase chemistry with caller-owned balances and reactions.
   * - Electrolyte bubble pressure
     - Production native, scoped
     - Fixed liquid composition and neutral vapor species; ions remain liquid-only.
   * - Reactive electrolyte bubble
     - Staged production native, scoped
     - Native speciation followed by native fixed-liquid electrolyte bubble pressure.
   * - IPOPT
     - Experimental opt-in
     - Residual-minimization refinement only; ``auto`` never selects IPOPT.

MEA benchmark scope
-------------------

``fit_mea_co2_h2o_electrolyte(...)`` is a fixed-composition native benchmark helper. It is useful for checking electrolyte pure-parameter regression plumbing, provenance, bounds, and native objective behavior.

It is not reactive bubble-pressure fitting. Fixed-liquid fugacity-coefficient residuals should not be interpreted as a full MEA absorber VLE regression objective, especially when fitted values hit bounds or residuals remain large.
