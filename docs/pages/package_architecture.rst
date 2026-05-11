Package Architecture
====================

``epcsaft`` is a package repo. The Python layer owns input validation,
user-facing objects, diagnostics, documentation examples, and workflow
orchestration. The equation-of-state runtime and package-owned phase,
chemical-equilibrium, and regression kernels are native C++ exposed through
``pybind11``.

Core Surfaces
-------------

Use these imports for new code:

* ``epcsaft.eos`` for ``ePCSAFTMixture`` and ``ePCSAFTState``. It also exports
  ``Mixture`` and ``State`` aliases for shorter new-user examples.
* ``epcsaft.equilibrium`` for neutral and electrolyte equilibrium helpers.
* ``epcsaft.electrolyte`` for electrolyte LLE and fixed-liquid bubble pressure.
* ``epcsaft.reactive`` for reactive speciation, reactive electrolyte bubble
  pressure, and reactive regression batch/context objects.
* ``epcsaft.regression`` for public fitting helpers.
* ``epcsaft.parameters`` for packaged parameter dataset loading.
* ``epcsaft.benchmarks`` for package-owned timing and smoke benchmarks.
* ``epcsaft.diagnostics`` for ``capabilities()`` and ``runtime_build_info()``.

Top-level imports remain stable for existing users. The organized modules are
navigation aids, not a breaking API move.

Repository Layout
-----------------

``src/epcsaft`` contains the package. ``tests`` contains package/API/native
contracts. ``data/reference`` is the canonical source-checkout reference-data
library. ``analyses`` contains paper-validation and analysis workflows, each
with local ``data`` and ``results`` folders.

Generated benchmark output, run payloads, build trees, and local graph or temp
outputs are not source artifacts.

Native-First Policy
-------------------

Package-owned regression, equilibrium, and speciation workflows must use native
runtime kernels for thermodynamic calculations. Python may batch rows, validate
inputs, manage continuation seeds, and report diagnostics, but it should not
silently become the production thermodynamic solver.

Native Build-Time Dependency Controls
-------------------------------------

``cmake`` options used by the package and workflow scripts are intentionally
explicit:

* ``EPCSAFT_ENABLE_CERES`` enables package-wide Ceres support for native residual
  and regression backends.
* ``EPCSAFT_USE_SYSTEM_CERES`` requests an installed Ceres package instead of
  the fallback fetch path.
* ``EPCSAFT_ENABLE_CPPAD`` enables package-wide CppAD support where native
  derivative workflows are active.
* ``EPCSAFT_USE_SYSTEM_CPPAD`` requests an installed CppAD package.

The public CLI (``scripts/build_epcsaft.py``) exposes equivalent switches:
``--enable-ceres``, ``--use-system-ceres``, ``--enable-cppad``, and
``--use-system-cppad``.
