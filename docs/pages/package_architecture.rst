Package Architecture
====================

``epcsaft`` is a package repo. The Python layer owns input validation,
user-facing objects, diagnostics, documentation examples, and workflow
orchestration. The equation-of-state runtime and package-owned phase,
chemical-equilibrium, and regression kernels are native C++ exposed through
``pybind11``.

Organization Boundary
---------------------

``epcsaft`` remains one installable distribution package. The project uses
clean internal subsystem boundaries instead of splitting EOS, equilibrium,
regression, native, data, or benchmark code into separate packages.

The target internal shape is:

.. code-block:: text

   src/epcsaft/
     eos/
     equilibrium/
     regression/
     native/
     data/

Subsystem Boundaries
--------------------

EOS harness
   Owns ``ePCSAFTMixture``, ``ePCSAFTState``, state construction, property
   evaluation wrappers, and the Python-facing equation-of-state contract. It
   validates user inputs and delegates thermodynamic calculations to the native
   runtime.

Equilibrium
   Owns phase-equilibrium, stability, bubble/dew, electrolyte LLE, and
   chemical-equilibrium orchestration. It may use Python for problem objects,
   request normalization, result shaping, and diagnostics, but production
   thermodynamic evaluations should route through the EOS/native boundary.
   The public ``mixture.equilibrium(kind=...)`` facade adapts non-reactive
   string requests into typed problem objects and shared route diagnostics;
   reactive-specialized routes remain explicit entrypoints with their own
   option checks.

Regression
   Owns fitting problem definitions, records, provenance validation, objective
   assembly, derivative diagnostics, and fit-result serialization. Public
   regression helpers remain Python-facing while expensive objective and
   derivative work should use native kernels when available.

Native
   Owns C++ kernels, pybind11 bindings, native capability reporting, and
   internal C++ helpers. Python code should call native functionality through
   the public runtime surfaces or thin package-owned adapters, not by reaching
   into build artifacts directly.

Data
   Owns packaged parameter datasets, dataset validation, reference-data loading
   contracts, and reusable package data. Analysis-local inputs belong under the
   relevant ``analyses/<category>/<id>/figures/<figure_id>/input`` tree instead of becoming hidden package
   dependencies.

Benchmarks
   Owns package-maintained timing, smoke, and regression benchmarks that protect
   runtime expectations. Benchmark execution lives under ``scripts/benchmarks``,
   may consume the public package API and packaged/reference data, and should
   not become a runtime package import dependency for normal users. The
   package-owned literature benchmark suite inventory also belongs there: it is
   the canonical classification surface for which literature anchors are already
   supported by tests versus still blocked on follow-up generic capability work.

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
* ``scripts/benchmarks`` for package-owned timing and smoke benchmarks.
* ``epcsaft.diagnostics`` for ``capabilities()`` and ``runtime_build_info()``.

Top-level imports remain stable for existing users. Benchmark execution helpers
are validation assets, not runtime thermodynamic APIs; import them from
``scripts.benchmarks`` rather than from the runtime package.

Import Policy
-------------

Public user code should import from the top-level package or from documented
subsystem modules:

* ``import epcsaft``
* ``from epcsaft import ePCSAFTMixture``
* ``from epcsaft.equilibrium import ...``
* ``from epcsaft.regression import ...``

Internal modules may share package-owned helpers when that keeps behavior
centralized, but subsystem code should avoid circular ownership. In particular,
benchmarks and docs may depend on public APIs, while core runtime modules must
not depend on benchmark entrypoints or generated analysis artifacts.

Compatibility Policy
--------------------

Current public imports must continue working across boundary cleanups:

.. code-block:: python

   import epcsaft

   epcsaft.ePCSAFTMixture
   epcsaft.solve_reactive_speciation
   epcsaft.fit_pure_neutral

Cleaner subsystem imports may be added over time, but large module moves must
land in small refactor PRs with stable facade imports and focused API tests.
Do not use package-boundary work as a reason to break downstream notebooks,
MEA/Li extraction consumers, or existing documented imports.

Optional Dependency Policy
--------------------------

The default install should keep the lightweight runtime usable. Heavy or
platform-sensitive dependencies belong behind optional dependency groups,
feature flags, or runtime capability checks. For example, Ipopt-dependent
workflows must fail with actionable diagnostics when native Ipopt is not
compiled or a route is outside the compiled native adapter surface, instead of
making the base package import fail.

Native build capabilities should be reported through ``capabilities()`` and
``runtime_build_info()`` so downstream projects can select supported workflows
without probing private modules.

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
inputs, and report diagnostics, but it should not
silently become the production thermodynamic solver.
