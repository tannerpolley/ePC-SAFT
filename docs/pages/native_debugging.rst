Native And Equation Debugging
=============================

Use this page when a change crosses the Python wrapper, pybind11 seam, or native ePC-SAFT equations.

Runtime flow
------------

The public API starts in ``src/epcsaft/epcsaft.py``. ``ePCSAFTMixture`` normalizes parameter payloads and owns a native ``NativeMixture``. ``ePCSAFTState`` validates ``T/x/P`` or ``T/x/rho`` inputs, constructs ``NativeState``, and delegates pressure, density, residual Helmholtz, fugacity, and activity-coefficient calls to the native extension.

The pybind11 boundary is ``src/epcsaft/bindings.cpp``. It exposes ``NativeArgs``, ``NativeMixture``, ``NativeState``, contribution-result structs, and native regression helpers through the private ``epcsaft._core`` module.

The native implementation lives under ``src/epcsaft/native``. High-traffic files are:

- ``epcsaft_density.cpp`` for pressure-to-density closure.
- ``epcsaft_ares.cpp`` for residual Helmholtz contribution totals.
- ``epcsaft_Z.cpp`` for compressibility factor and pressure from density.
- ``epcsaft_mu.cpp`` and ``epcsaft_fugcoef.cpp`` for residual chemical potential and fugacity.
- ``epcsaft_activity.cpp`` for activity, osmotic, and solvation outputs.

Validation commands
-------------------

Use the normal build path first:

.. code-block:: powershell

   uv run python scripts/build_epcsaft.py
   uv run python scripts/codex_doctor.py
   uv run python run_pytest.py --confidence -q

For C++ iteration after the build tree is configured:

.. code-block:: powershell

   uv run python scripts/build_epcsaft.py --build-only --parallel 10
   uv run python run_pytest.py --runtime -q

For method-speed checks:

.. code-block:: powershell

   uv run python run_pytest.py --profile -q

Equation traceability
---------------------

``docs/latex/equations.tex`` is the source of truth for equation text. ``docs/equations.md`` and ``docs/equations_registry.yaml`` are generated navigation views.

Native owner comments use ``// EqID: <id>`` near the implementing C++ function. When touching equation code, keep the EqID comment close to the function that owns the expression and run:

.. code-block:: powershell

   uv run python scripts/sync_equation_registry.py --check --strict-traceability
   uv run python run_pytest.py tests/test_equation_registry.py -q

See :doc:`equation_traceability` for the EqID classification and owner-comment checklist.

Debugging checklist
-------------------

- Reproduce the behavior through a public ``ePCSAFTMixture`` / ``ePCSAFTState`` call before debugging private native functions.
- Compare pressure-created and density-created states when investigating density closure.
- Request contribution terms with ``return_contribution_terms=True`` when debugging residual Helmholtz, compressibility factor, chemical potential, or fugacity totals.
- Use ``tests/test_native_runtime_contracts.py`` for fast neutral and ionic contribution-map regression coverage.
