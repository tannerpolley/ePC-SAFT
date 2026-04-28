Equation Traceability Checklist
===============================

Use this checklist when changing equation text, native equation code, or EqID ownership.

Classify EqIDs
--------------

``docs/latex/equations.tex`` is the only editable equation metadata source. Every EqID must be in one of these states:

- Implemented: the equation has at least one nearby native owner comment, ``// EqID: <id>``.
- Documentation-only: the equation is reference material, notation, derivation, or explanatory context with no direct native owner.

Do not hand-edit ``docs/equations.md`` or ``docs/equations_registry.yaml``. Regenerate them with ``scripts/sync_equation_registry.py``.

Place Owner Comments
--------------------

Place each ``// EqID: <id>`` comment immediately above the C++ function or expression block that owns the equation. Prefer the narrowest owner that future agents can inspect quickly, such as a contribution helper, density closure function, or activity-coefficient conversion block.

Use ``Documentation-only`` only when there is no direct implemented owner. Do not use it to hide missing traceability for active formulas.

Validate
--------

For equation work, run:

.. code-block:: powershell

   uv run python scripts/sync_equation_registry.py
   uv run python scripts/sync_equation_registry.py --check --strict-traceability
   uv run python run_pytest.py tests/test_equation_registry.py -q

Before handoff after native/equation changes, also run:

.. code-block:: powershell

   uv run python run_pytest.py --native -q
   uv run python run_pytest.py --runtime -q
   uv run python run_pytest.py --confidence -q

Equation-Family Coverage Matrix
-------------------------------

Native/equation coverage is tracked by equation family plus targeted edge-case
contracts. Do not require one test per EqID; require every implemented runtime
family to have at least one deterministic contract that exercises the public
Python API through the pybind/native path.

Current coverage families:

- density closure / root solving
- residual Helmholtz energy
- compressibility factor
- temperature derivative
- composition derivative
- residual chemical potential
- fugacity coefficient
- Debye-Huckel / Born / ionic activity coefficient
- reference-state and density cache behavior
