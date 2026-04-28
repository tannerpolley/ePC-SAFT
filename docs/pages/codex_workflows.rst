Codex Workflow Guide
====================

This is the operational guide for future Codex agents and maintainers working from this source tree. Use these commands before inventing new build, test, or package workflows.

Default new-agent sequence
--------------------------

Start every new Codex thread with this sequence:

.. code-block:: powershell

   uv sync --no-install-project
   uv run python scripts/build_epcsaft.py
   uv run python scripts/codex_doctor.py
   uv run python run_pytest.py --confidence -q

This is the expected healthy baseline. It creates the uv environment, builds the in-place pybind11 ``epcsaft._core`` extension, verifies imports and tool paths, then runs the standard confidence suite.

Use ``uv run python run_pytest.py ...`` for repo validation. Direct ``uv run python -m pytest ...`` works, but the wrapper sets ``src`` on the import path and uses a per-run pytest temp directory that is safer for Codex and Windows runs.

Command matrix
--------------

.. list-table::
   :header-rows: 1
   :widths: 28 42 30

   * - Situation
     - Command
     - Use when
   * - First setup or uncertain state
     - ``uv sync --no-install-project`` then ``uv run python scripts/build_epcsaft.py`` then ``uv run python scripts/codex_doctor.py``
     - Starting a fresh thread, after dependency changes, or after a failed import.
   * - Handoff validation
     - ``uv run python run_pytest.py --confidence -q``
     - Before claiming repo runtime confidence. This includes the generic slice plus native runtime contracts.
   * - Fast generic validation
     - ``uv run python run_pytest.py --generic -q``
     - Quick checks for Python/runtime/regression API changes when native contract coverage is not required yet.
   * - Python API work
     - ``uv run python run_pytest.py --api -q``
     - Public wrapper, parameter-template, or regression API edits.
   * - Native or density/equation work
     - ``uv run python scripts/build_epcsaft.py --build-only --parallel 10`` then ``uv run python run_pytest.py --runtime -q``
     - C++ iteration after ``build/dev`` is already configured.
   * - Native contract only
     - ``uv run python run_pytest.py --native -q``
     - Fast check for pressure-vs-density and contribution-map contracts.
   * - Quick method-speed check
     - ``uv run python run_pytest.py --profile -q``
     - Runtime-only profiling. The wrapper enables the required performance environment flag.
   * - Full method-speed check
     - ``uv run python run_pytest.py --profile-full -q -s``
     - Comprehensive runtime, MIAC, and regression profiling before making broad speed claims. This can take about a minute locally; allow at least 120 seconds.
   * - Package boundary
     - ``uv run python scripts/build_dist.py``
     - Wheel/sdist and smoke-import validation.
   * - Installed/source diagnostic
     - ``uv run python -m epcsaft``
     - Confirm package and ``epcsaft._core`` paths.

Build rules
-----------

The canonical native build command is:

.. code-block:: powershell

   uv run python scripts/build_epcsaft.py

Use ``--build-only --parallel 10`` only after the CMake tree already exists. Use ``--configure-only`` when you need to refresh CMake configuration without compiling.

Do not use ``--clean`` for routine validation. ``uv run python scripts/build_epcsaft.py --clean`` is a repair action for stale CMake state or stale/locked ``_core`` artifacts. If Windows reports that ``_core*.pyd`` is locked, stop Python REPLs, tests, IDE run configurations, or Codex sub-agents that imported ``epcsaft._core`` before retrying.

Parallel agent safety
---------------------

The dev build tree and temp/profile outputs under ``build/`` are shared disposable state. In multi-agent sessions, keep native rebuild, clean, and repair coordination on the main thread unless that work is explicitly assigned.

- Do not run clean or repair actions while tests, REPLs, IDE run configurations, or other agents may import ``epcsaft._core``.
- Prefer one native builder at a time for ``build/dev`` and the in-place ``_core`` extension.
- Let sub-agents run focused test slices for their lane, and reserve full build, doctor, and ``--confidence`` validation for coordinated handoff checks.
- Use ``uv run python run_pytest.py --profile -q`` for quick runtime-only speed claims. Use ``uv run python run_pytest.py --profile-full -q -s`` before broad method-speed claims, allow at least 120 seconds, then read ``build/runtime_profile/*.md`` before reporting conclusions.

Test selection rules
--------------------

Use the smallest relevant test first, then run ``--confidence`` before handoff.

- Python wrapper/API changes: ``uv run python run_pytest.py --api -q`` first, then ``uv run python run_pytest.py --confidence -q``.
- Native/equation changes: ``uv run python scripts/build_epcsaft.py --build-only --parallel 10`` first, then ``uv run python run_pytest.py --runtime -q``, then ``uv run python run_pytest.py --confidence -q``.
- Equation traceability changes: ``uv run python scripts/sync_equation_registry.py --check --strict-traceability`` then ``uv run python run_pytest.py tests/test_equation_registry.py -q``.
- Performance claims: ``uv run python run_pytest.py --profile -q -s`` is the quick runtime-only profile; use ``uv run python run_pytest.py --profile-full -q -s`` only for broad speed claims. Read the generated ``build/runtime_profile/*.md`` reports. Do not rely on skipped profile tests or code inspection alone.

``--profile`` is the quick runtime-only profile. ``--profile-full`` runs runtime, MIAC, and regression profiles and is the preferred evidence path for comprehensive speed reviews; use a timeout of at least 120 seconds.
- Packaging changes: ``scripts/build_dist.py``.

Keep generated plot/gallery and generated CSV workflows out of normal validation unless the task explicitly asks for them.

Use ``uv run python run_pytest.py --list-slices`` when you need to inspect what each named slice runs before choosing a validation command.

For parallel Codex sessions, leave the default repo-local temp behavior alone unless it causes contention. When running concurrent pytest lanes, set ``EPCSAFT_PYTEST_TEMP_ROOT`` to an external temp root for the extra lanes so each run gets an isolated ``pytest-temp`` child.

Runtime speed rule
------------------

For repeated runtime calls, build ``ePCSAFTMixture`` and ``ePCSAFTState`` once and reuse them inside hot loops. The quick profile report compares reused-state activity-coefficient calls against full rebuild calls and flags the ratio when rebuilds dominate runtime.

Troubleshooting
---------------

Run ``uv run python scripts/codex_doctor.py`` whenever imports, tool paths, or ``_core`` state are unclear. It reports the active Python, git ref, uv/cmake/ninja paths, ``epcsaft`` import path, ``epcsaft._core`` path, and the next recommended command.

If ``scripts/build_epcsaft.py`` appears slow, wait for the configured timeout before treating it as broken. Full configure/build can take far longer than the fast rebuild path; incremental ``--build-only --parallel 10`` is the intended C++ edit loop.
