Codex Workflow Guide
====================

This is the operational guide for future Codex agents and maintainers working from this source tree. Use these commands before inventing new build, test, or package workflows.

Default new-agent sequence
--------------------------

Start every new Codex thread with this sequence:

.. code-block:: powershell

   uv sync --no-install-project
   uv run python scripts/build_epcsaft.py
   uv run python scripts/codex_check.py quick

This is the expected healthy baseline. It creates the uv environment, builds the in-place pybind11 ``epcsaft._core`` extension, verifies imports/tool paths, generated-output state, and plot manifest presence/schema through doctor, then runs the fast contract suite. The default test route intentionally samples representative API, native, regression, equilibrium, and workflow contracts instead of running full equilibrium/regression reproductions or generated plot production. Use ``uv run python scripts/codex_check.py confidence`` before handoff when extra native runtime contracts and plot manifest path validation should be included.

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
     - ``uv run python scripts/codex_check.py confidence``
     - Before claiming repo runtime confidence. This includes doctor, the confidence slice, and the plot manifest check.
   * - Fast generic validation
     - ``uv run python scripts/codex_check.py quick``
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
   * - Electrolyte LLE confidence
     - ``uv run python run_pytest.py --equilibrium-confidence -q -s``
     - Bounded Khudaida fixture plus cached fixed-phase residual contract. Native confidence solving and full report generation remain explicit opt-ins.
   * - Docs check
     - ``uv run python scripts/codex_check.py docs``
     - Validate the tracked plot manifest and build Sphinx HTML under ``build/docs-html``.
   * - Plot asset rebuild
     - ``uv run python scripts/codex_check.py plots``
     - Opt-in generated plot producers, manifest refresh, and static asset report.
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

LaTeX and Overleaf mirror
-------------------------

``docs/latex`` is normal tracked repo content and is the source of truth for equation-heavy LaTeX files. It is not a Git submodule.

Use this once to create or validate the external Overleaf checkout:

.. code-block:: powershell

   .\scripts\setup_latex_mirror.ps1

After LaTeX edits are committed or ready to publish, mirror the current ``docs/latex`` tree to Overleaf:

.. code-block:: powershell

   .\scripts\sync_latex_mirror.ps1

The mirror lives at ``C:\Users\Tanner\Documents\git\LaTeX-Projects\ePC-SAFT-LaTeX`` and owns the Overleaf Git remote. The sync script copies the current LaTeX source tree and intentional top-level artifacts; generated ``docs/latex/out`` build products remain ignored in this repo.

Parallel agent safety
---------------------

The dev build tree and temp/profile outputs under ``build/`` are shared disposable state. In multi-agent sessions, keep native rebuild, clean, and repair coordination on the main thread unless that work is explicitly assigned.

- Do not run clean or repair actions while tests, REPLs, IDE run configurations, or other agents may import ``epcsaft._core``.
- Prefer one native builder at a time for ``build/dev`` and the in-place ``_core`` extension.
- Let sub-agents run focused test slices for their lane, and reserve full build, doctor, and ``--confidence`` validation for coordinated handoff checks.
- Use ``uv run python run_pytest.py --profile -q`` for quick runtime-only speed claims. Use ``uv run python run_pytest.py --profile-full -q -s`` before broad method-speed claims, allow at least 120 seconds, then read ``build/runtime_profile/*.md`` before reporting conclusions.

Project-local Git worktrees
---------------------------

Use ``scripts/create_codex_worktree.ps1`` from the primary checkout instead of raw ``git worktree add`` when a Codex agent needs a project-local worktree under ``.worktrees/``. The helper creates ``.worktrees/<name>`` and registers the new checkout path as a Git ``safe.directory`` so future Git commands inside that worktree do not fail on Windows because Codex runs as a sandbox user.

.. code-block:: powershell

   .\scripts\create_codex_worktree.ps1 -Name equilibrium-v3 -Branch codex/equilibrium-v3

``.worktrees/`` must stay ignored in ``.gitignore`` before using the helper. The ``safe.directory`` registration is intentionally path-specific and global to the current Windows user; Codex may still ask for approval because that is a persistent Git trust change. Use ``-SkipSafeDirectory`` only when you plan to use per-command ``git -c safe.directory=<path> ...`` overrides instead.

Test selection rules
--------------------

Use the smallest relevant test first, then run ``scripts/codex_check.py confidence`` before handoff. Use ``uv run python run_pytest.py --all -q`` only when you explicitly need the exhaustive historical suite.

- Python wrapper/API changes: ``uv run python run_pytest.py --api -q`` first, then ``uv run python run_pytest.py --confidence -q``.
- Native/equation changes: ``uv run python scripts/build_epcsaft.py --build-only --parallel 10`` first, then ``uv run python run_pytest.py --runtime -q``, then ``uv run python run_pytest.py --confidence -q``.
- Equation traceability changes: ``uv run python scripts/sync_equation_registry.py --check --strict-traceability`` then ``uv run python run_pytest.py tests/native/test_equation_registry.py -q``.
- Performance claims: ``uv run python run_pytest.py --profile -q -s`` is the quick runtime-only profile; use ``uv run python run_pytest.py --profile-full -q -s`` only for broad speed claims. Read the generated ``build/runtime_profile/*.md`` reports. Do not rely on skipped profile tests or code inspection alone.
- Plot/manifest changes: ``uv run python scripts/build_plot_manifest.py --check`` first. Use ``uv run python scripts/codex_check.py plots`` only when regenerating local plot outputs is explicitly part of the task.

``--profile`` is the quick runtime-only profile. ``--profile-full`` runs runtime, MIAC, and regression profiles and is the preferred evidence path for comprehensive speed reviews; use a timeout of at least 120 seconds.
- Packaging changes: ``scripts/build_dist.py``.

Keep generated plot assets and generated CSV workflows out of normal validation unless the task explicitly asks for them. Use ``--plots`` or ``scripts/codex_check.py plots`` for the opt-in plot asset slice.

Use ``uv run python run_pytest.py --list-slices`` when you need to inspect what each named slice runs before choosing a validation command.

For parallel Codex sessions, leave the default repo-local temp behavior alone unless it causes contention. When running concurrent pytest lanes, set ``EPCSAFT_PYTEST_TEMP_ROOT`` to an external temp root for the extra lanes so each run gets an isolated ``pytest-temp`` child.

Runtime speed rule
------------------

For repeated runtime calls, build ``ePCSAFTMixture`` and ``ePCSAFTState`` once and reuse them inside hot loops. The quick profile report compares reused-state activity-coefficient calls against full rebuild calls and flags the ratio when rebuilds dominate runtime.

Troubleshooting
---------------

Run ``uv run python scripts/codex_doctor.py`` whenever imports, tool paths, ``_core`` state, generated-output tracking, or plot manifests are unclear. It reports the active Python, git ref, uv/cmake/ninja paths, ``epcsaft`` import path, ``epcsaft._core`` path, required native symbol presence, generated artifact state, and the next recommended command.

If ``scripts/build_epcsaft.py`` appears slow, wait for the configured timeout before treating it as broken. Full configure/build can take far longer than the fast rebuild path; incremental ``--build-only --parallel 10`` is the intended C++ edit loop.
