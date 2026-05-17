Development Workflow Guide
==========================

This is the operational guide for maintainers working from this source tree. Use these commands before inventing new build, test, or package workflows.

Default source-checkout sequence
--------------------------------

Start every fresh source checkout with this sequence:

.. code-block:: powershell

   uv sync --no-install-project
   uv run python scripts/dev/build_epcsaft.py
   uv run python scripts/dev/doctor.py
   uv run python scripts/dev/validate_project.py quick

This is the expected healthy baseline. It creates the uv environment, builds the in-place pybind11 ``epcsaft._core`` extension, verifies imports/tool paths and generated-output state through doctor, then runs the fast contract suite. The default test route intentionally samples representative API, native, regression, equilibrium, and workflow contracts instead of running full equilibrium/regression reproductions or generated plot production. Use ``uv run python scripts/dev/validate_project.py confidence`` before release or broad runtime claims when extra native runtime contracts should be included.
The current development and CI smoke baseline is Python 3.13, while ``pyproject.toml`` still declares package compatibility with Python ``>=3.9``.

Use ``uv run python run_pytest.py ...`` for repo validation. Direct ``uv run python -m pytest ...`` works, but the wrapper sets ``src`` on the import path and uses a per-run pytest temp directory that is safer for Windows and parallel local runs.

Command matrix
--------------

.. list-table::
   :header-rows: 1
   :widths: 28 42 30

   * - Situation
     - Command
     - Use when
   * - First setup or uncertain state
     - ``uv sync --no-install-project`` then ``uv run python scripts/dev/build_epcsaft.py`` then ``uv run python scripts/dev/doctor.py``
     - Starting a fresh thread, after dependency changes, or after a failed import.
   * - Handoff validation
     - ``uv run python scripts/dev/validate_project.py confidence``
     - Before claiming repo runtime confidence. This includes doctor and the confidence slice.
   * - Fast generic validation
     - ``uv run python scripts/dev/validate_project.py quick``
     - Quick checks for Python/runtime/regression API changes when native contract coverage is not required yet.
   * - Python API work
     - ``uv run python run_pytest.py --api -q``
     - Public wrapper, parameter-template, or regression API edits.
   * - Equilibrium/speciation workflows
     - ``uv run python run_pytest.py --equilibrium-api -q``
     - Fast representative check for neutral equilibrium, electrolyte LLE, reactive speciation, reactive electrolyte bubble contracts, derivative-backend contracts, and capability reporting.
   * - Native or density/equation work
     - ``uv run python scripts/dev/build_epcsaft.py --build-only --parallel 10`` then ``uv run python run_pytest.py --runtime -q``
     - C++ iteration after ``build/dev`` is already configured.
   * - Native contract only
     - ``uv run python run_pytest.py --native -q``
     - Fast check for pressure-vs-density and contribution-map contracts.
   * - Electrolyte LLE confidence
     - ``uv run python run_pytest.py --equilibrium-confidence -q -s``
     - Bounded Khudaida fixture plus cached fixed-phase residual contract. Native confidence solving and full report generation remain explicit opt-ins.
   * - Docs check
     - ``uv run python scripts/dev/validate_project.py docs``
     - Build Sphinx HTML under ``build/docs-html``.
   * - Quick method-speed check
     - ``uv run python run_pytest.py --profile -q``
     - Runtime-only profiling. The wrapper enables the required performance environment flag.
   * - Neutral equilibrium benchmark
     - ``uv run python scripts/benchmarks/benchmark_neutral_equilibrium.py --warmup 20 --repeat 100``
     - Measure neutral state, TP flash, and seeded neutral LLE without any FeOs dependency.
   * - Literature benchmark inventory
     - ``uv run python scripts/benchmarks/benchmark_literature_suite.py``
     - Review the package-owned literature benchmark scope, including which issue anchors are already supported by tests and which still require follow-up work.
   * - Full method-speed check
     - ``uv run python run_pytest.py --profile-full -q -s``
     - Comprehensive runtime, MIAC, and regression profiling before making broad speed claims. This can take about a minute locally; allow at least 120 seconds.
   * - Package boundary
     - ``uv run python scripts/dev/build_dist.py``
     - Wheel/sdist and smoke-import validation. Isolated package builds default to serial native compilation to avoid Windows Ceres memory spikes; use ``--parallel N`` only when the machine has enough headroom.
   * - Installed/source diagnostic
     - ``uv run python -m epcsaft``
     - Confirm package and ``epcsaft._core`` paths.

Build rules
-----------

The canonical local native build command is:

.. code-block:: powershell

   uv run python scripts/dev/build_epcsaft.py

That command uses ``--profile fast`` by default: Ceres and CppAD are enabled, and Ipopt is disabled unless explicitly requested. Ceres is required for native regression builds, and CppAD is required for derivative-capable package builds.

Wheel/editable/path installs go through the PEP 517/scikit-build backend and use the same required Ceres/CppAD policy. Repeated clean builds can still avoid rebuilding Ceres by using a prebuilt system Ceres package.

.. list-table::
   :header-rows: 1
   :widths: 16 22 22 40

   * - Profile
     - Native options
     - Default Windows parallelism
     - Use when
   * - ``fast``
     - Ceres ON, CppAD ON
     - ``4``
     - Normal source-checkout setup, C++ iteration, reactive/speciation work, and most validation.
   * - ``full``
     - Ceres ON, CppAD ON
     - ``4``
     - Alias for the required Ceres + CppAD dependency profile.
   * - ``ipopt``
     - Ceres ON, CppAD ON, Ipopt ON
     - ``4``
     - Native Ipopt adapter development or validation with a system Ipopt package.

Use ``--build-only --parallel 10`` only after the CMake tree already exists. ``--build-only`` does not reconfigure profile flags; it builds whatever ``build/dev/CMakeCache.txt`` already says. Use ``--configure-only`` when you need to refresh CMake configuration without compiling. For a new ``build/dev`` tree, ``scripts/dev/build_epcsaft.py`` now prefers Ninja when ``ninja`` is available on ``PATH`` because it is usually faster than MinGW Makefiles for repeated local rebuilds. Existing CMake trees keep their original generator; doctor reports ``build_generator_recommendation`` when ``uv run python scripts/dev/build_epcsaft.py --clean --generator ninja`` is the appropriate one-time migration from an older MinGW tree.

Every native build writes ``build/dev/build_epcsaft.log`` and finishes with an ``epcsaft._core`` import check when compilation runs. Use ``uv run python scripts/dev/build_epcsaft.py --status`` when you need a non-mutating check of the configured generator, Ceres/CppAD flags, system-Ceres/Ceres_DIR state, importable ``_core`` artifacts, stale ``.ninja_lock`` state, last Ninja target, and live repo-owned build processes. This is the safest first check when an IDE run or interrupted terminal build appears hung.

For IDE run configurations, keep commands explicit instead of relying on one overloaded target:

- Native build status: ``uv run python scripts/dev/build_epcsaft.py --status``
- Native incremental build: ``uv run python scripts/dev/build_epcsaft.py --build-only --parallel 10``
- Native configure/build: ``uv run python scripts/dev/build_epcsaft.py --profile fast``
- Clean Ceres + CppAD proof: ``uv run python scripts/dev/build_epcsaft.py --clean --profile full --parallel 4``
- Native Ipopt proof with an install root: ``uv run python scripts/dev/build_epcsaft.py --clean --profile ipopt --ipopt-root C:\path\to\Ipopt``
- Native Ipopt proof with an ``IpoptConfig.cmake`` directory: ``uv run python scripts/dev/build_epcsaft.py --clean --profile ipopt --ipopt-dir C:\path\to\lib\cmake\Ipopt``

Do not use ``--clean`` for routine validation. ``uv run python scripts/dev/build_epcsaft.py --clean`` is a repair action for stale CMake state or stale/locked ``_core`` artifacts. A clean dev build deletes the reusable CMake tree, so Ceres source/configuration/build work under ``build/dev/_deps`` may run again unless you use a prebuilt system Ceres package. If Windows reports that ``_core*.pyd`` is locked, stop Python REPLs, tests, IDE run configurations, or parallel workers that imported ``epcsaft._core`` before retrying. If ``--status`` reports a stale Ninja lock, inspect the listed process ids and stop only repo-owned build processes before retrying.

If Ceres becomes part of a repeated local workflow, build or install Ceres once outside ``build/dev`` and use the system-Ceres path instead of vendoring it through ``FetchContent`` on every clean full-profile configure. The supported command shape is:

.. code-block:: powershell

   uv run python scripts/dev/build_system_ceres.py --parallel 4
   uv run python scripts/dev/build_epcsaft.py --profile full --use-system-ceres --ceres-dir C:\path\to\lib\cmake\Ceres

``--ceres-dir`` should point at the directory containing ``CeresConfig.cmake``. Ceres' own CMake documentation supports consuming either an installed Ceres package or an exported Ceres build directory through ``find_package(Ceres)`` and ``Ceres::ceres``.

The same prebuilt Ceres package can accelerate downstream path installs without changing package defaults:

.. code-block:: powershell

   $env:EPCSAFT_PEP517_CERES_DIR = "C:\path\to\lib\cmake\Ceres"
   $env:EPCSAFT_PEP517_USE_SYSTEM_CERES = "1"
   $env:EPCSAFT_PEP517_BUILD_DIR = "$PWD\.uv-cache\epcsaft-build"
   uv sync --reinstall-package epcsaft

Without those environment variables, wheel/editable/path installs still build the default Ceres-enabled package through ``FetchContent``.

LaTeX and Overleaf mirror
-------------------------

``docs/latex`` is normal tracked repo content and is the source of truth for equation-heavy LaTeX files. It is not a Git submodule.

Use this once to create or validate the external Overleaf checkout:

.. code-block:: powershell

   .\scripts\docs\setup_latex_mirror.ps1

After LaTeX edits are committed or ready to publish, mirror the current ``docs/latex`` tree to Overleaf:

.. code-block:: powershell

   .\scripts\docs\sync_latex_mirror.ps1

The mirror lives at ``C:\Users\Tanner\Documents\git\LaTeX-Projects\ePC-SAFT-LaTeX`` and owns the Overleaf Git remote. The sync script copies the current LaTeX source tree and intentional top-level artifacts; generated ``docs/latex/out`` build products remain ignored in this repo.

Parallel worker safety
----------------------

The dev build tree and temp/profile outputs under ``build/`` are shared disposable state. In parallel sessions, coordinate native rebuild, clean, and repair work so only one process owns the native extension at a time.

- Do not run clean or repair actions while tests, REPLs, IDE run configurations, or other workers may import ``epcsaft._core``.
- Prefer one native builder at a time for ``build/dev`` and the in-place ``_core`` extension.
- Let parallel workers run focused test slices for their lane, and reserve full build, doctor, and ``--confidence`` validation for coordinated handoff checks.
- Use ``uv run python run_pytest.py --profile -q`` for quick runtime-only speed claims. Use ``uv run python run_pytest.py --profile-full -q -s`` before broad method-speed claims, allow at least 120 seconds, then read ``build/runtime_profile/*.md`` before reporting conclusions.

Project-local Git worktrees
---------------------------

Use ``scripts/dev/create_dev_worktree.ps1`` from the primary checkout instead of raw ``git worktree add`` when a contributor needs a project-local worktree under ``.worktrees/``. The helper creates ``.worktrees/<name>`` and registers the new checkout path as a Git ``safe.directory`` so future Git commands inside that worktree do not fail on Windows when the checkout is accessed by tools running under different user contexts.

.. code-block:: powershell

   .\scripts\dev\create_dev_worktree.ps1 -Name equilibrium-v3 -Branch feature/equilibrium-v3

``.worktrees/`` must stay ignored in ``.gitignore`` before using the helper. The ``safe.directory`` registration is intentionally path-specific and global to the current Windows user. Use ``-SkipSafeDirectory`` only when you plan to use per-command ``git -c safe.directory=<path> ...`` overrides instead.

Test selection rules
--------------------

Use the smallest relevant test first, then run ``uv run python scripts/dev/validate_project.py confidence`` before release, merge, or broad runtime claims. Use ``uv run python run_pytest.py --all -q`` only when you explicitly need the exhaustive historical suite.

- Python wrapper/API changes: ``uv run python run_pytest.py --api -q`` first, then ``uv run python run_pytest.py --confidence -q``.
- Native/equation changes: ``uv run python scripts/dev/build_epcsaft.py --build-only --parallel 10`` first, then ``uv run python run_pytest.py --runtime -q``, then ``uv run python run_pytest.py --confidence -q``.
- Equation traceability changes: ``uv run python scripts/docs/sync_equation_registry.py --check --strict-traceability`` then ``uv run python run_pytest.py tests/native/contracts/test_equation_registry.py -q``.
- Performance claims: ``uv run python run_pytest.py --profile -q -s`` is the quick runtime-only profile; use ``uv run python run_pytest.py --profile-full -q -s`` only for broad speed claims. Read the generated ``build/runtime_profile/*.md`` reports. Do not rely on skipped profile tests or code inspection alone.
- Plot asset changes: run the owning ``analyses/<category>/<short_id>/scripts`` coordinator or the figure-local ``analyses/<category>/<short_id>/figures/<figure_id>/scripts`` entrypoint, plus any targeted opt-in test under ``analyses/package_validation/package_plot_smokes/tests``, only when regenerating local plot outputs is explicitly part of the task.

``--profile`` is the quick runtime-only profile. ``--profile-full`` runs runtime, MIAC, and regression profiles and is the preferred evidence path for comprehensive speed reviews; use a timeout of at least 120 seconds.
- Packaging changes: ``uv run python scripts/dev/build_dist.py``. The command defaults to ``--parallel 1`` for isolated PEP 517 builds; raise it only after confirming Ceres builds are not memory-bound.

Keep generated plot assets and generated CSV workflows out of normal validation unless the task explicitly asks for them. There is no named plot validation slice; target the owning script or test file directly when plot output work is in scope.

Use ``uv run python run_pytest.py --list-slices`` when you need to inspect what each named slice runs before choosing a validation command.

For parallel sessions, leave the default repo-local temp behavior alone unless it causes contention. When running concurrent pytest lanes, set ``EPCSAFT_PYTEST_TEMP_ROOT`` to an external temp root for the extra lanes so each run gets an isolated ``pytest-temp`` child.

Runtime speed rule
------------------

For repeated runtime calls, build ``ePCSAFTMixture`` and ``ePCSAFTState`` once and reuse them inside hot loops. The quick profile report compares reused-state activity-coefficient calls against full rebuild calls and flags the ratio when rebuilds dominate runtime.

Neutral equilibrium benchmark
-----------------------------

Use the package-owned benchmark harness when the claim is specifically about neutral equilibrium throughput rather than the broader runtime profile suite.

.. code-block:: powershell

   uv run python scripts/benchmarks/benchmark_neutral_equilibrium.py --warmup 20 --repeat 100
   uv run python scripts/benchmarks/benchmark_neutral_equilibrium.py --case tp_flash --warmup 20 --repeat 200
   uv run python scripts/benchmarks/benchmark_neutral_equilibrium.py --warmup 20 --repeat 100 --json build/benchmarks/neutral_equilibrium.json
   uv run python scripts/benchmarks/benchmark_neutral_equilibrium.py --warmup 20 --repeat 100 --baseline-json build/benchmarks/neutral_equilibrium_baseline_issue43.json

The harness benchmarks these package-owned neutral cases:

- ``neutral_state``
- ``tp_flash`` for methane/ethane/propane at ``T=220 K``, ``P=1e5 Pa``, ``z=[0.1, 0.3, 0.6]``
- ``lle_seeded`` for the methanol/cyclohexane seeded LLE fixture

Each case reports a deterministic fingerprint plus medians, spread metrics,
failures, and whether a neutral fast path was used. This harness does not
require FeOs and should remain the local performance guardrail for
issue-driven neutral-equilibrium work.

Literature benchmark inventory
------------------------------

Use the package-owned literature suite inventory when the question is coverage
and benchmark scope rather than runtime speed.

.. code-block:: powershell

   uv run python scripts/benchmarks/benchmark_literature_suite.py
   uv run python scripts/benchmarks/benchmark_literature_suite.py --case figiel_2025_ssm_ds_born
   uv run python scripts/benchmarks/benchmark_literature_suite.py --json build/benchmarks/literature_suite.json

The inventory reports each issue-scope literature anchor with a classification
such as ``already_supported_with_tests`` or ``blocker_requires_followup`` plus
the owning package/test surfaces. Use it to keep benchmark claims honest and to
avoid silently treating blocked literature routes as complete.

Troubleshooting
---------------

Run ``uv run python scripts/dev/doctor.py`` whenever imports, tool paths, ``_core`` state, or generated-output tracking are unclear. It reports the active Python, git ref, uv/cmake/ninja paths, ``epcsaft`` import path, ``epcsaft._core`` path, required native symbol presence, generated artifact state, and the next recommended command.

If ``scripts/dev/build_epcsaft.py`` appears slow, run ``uv run python scripts/dev/build_epcsaft.py --status`` first. If the status output shows a stale Ninja lock and live repo-owned build processes, resolve those processes before retrying. If the status output is clean, check whether ``build/dev/CMakeCache.txt`` reports ``CMAKE_GENERATOR:INTERNAL=MinGW Makefiles``. A clean one-time switch to Ninja can materially reduce rebuild overhead on Windows systems where Ninja is already installed. Clean Ceres configure/builds can still take longer than incremental rebuilds; ``--build-only --parallel 10`` is the intended C++ edit loop after the tree is configured.
