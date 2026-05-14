# Local Codex Instructions (Machine-Local, Do Not Commit)

## Startup Reads
- Read `C:\Users\Tanner\Documents\git\ePC-SAFT\docs\.codex-journal\user_preferences.md`.
- Read `C:\Users\Tanner\Documents\git\ePC-SAFT\docs\roadmaps\FULL_ROADMAP.md` before planning, coding, reviewing, or merging any new repo task. Treat it as the authoritative package-context and completion-standard roadmap unless the user explicitly narrows scope for the current task.
- Use the user-level `chemical-engineer` skill for ePC-SAFT thermodynamics, phase/chemical equilibrium, equation tracing, Python/pybind/native seams, and scientific validation work. Load only the specific user-level chemical-engineer references needed for the task.
- Use the user-level `jetbrains` skill when IntelliJ-backed semantic navigation, IDE diagnostics, or safe semantic refactors would materially improve correctness or speed for `ePC-SAFT` work.

## Skill Use
- Always use the user-level `grill-me` skill when the user asks to be grilled, asks to stress-test a plan or design, or explicitly prompts for `grill-me`.
- If the current thread has not loaded `grill-me`, say that briefly and follow the installed skill's behavior manually: interview the user one question at a time, walk decision branches in dependency order, and include a recommended answer with each question.

## Memory Policy
- `user_preferences.md`: durable user preferences only.
- Do not write memory for routine Q&A, small one-off edits, placeholder notes, or facts already recorded.
- Keep memory entries concise, date-stamped, deduplicated, and free of secrets.

## Repo Workflow
- Current backend: uv-managed Python, direct CMake dev builds, pybind11 `_core`, native C++ equations, pure Python public API wrappers.
- Public repo tools intentionally use developer-neutral names. Do not add new tracked files, tests, scripts, or docs with Codex-specific names unless they are local-only agent instructions.
- Best new-agent workflow: `uv sync --no-install-project`, then `uv run python scripts/build_epcsaft.py`, then `uv run python scripts/doctor.py`, then `uv run python scripts/validate_project.py quick`.
- Preferred setup: `uv sync --no-install-project`.
- Preferred native build: `uv run python scripts/build_epcsaft.py`; use `uv run python scripts/build_epcsaft.py --build-only --parallel 10` only after `build/dev` is already configured.
- Preferred tests: `uv run python run_pytest.py <focused-test-targets> -q`. Prefer this wrapper over direct `pytest` for Codex/Windows work because it sets the source path and manages per-run pytest temp state.
- Preferred high-level validation: `uv run python scripts/validate_project.py quick` for normal fast validation; `uv run python scripts/validate_project.py confidence` before handoff when native runtime confidence matters; `uv run python scripts/validate_project.py docs` for Sphinx.
- Preferred doctor: `uv run python scripts/doctor.py`.
- Preferred distribution check: `uv run python scripts/build_dist.py` or `uv build` when specifically testing the packaging boundary.
- Use `.codex\environments\environment.toml` actions when available; they are aligned to the uv/CMake/pybind workflow.
- Prefer normal native builds. Treat `scripts/build_epcsaft.py --clean`, `Repair Native Build (Clean)`, and `Clean Build Artifacts` as coordinated repair actions, not routine validation.
- Coordinate native rebuilds through the main thread when multiple agents or processes may be active. Prefer one `_core` builder at a time.
- Do not run `_core`-deleting clean/repair actions in parallel with tests, Python REPLs, IDE run configurations, or sub-agents that may import `epcsaft._core`.
- If `_core*.pyd` is locked, stop the importing Python/test/IDE/Codex process, rerun the normal native build, then run doctor.
- Do not reintroduce Conda, Cython, setuptools editable installs, `setup.py build_ext`, or `tests/test_cython.py` as the normal workflow.
- Treat `build/dev`, `build/temp`, `build/pytest-temp`, `build/runtime_profile`, and `build/uv-cache` as shared disposable generated state; do not track build artifacts, profile reports, or `_core*.pyd`.

## Sandbox Notes
- Do not use plot gallery server/index workflows in this repo. Plot assets remain source-owned artifacts; any gallery app should discover them externally.
- For long-running localhost servers, prefer foreground actions. Do not expect a sandboxed `Start-Process` child to keep running after the shell command returns; request an escalated background start only when the agent must own a persistent server.
- For live process checks in the sandbox, prefer `Get-Process`, `Get-NetTCPConnection`, and `Invoke-WebRequest` with short explicit timeouts. Avoid WMI/CIM command-line inspection (`Get-CimInstance Win32_Process`) unless escalated; it can be denied by the Windows restricted token.
- For closed-port checks, prefer `Invoke-WebRequest -TimeoutSec` or a short .NET TCP probe over `Test-NetConnection`, which can outlive short command timeouts and emit noisy warnings.

## Git Sandbox Rules
- `workspace-write` keeps protected paths such as `.git/` read-only. Do not first try `.git`-writing commands in the sandbox and then retry after failure.
- Run normal local git write commands with sandbox escalation on the first attempt: `git add`, `git commit`, `git fetch`, `git restore --staged`, `git rm --cached`, `git switch -c <branch>`, and `git checkout -b <branch>`.
- User-level prefix rules in `C:\Users\Tanner\.codex\rules\default.rules` already allow common local git workflow commands, but rules do not make `.git/` writable if the command still runs inside the sandbox. Use the narrow matching git escalation approval; do not request broad shell or danger-full-access approval for routine git staging/commits.
- If a git write still reports `.git/index.lock: Permission denied`, check `Test-Path -LiteralPath '.git\index.lock'`. If the file is absent, treat it as sandbox/approval friction, not a stale lock.
- Keep `git push`, `git pull`, `git merge`, `git rebase`, non-branch-creation `git checkout`/`git switch`, `git stash`, and tags explicit/prompted. Never run destructive commands such as `git reset --hard` or `git clean` unless the user explicitly requests them.

## Documentation And Skill Maintenance
- Keep user-facing docs updated when public API behavior, install/build/test workflow, or package layout changes.
- Keep `.codex\environments\*.toml` and `.codex\environments\README.md` aligned with current repo commands whenever build, test, docs, plot, package, or script names change. Remove stale actions immediately when workflows are deleted, especially plot/gallery/server actions and old test paths.
- After workflow changes, sanity-check local environment configs with TOML parsing and a stale-command search before handoff; at minimum scan for removed scripts, removed flags, old test paths, and obsolete docs/gallery references.
- Keep the user-level `chemical-engineer` skill current when equation ownership, native contribution structure, Python/native seams, or EoS theory guidance changes.
- Do not duplicate detailed equation theory in AGENTS.md; put durable science/code navigation details in the user-level `chemical-engineer` skill references.
- Prefer concise docs that explain the current workflow; remove obsolete workaround history unless it is still needed to avoid a recurring pitfall.

## Plot Math Notation
- Any new or regenerated project plot should use proper math notation wherever a standard thermodynamic symbol exists. Prefer Matplotlib/Plotly LaTeX-style labels over plain text such as `h-res`, `A-res`, `rho`, `fugacity coefficient`, or `activity coefficient`.
- Use labels such as `r"$A^{res}$"`, `r"$h^{res}$"`, `r"$g^{res}$"`, `r"$s^{res}$"`, `r"$\rho$"`, `r"$Z$"`, `r"$\phi_i$"` or `r"$\varphi_i$"` for fugacity coefficient, `r"$\gamma_i$"` / `r"$\gamma_{\pm}$"` for activity coefficients, `r"$\mu_i^{res}$"` for residual chemical potential, and `r"$\epsilon_r$"` for relative permittivity when applicable.
- Keep units outside the math expression when possible, for example `r"$\rho$ / mol m$^{-3}$"` or `r"$P$ / Pa"`. If a paper figure uses a specific published notation, follow the paper unless it conflicts with package-wide clarity.
- Apply the same convention to titles, axes, legends, colorbars, hover labels, and CSV-backed interactive plot labels when those labels are generated by repo code.

## Sub-Agent Policy
- Standing user preference: use repo sub-agents for non-trivial ePC-SAFT work when the task has clear non-blocking slices.
- Main agent acts as orchestrator for cross-cutting work: keep final decisions, integration, sandbox escalation decisions, and immediately blocking implementation on the main thread.
- Main agent owns `_core` rebuild, clean, and repair coordination. Sub-agents may inspect files and run focused tests, but should not rebuild `_core` or run `_core`-deleting repair/cleanup commands unless explicitly assigned and coordinated.
- Prefer Spark owner agents for sidecar exploration, review, validation, and bounded implementation because they are cheaper and focused.
- Do not use sub-agents for simple Q&A, tiny one-file edits, mechanical text edits, or when the next step is blocked on the delegated result.
- When using sub-agents, assign explicit ownership by file/module area and avoid overlapping writes.
- Use multiple sub-agents early only when their work can run in parallel without blocking the main thread.

## Repo Owner Agents
- `build_packaging_owner`: package/build workflow owner for uv, CMake, pybind, scikit-build, wheels/sdists, and build scripts.
- `native_equation_owner`: read-only native math/equation reviewer for C++ equation logic, contribution accounting, density/phase closure, and equation-doc consistency.
- `python_api_test_owner`: Python API/test owner for public wrappers, pytest workflow, focused tests, and validation against `epcsaft._core`.
- `command_runner`: validation-only runner for non-blocking doctor/build/test/smoke commands. It should not edit files.

## Routing Playbooks
- Build/package changes: delegate build/package review or bounded edits to `build_packaging_owner`; delegate independent command validation to `command_runner`.
- Native/equation changes: delegate correctness review to `native_equation_owner`; delegate focused API/test coverage to `python_api_test_owner`.
- Python API/runtime changes: delegate tests and API risk checks to `python_api_test_owner`; involve `native_equation_owner` if behavior crosses into native equation logic.
- Cross-layer changes: use `native_equation_owner` for correctness risk, `python_api_test_owner` for API/test evidence, `build_packaging_owner` for build/package impact, and `command_runner` for parallel validation.
- Branch review/change audit: use owner agents by concern area, then have the main thread synthesize findings and decide edits.

## IntelliJ-Backed Tooling
- Load and follow the user-level `jetbrains` skill for semantic code work, public API usage tracing, call/type hierarchy work, safe semantic refactors, or post-edit IDE diagnostics.
- Propose IntelliJ-backed tooling unprompted when it is clearly the more efficient or safer option for `ePC-SAFT` work.
- Ask before launching or focusing IntelliJ for `C:\Users\Tanner\Documents\git\ePC-SAFT`.
- After approval, keep automation minimal: launch or focus IntelliJ for this repo, wait for indexing readiness, then use MCP tools.
- Treat `intellij-index` as an IDE-backed index over open IntelliJ projects, not as a generic filesystem indexer.
- Call `ide_index_status` at most once near the start of semantic work, pass `project_path`, and fall back without repeated retries if MCP is unavailable or indexing is not ready.
- Prefer IntelliJ-backed semantic navigation and diagnostics when symbol meaning matters:
  - definitions, references, implementations, super methods, call hierarchy, type hierarchy, indexed file/class lookup
  - file analysis, build errors, and test result diagnostics
- Safe semantic refactors are encouraged when they are clearly preferable to manual edits:
  - `ide_refactor_rename`
  - `ide_move_file`
  - `ide_refactor_safe_delete`
- Continue to prefer shell and `rg` for plain text search, config files, generated files, repo-wide mechanical checks, and files outside the IntelliJ project.
- Pass `project_path` when multiple IntelliJ projects may be open.
- Do not manually edit `.idea/workspace.xml` or `.idea/*.iml` through shell write commands. Use `scripts/dev/configure_jetbrains_project.py` for deterministic JetBrains metadata cleanup when that path exists.
- Concrete trigger examples for preferring IntelliJ-backed tooling:
  - tracing wrapper-to-`_core` call paths
  - checking where a public API symbol is used
  - following override or implementation chains
  - verifying whether code is truly unused before deletion
  - using semantic rename, move, or safe delete instead of text-based edits

## Git Commit Policy
- For tracked repo changes, finish with local git commits unless explicitly told not to; do not push unless explicitly requested.
- `docs/latex` is normal tracked content in this repo, not a Git submodule.
- When LaTeX files should be sent to Overleaf, use `scripts\setup_latex_mirror.ps1` once to create or validate `C:\Users\Tanner\Documents\git\ePC-SAFT-LaTeX`, then run `scripts\sync_latex_mirror.ps1` from the main repo. That external mirror owns the Overleaf Git remote and performs the mirror commit/push.
- Do not run `git -C docs/latex ...` as an Overleaf submodule workflow; `docs/latex` has no nested Git checkout.
- Before handoff, review `git status --short` and do not leave modified tracked files uncommitted.
- If the worktree is already dirty, include all modified tracked files in the commit plan unless the user says otherwise.
- Split commits by logical concern when tracked changes cover different areas.
- Do not force ignored local-only files such as `AGENTS.md`, `.codex/**`, or `docs/.codex-journal/**` into git commits.
- When staging or committing, use the first-attempt escalation rule from `Git Sandbox Rules` so the workflow does not waste a failed sandboxed `.git` write.

## Scope
- These instructions and `.codex/**` files are machine-local for this repo and should remain untracked unless the user explicitly asks otherwise.

## ePC-SAFT Cross-Repo Integration
- This repo is the upstream source of truth for the `epcsaft` package. It owns public APIs, native kernels, solver/result contracts, package tests, docs, and release behavior.
- Downstream repos such as `Lithium_Extraction`, `MEA-Thermodynamics`, and `MEA-Absorption-Column` own their case-study data, project-specific scientific thresholds, process models, and reports.
- Use the user-level `epcsaft-cross-repo` skill for upstream/downstream contract work, downstream-consumer policy, package feedback loops, and cross-repo handoffs.
- Do not absorb full downstream studies into normal upstream validation. Reduce repeated downstream friction to compact public-API reproductions before adding upstream tests or fixtures.
- Keep `capabilities()` honest. Do not broaden package capability claims unless package validation proves the supported path.
