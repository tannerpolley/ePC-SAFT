# Subagent Briefs - Issue #120 Cleanup

Use these briefs when activating GoalBuddy tasks. Agents are not alone in the codebase: do not revert unrelated edits, do not overlap writes with another worker, and coordinate native rebuilds through the main PM thread.

## Required Shared Context

- Read `docs/roadmaps/FULL_ROADMAP.md` and GitHub issue #120 before changing files.
- Treat issue #120 and comment `4456462579` as the execution sequence.
- Keep the cleanup conservative: relocation, de-duplication, guardrail hardening, and public-surface clarification only.
- Do not change thermodynamic equations, solver algorithms, regression math, benchmark targets, tolerances, parameter values, or scientific outputs.
- Avoid exact banned backend and derivative-status tokens in committed text.

## PM

- Owns branch hygiene, `state.yaml`, commits, final cleanup hook, and any push or PR decision.
- Runs or delegates baseline and final validation.
- Coordinates `_core` rebuilds and any clean native build so tests and agents do not hold the extension.

## Scout

- Best for Phase 0 inventories and import-surface mapping.
- Must produce file-path evidence and exact commands, not generic summaries.
- Should map risks for benchmark relocation, diagnostics relocation, native scaffolding, staged/optional modules, broad catches, and public API exposure.

## Judge

- Gates Phase 0 before implementation.
- Reviews ambiguous compatibility-shim decisions and native-regression Path A versus Path B.
- Final audit must answer whether every #120 deliverable exists, validation passed, and no scientific behavior changed.

## Worker

- Executes one queued phase at a time with the `allowed_files`, `verify`, and `stop_if` on that task.
- Stops instead of widening scope when imports, scientific behavior, native build ownership, or public API compatibility are ambiguous.
- Leaves receipts with changed files, commands, summary, and any remaining blockers.

## ePC-SAFT Owner Agents

- `build_packaging_owner`: use for build scripts, CMake, pybind, validation workflow, distribution or docs validation risk.
- `native_equation_owner`: use as read-only reviewer if a cleanup could accidentally touch EOS/property/Helmholtz contribution semantics.
- `native_solver_backend_owner`: use for native regression scaffolding, Ceres, CppAD, solver/result contracts, and native backend risk.
- `python_api_test_owner`: use for public API, compatibility shim, runtime import, pytest, and wrapper behavior risk.
- `command_runner`: use for validation-only runs and exact failure capture.

## JetBrains Use

Use JetBrains semantic tools when symbol ownership, references, call hierarchy, safe delete, or rename safety materially improves correctness. If the IDE MCP is unavailable or indexing is not ready, fall back to `rg`, repo tests, and focused file inspection.
