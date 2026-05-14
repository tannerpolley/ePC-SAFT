# IntelliJ-Backed Navigation Policy for `ePC-SAFT`

## Summary
Establish a repo-supported workflow where agents proactively suggest IntelliJ-backed tooling when it is clearly the more efficient option, while keeping project opening user-approved and automation minimal. The goal is to make semantic navigation, IDE diagnostics, and safe semantic refactors part of the normal toolbox for `ePC-SAFT`, instead of leaving agents on shell and `rg` by default even when the IDE index would be better.

## Key Changes
- Update this repo’s `AGENTS.md` to state that agents should propose IntelliJ-backed tooling unprompted when it is clearly the better option for `ePC-SAFT` work.
- Clarify in the user-level `TOOL_AND_REFACTORING.md` that read-only IntelliJ semantic navigation and IDE diagnostics are encouraged when symbol meaning matters, while semantic refactor tools are allowed when they are clearly the best fit.
- Keep the approval boundary explicit:
  - Agents may suggest opening IntelliJ.
  - Agents must ask before launching or focusing IntelliJ for the repo.
  - After approval, automation stays minimal: launch or focus IntelliJ for `C:\Users\Tanner\Documents\git\ePC-SAFT`, wait for indexing readiness, then use MCP tools.
- Explicitly promote these tool categories once the project is open and indexed:
  - Semantic navigation: definitions, references, implementations, super methods, call hierarchy, type hierarchy, indexed file/class lookup.
  - IDE diagnostics: file analysis, build errors, test result diagnostics.
  - Safe semantic refactors: rename, move, and safe delete when they are clearly preferable to manual edits.
- Add concrete trigger examples in policy text so agents do not miss the opportunity:
  - tracing wrapper-to-`_core` call paths
  - checking where a public API symbol is used
  - following override or implementation chains
  - verifying whether code is truly unused before deletion
  - using semantic rename/move/delete instead of text-based edits

## Interface And Behavior Changes
- No product/API change to `intellij-index` itself is assumed.
- The behavior change is policy-level:
  - `intellij-index` is treated as an IDE-backed index over open IntelliJ projects, not as a generic filesystem indexer.
  - Agents should pass `project_path` when multiple projects may be open.
  - Agents should continue using shell and `rg` for plain text search, config files, generated files, repo-wide mechanical checks, and files outside the IntelliJ project.

## Test Plan
- Review the updated policy text for consistency between repo-local and user-level guidance.
- Validate the decision flow against representative `ePC-SAFT` tasks:
  - plain text/config lookup still stays on shell/`rg`
  - semantic definition/reference tracing leads to an IntelliJ suggestion
  - safe rename/move/delete is allowed only when clearly the best fit
  - IntelliJ opening remains suggestion-first and approval-gated
- Confirm the documented flow matches real tool behavior:
  - project must be open in IntelliJ
  - indexing readiness is checked before relying on IDE-backed tools
  - no brittle UI choreography is required beyond launch/focus and readiness waiting

## Assumptions
- The IntelliJ MCP remains available only through an active IntelliJ instance with the repo open and indexed.
- The initial rollout is intentionally limited to `C:\Users\Tanner\Documents\git\ePC-SAFT`.
- “Clearly better” means materially more efficient or safer than shell-based navigation for the task at hand, not merely available.
- Minimal automation means launch/focus plus readiness waiting only; no scripted clicking or typing inside IntelliJ by default.
