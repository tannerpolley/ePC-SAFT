# GitHub Issue #96 Scope

- Issue URL: https://github.com/tannerpolley/ePC-SAFT/issues/96
- Title: Downstream integration smoke tests
- State when read: OPEN
- Branch: codex/downstream-integration-smokes
- Dependencies: F, G, I, J, K

## Purpose

Prove downstream projects can use generic package APIs without private workaround code.

## In Scope

- MEA-Thermodynamics smoke
- Lithium_Extraction smoke
- MEA-Absorption-Column smoke
- Generic problem construction
- Generic outputs consumed downstream
- No copied EOS implementation

## Out Of Scope

- Downstream-application-specific public APIs in `epcsaft`
- Downstream metrics computed inside `epcsaft` package APIs

## Required Policy

- No finite difference
- CppAD for explicit algebraic derivatives
- Analytic formulas where exact and validated
- Implicit sensitivities for solved states
- No `backend_unavailable` for required workflows
- `backend_unavailable` only for explicitly out-of-scope workflows
- No application-specific public APIs

## Validation Named By Issue

- `uv run python scripts/validate_project.py quick`
- `uv run python scripts/validate_project.py docs`
