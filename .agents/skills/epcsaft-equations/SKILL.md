---
name: epcsaft-equations
description: Use for ePC-SAFT equation-of-state theory, Helmholtz residual contributions, equation-to-code mapping, Python/pybind/native math seams, and consistency work involving docs/latex/equations.tex, generated equation docs, equation registry entries, or native contribution code.
---

# ePC-SAFT Equations

Use this skill to work on the scientific and mathematical side of this repo: ePC-SAFT theory, residual Helmholtz contribution accounting, equation ownership, and implementation tracing from equations to Python/pybind/C++ code.

## Source Of Truth

- Treat `docs/latex/equations.tex` as the equation source of truth.
- Treat `docs/equations.md` and `docs/equations_registry.yaml` as generated views used for navigation and consistency checks.
- Do not update generated equation docs directly unless the task explicitly asks for generated output repair.

## Implementation Chain

Trace behavior in this order:

1. Public Python API in `src/epcsaft/__init__.py`, `src/epcsaft/epcsaft.py`, `src/epcsaft/parameters.py`, and `src/epcsaft/regression.py`.
2. Python/pybind seam in `src/epcsaft/bindings.cpp`, which exposes the private extension `epcsaft._core`.
3. Native C++ owners in `src/epcsaft/native/**`.
4. Tests in `tests/**` and paper-validation helpers only when they are relevant to the behavior under review.

## References

Load only the reference needed for the task:

- `references/theory-overview.md`: ePC-SAFT as a Helmholtz-residual EOS and the major contribution families.
- `references/code-map.md`: where equation concepts live in the Python, pybind, and native code.
- `references/equation-workflow.md`: procedure for equation consistency, EqID ownership, and owner-agent routing.
- `references/docs-maintenance.md`: when and how to keep the skill, equation docs, and user-facing docs synchronized.

## Working Rules

- Preserve public API contracts unless the user explicitly asks for an API change.
- Keep equation edits synchronized through `docs/latex/equations.tex` and the repo's equation sync workflow.
- Update this skill when the Python/pybind/native implementation chain, equation ownership, or contribution structure changes.
- For native/equation review, prefer `native_equation_owner`; for focused Python API/test coverage, prefer `python_api_test_owner`.
- Report uncertainty with exact paths, symbols, EqIDs, and the narrowest validation command.
