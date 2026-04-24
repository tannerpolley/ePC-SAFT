# Equation Workflow

## Equation Consistency Procedure

1. Start from `docs/latex/equations.tex`; it is the source of truth for documented equations and EqID metadata.
2. Use `docs/equations.md` for quick human-readable navigation.
3. Use `docs/equations_registry.yaml` for generated EqID metadata, labels, descriptions, and implementation-owner hints.
4. Search native code for the EqID, symbol name, or relevant contribution function.
5. If no native owner exists, report it as missing ownership rather than guessing.
6. Identify the narrowest Python/API test or validation command that should cover the behavior.

## Do Not Treat Generated Docs As Source

Do not manually edit `docs/equations.md` or `docs/equations_registry.yaml` as the primary fix for equation drift. Update `docs/latex/equations.tex` or the code owner comments/workflow that generates those files, then regenerate with the repo equation sync script when the task calls for it.

## Owner-Agent Routing

- Use `native_equation_owner` for read-only native C++ and equation review, especially when tracing `docs/latex/equations.tex` to `src/epcsaft/native/**`.
- Use `python_api_test_owner` to identify or update focused Python API tests and `run_pytest.py` validation for equation-facing behavior.
- Use `build_packaging_owner` only if the equation/code task affects pybind/CMake extension wiring, native source lists, wheels, sdists, or uv/CMake workflow scripts.

## Reporting Format

For equation or math findings, report:

- EqID or equation label when available.
- Source path and line or section in `docs/latex/equations.tex`.
- Native owner path and function/symbol, or state that no owner was found.
- Affected public Python method if applicable.
- Likely failure mode.
- Narrow validation command.
