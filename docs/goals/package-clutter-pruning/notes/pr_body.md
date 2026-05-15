# Summary

Closes issue #120 by pruning package clutter before the next vertical
implementation lane.

- Moves benchmark execution logic out of the installable package while keeping
  compatibility import shims.
- Moves equilibrium confidence and thermodynamic diagnostic implementations into
  workflow validation helpers while keeping package shims.
- Removes the unused native regression scaffold and CMake hooks.
- Moves the IPOPT implementation behind a private optional-backend package while
  keeping `epcsaft.ipopt_backend` as the public compatibility alias.
- Narrows broad exception handling, rewrites ambiguous package-facing wording,
  and documents the public API surface plus future module split work.
- Repairs workflow validation targets after the diagnostics test move.

# Validation

- `uv run python scripts/dev/validate_project.py quick` - pass, 31 tests.
- `uv run python scripts/dev/validate_project.py docs` - pass, Sphinx build.
- `uv run python run_pytest.py tests/api/package -q` - pass, 6 tests.
- `uv run python run_pytest.py tests/api/runtime -q` - pass, 43 tests.
- `uv run python run_pytest.py tests/workflows/repo -q` - pass, 53 tests.
- `git diff --check` - pass.
- Tracked artifact scan - pass with justified `tests/workflows/build/*` matches.

# Scope

This PR is organization cleanup. It does not intentionally change equation
forms, solver algorithms, regression objectives, benchmark target values, or
public runtime semantics. Compatibility imports remain where the current tests
and docs exercise them.
