# T016 Native AD Scaffold Receipt

## Result

Added a slice-scoped native AD scaffold under `src/epcsaft/native/autodiff/`:

- `ad_scalar.h` provides the optional CppAD scalar type and templated residual primitives.
- `ad_derivative_checks.h/.cpp` provides an internal derivative-check surface for scaled, pressure-log, and reaction-log residual primitives.
- `_core._native_autodiff_derivative_checks()` exposes the check payload for native tests.
- Default builds without CppAD return `backend_unavailable` and `Backend_unavailable_used=false`.
- CppAD-enabled builds return real CppAD derivatives for the check residuals and still report `Backend_unavailable_used=false`.

This is intentionally not a whole-EOS AD conversion. It is the first native substrate needed by the supported production regression slice.

## CppAD Build Finding

The bundled FetchContent CppAD setup needed generated `cppad/configure.hpp` and `cppad_lib/temp_file.cpp` for `_core` linking. CMake now generates the config header and links the required CppAD support source when bundled CppAD is enabled.

## Validation

- `uv run python scripts/build_epcsaft.py --build-only --parallel 10`
  - pass; default dev build, CppAD disabled
- `uv run python run_pytest.py tests/native/test_cppad_eos_derivatives.py tests/native/test_cppad_reactive_speciation_derivatives.py tests/native/test_cppad_bubble_derivatives.py -q`
  - pass; `3 passed`
- `uv run ruff check tests/native/test_cppad_eos_derivatives.py tests/native/test_cppad_reactive_speciation_derivatives.py tests/native/test_cppad_bubble_derivatives.py`
  - pass
- `uv run black --check tests/native/test_cppad_eos_derivatives.py tests/native/test_cppad_reactive_speciation_derivatives.py tests/native/test_cppad_bubble_derivatives.py`
  - pass
- `cmake -S . -B build/temp/cppad-ad-check ... -DEPCSAFT_ENABLE_CPPAD=ON`
  - pass
- `cmake --build build/temp/cppad-ad-check --target _core --parallel 10`
  - pass
- `uv run python -c "import sys; sys.path.insert(0, r'build/temp/cppad-ad-check'); import _core; print(_core._native_autodiff_derivative_checks())"`
  - pass; reports `cppad_compiled=True`, `cppad_used=True`, `derivative_backend='cppad'`, `Backend_unavailable_used=False`, `max_abs_error=0.0`

## Remaining Scope

T017 still has to wire real native thermodynamic row objects/evaluator for reactive speciation and reactive electrolyte bubble rows. The AD scaffold only proves the optional CppAD substrate compiles and can provide native derivatives without Backend unavailables.

