# T019 Native Ceres Thermodynamic Fit Receipt

## Result

Implemented an honest first native thermodynamic Ceres fit slice:

- Added `fit_native_thermo_regression(...)` in native C++ and exposed it through pybind/Python.
- C++ now applies serialized parameter vectors inside the native hot loop.
- Supported production fit slice is currently reactive-speciation rows with reaction log-equilibrium-constant parameters, ideal mole-fraction standard states, and speciation targets.
- Ceres owns parameter iteration when `EPCSAFT_ENABLE_CERES=ON`.
- Jacobians for that slice use analytic/implicit sensitivities from the converged native speciation residual Jacobian. No Python objective loop and no backend-unavailable derivative path is used.
- Unsupported row/parameter derivative combinations return canonical `backend_unavailable`.

## Scope And Limitations

This is not the full issue #53 completion slice yet. It does not yet provide Ceres derivatives for reactive electrolyte bubble pressure rows, Born/SSM+DS parameters, `k_ij`, or activity-coupled speciation standard states. Those combinations are intentionally reported as `backend_unavailable` until analytic/CppAD/implicit derivatives exist.

## Evidence

- Default native build passed:
  - `uv run python scripts/build_epcsaft.py --build-only --parallel 10`
- Focused default tests passed:
  - `uv run python run_pytest.py tests/native/test_cppad_reactive_speciation_derivatives.py tests/native/test_cppad_bubble_derivatives.py tests/native/test_native_ceres_thermodynamic_regression.py tests/native/test_native_ceres_reactive_pressure_speciation.py tests/native/test_native_regression_types.py -q`
  - Result: `14 passed in 5.57s`
- Lint/format passed for touched Python surfaces:
  - `uv run ruff check src/epcsaft/native_regression.py src/epcsaft/__init__.py tests/native/test_cppad_reactive_speciation_derivatives.py tests/native/test_native_ceres_thermodynamic_regression.py tests/native/test_native_regression_types.py`
  - `uv run black --check src/epcsaft/native_regression.py src/epcsaft/__init__.py tests/native/test_cppad_reactive_speciation_derivatives.py tests/native/test_native_ceres_thermodynamic_regression.py tests/native/test_native_regression_types.py`
- Ceres-enabled temp build passed after adding a MinGW-only Ceres internal compile workaround:
  - `cmake -S . -B build/temp/ceres-fit-check -G Ninja -Dpybind11_DIR=... -DPython_EXECUTABLE=... -DEPCSAFT_ENABLE_CERES=ON -DEPCSAFT_DEV_INPLACE=OFF`
  - `cmake --build build/temp/ceres-fit-check --target _core --parallel 10`
- Direct Ceres-enabled smoke proved objective decrease:
  - `optimizer_backend=ceres`
  - `derivative_backend=analytic_implicit`
  - `initial_cost=0.05174377188057001`
  - `final_cost=1.1771098705162825e-19`
  - message includes `native_hot_loop=true; python_objective_used=false; Backend_unavailable_used=false`

## Files

- `CMakeLists.txt`
- `src/epcsaft/native/epcsaft_electrolyte.h`
- `src/epcsaft/native/regression/regression_types.h`
- `src/epcsaft/native/regression/regression_types.cpp`
- `src/epcsaft/native/regression/thermo_regression.h`
- `src/epcsaft/native/regression/thermo_regression.cpp`
- `src/epcsaft/native/regression/implicit_sensitivity.cpp`
- `src/epcsaft/bindings.cpp`
- `src/epcsaft/native_regression.py`
- `src/epcsaft/__init__.py`
- `src/epcsaft/runtime.py`
- `tests/native/test_cppad_reactive_speciation_derivatives.py`
- `tests/native/test_native_ceres_thermodynamic_regression.py`
- `tests/native/test_native_regression_types.py`

