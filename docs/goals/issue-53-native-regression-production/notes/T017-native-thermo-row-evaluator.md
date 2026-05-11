# T017 Native Thermodynamic Row Evaluator Receipt

## Result

Added the first native thermodynamic regression row evaluator:

- `src/epcsaft/native/regression/thermo_regression.h`
- `src/epcsaft/native/regression/thermo_regression.cpp`
- `_core._evaluate_native_thermo_regression_rows(...)`
- `epcsaft.evaluate_native_thermo_regression_rows(...)`

Supported row modes in this slice:

- `reactive_speciation`: calls `chemical_equilibrium_native(...)` in C++ and packs fixed-shape `speciation`, `reaction`, and `activity` residual targets.
- `reactive_electrolyte_bubble`: calls `electrolyte_bubble_pressure_native(...)` in C++ and packs fixed-shape pressure or vapor-composition targets.

Unsupported row modes receive fixed-shape penalty residuals instead of falling back to Python.

## Important Limit

This task adds native thermodynamic row evaluation, not Ceres parameter iteration and not implicit sensitivities. `reactive_speciation` can still report a finite-difference Jacobian internally for nonideal standard states until T018/T020 harden derivative policy. The T017 speciation test uses an ideal mole-fraction standard state so the native row evaluator reports `analytic`.

## Validation

- `uv run python scripts/build_epcsaft.py --build-only --parallel 10`
  - pass
- `uv run python run_pytest.py tests/native/test_native_ceres_thermodynamic_regression.py tests/native/test_native_ceres_reactive_pressure_speciation.py -q`
  - pass; `3 passed`
- `uv run ruff check src/epcsaft/native_regression.py src/epcsaft/__init__.py tests/native/test_native_ceres_thermodynamic_regression.py tests/native/test_native_ceres_reactive_pressure_speciation.py`
  - pass
- `uv run black --check src/epcsaft/native_regression.py src/epcsaft/__init__.py tests/native/test_native_ceres_thermodynamic_regression.py tests/native/test_native_ceres_reactive_pressure_speciation.py`
  - pass
