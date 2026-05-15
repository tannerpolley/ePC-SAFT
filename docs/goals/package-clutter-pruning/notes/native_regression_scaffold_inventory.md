# Native Regression Scaffold Inventory

Phase 4 source scope from issue #120:

- `src/epcsaft/native/regression/ceres_backend.cpp`
- `src/epcsaft/native/regression/ceres_backend.h`
- `src/epcsaft/native/regression/parameter_map.h`
- `src/epcsaft/native/regression/regression_problem.cpp`
- `src/epcsaft/native/regression/regression_problem.h`
- `src/epcsaft/native/regression/target_rows.h`
- `src/epcsaft/native/regression/thermo_evaluator.h`

## Classification

The folder contained tiny contract/placeholder declarations rather than the
compiled regression implementation. The real native regression behavior remains
in `src/epcsaft/native/epcsaft_regression.cpp`, with public native structs and
function declarations in `src/epcsaft/native/epcsaft_electrolyte.h`.

| Removed file | Prior role | Runtime dependency |
| --- | --- | --- |
| `ceres_backend.cpp` / `ceres_backend.h` | One compile-flag status helper for `_native_ceres_smoke`. | Replaced by a direct compile-flag check in `src/epcsaft/bindings.cpp`. |
| `regression_problem.cpp` / `regression_problem.h` | Unreferenced contract struct and helper. | None. |
| `parameter_map.h` | Unreferenced three-name parameter map. | None. |
| `target_rows.h` | Unreferenced target-row enum. | None. |
| `thermo_evaluator.h` | Unreferenced contract struct. | None. |

## Build Ownership

`CMakeLists.txt` no longer glob-builds `src/epcsaft/native/regression/*.cpp` or
adds that folder to native include paths. Native regression compilation remains
owned by the top-level native source glob through `src/epcsaft/native/epcsaft_regression.cpp`.

## Acceptance Notes

No regression objective math, solver algorithm, target definitions, parameter
values, or public scientific outputs were changed. This phase removed scaffold
that could be mistaken for a real native-regression architecture while keeping
the existing Ceres smoke metadata behavior.
