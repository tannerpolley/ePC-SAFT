# T003 Build And Dependency Plumbing Receipt

Date: 2026-05-11

## Result

Done.

## Changes

- Broadened `EPCSAFT_ENABLE_CERES` wording to package-wide native Ceres support.
- Added `EPCSAFT_ENABLE_CPPAD` and `EPCSAFT_USE_SYSTEM_CPPAD`.
- Added CppAD FetchContent include plumbing for the default non-system path.
- Added `EPCSAFT_CERES_AVAILABLE` and `EPCSAFT_CPPAD_AVAILABLE` CMake cache keys for downstream capability reporting.
- Added `cmake/FindEigen3.cmake` so Ceres can consume the repo's `includeigen` Eigen headers during FetchContent configuration.
- Set Ceres FetchContent to `MINIGLOG=ON` so the default local source-build path does not require a system `glog`.
- Added `src/epcsaft/native/regression/*.cpp` and the matching include directory to the native build.
- Added `scripts/build_epcsaft.py` flags:
  - `--enable-ceres`
  - `--use-system-ceres`
  - `--enable-cppad`
  - `--use-system-cppad`
- Added doctor/runtime native dependency reporting for Ceres and CppAD.
- Updated runtime capability truth so the current Python-batched reactive regression path no longer reports itself as Issue #53 native production.
- Updated focused tests and package architecture/diagnostics docs for the new dependency controls.

## Validation

- `uv run python -m py_compile src/epcsaft/runtime.py scripts/build_epcsaft.py scripts/doctor.py`
- `uv run python run_pytest.py tests/api/test_runtime.py tests/workflows/test_build_epcsaft_script.py tests/workflows/test_workflow_entrypoints.py -q`
  - 51 passed.
- `uv run python run_pytest.py tests/api/test_runtime.py tests/workflows/test_build_epcsaft_script.py tests/workflows/test_workflow_entrypoints.py tests/workflows/test_run_pytest.py::test_doctor_recommends_ninja_migration_for_mingw_build_tree tests/workflows/test_run_pytest.py::test_doctor_does_not_recommend_ninja_migration_when_already_ninja -q`
  - 53 passed.
- `uv run python scripts/build_epcsaft.py --build-only --parallel 10`
  - Passed after CMake reconfigure; final full native rebuild completed in about 34 seconds.
- `uv run python scripts/doctor.py`
  - Passed; Ceres/CppAD reported disabled and available false in the default build.
- `uv run ruff check src/epcsaft/runtime.py scripts/build_epcsaft.py scripts/doctor.py tests/api/test_runtime.py tests/workflows/test_build_epcsaft_script.py tests/workflows/test_workflow_entrypoints.py tests/workflows/test_run_pytest.py`
  - Passed.
- `uv run black --check src/epcsaft/runtime.py scripts/build_epcsaft.py scripts/doctor.py tests/api/test_runtime.py tests/workflows/test_build_epcsaft_script.py tests/workflows/test_workflow_entrypoints.py tests/workflows/test_run_pytest.py`
  - Passed.
- Isolated configure probe:
  - `cmake -S . -B build/temp/cppad-config-check ... -DEPCSAFT_ENABLE_CPPAD=ON`
  - Passed.
- Isolated configure probe:
  - `cmake -S . -B build/temp/ceres-config-check ... -DEPCSAFT_ENABLE_CERES=ON`
  - Passed after Eigen shim and `MINIGLOG=ON`.
- Isolated configure probe:
  - `cmake -S . -B build/temp/ceres-cppad-config-check ... -DEPCSAFT_ENABLE_CERES=ON -DEPCSAFT_ENABLE_CPPAD=ON`
  - Passed.

## Notes

- The Ceres/CppAD switches are build/capability plumbing only. No production regression solver/backend behavior is implemented by T003.
- A broad `ruff check src tests scripts docs/pages/*.rst` was intentionally not used as evidence because it linted `.rst` as Python and exposed unrelated pre-existing script lint outside this tranche.
