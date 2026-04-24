# Codex Implementation Record: uv + pybind11 Backend Migration

This branch replaces the setuptools/Cython development loop with a `uv` + direct-CMake workflow and ports the Python-facing runtime to pure Python wrappers backed by `epcsaft._core`.

## Current State

- Package layout remains `src/epcsaft`.
- Native equation code remains under `src/epcsaft/native`.
- Cython files and setuptools build metadata are removed from the runtime/build path.
- `pyproject.toml` now uses `scikit-build-core` for packaging and `tool.uv.package = false` so the default `uv sync --no-install-project` loop does not invoke PEP 517 hooks.
- `uv.toml` routes uv's cache into `build/uv-cache` so sandboxed Windows runs do not touch `%LOCALAPPDATA%\uv\cache`.
- `CMakeLists.txt` builds `epcsaft._core` with pybind11 and copies it into `src/epcsaft` when `EPCSAFT_DEV_INPLACE=ON`.

## Default Commands

```powershell
uv sync --no-install-project
uv run python scripts\build_epcsaft.py --clean
uv run python scripts\codex_doctor.py
uv run python run_pytest.py tests\test_runtime.py -q
```

Use `uv build` only for package artifacts and CI-style checks.
