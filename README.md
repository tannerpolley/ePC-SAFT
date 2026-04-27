# ePC-SAFT

`epcsaft` is a Python package for the ePC-SAFT equation of state with association and electrolyte terms. The public API is Python; the performance-critical thermodynamic runtime is native C++ exposed through a private pybind11 module, `epcsaft._core`.

## Install

Most users should install from PyPI:

```bash
pip install epcsaft
```

The package includes a compiled C++ extension. Wheels are preferred for end users; source builds require a working C++ toolchain, CMake, and Ninja.

## Development

This repository uses `uv` for Python environment management and direct CMake for the local native build loop.

`uv.toml` routes uv's cache into `build/uv-cache` so Codex/Windows sandbox runs do not touch `%LOCALAPPDATA%\uv\cache`.

```powershell
uv sync --no-install-project
uv run python scripts\build_epcsaft.py
uv run python scripts\codex_doctor.py
uv run python run_pytest.py --confidence -q
```

Direct pytest also works, for example `uv run python -m pytest tests\test_runtime.py -q`, but `uv run python run_pytest.py ...` is preferred for Codex and Windows runs because it manages pytest temporary directories more predictably. Set `EPCSAFT_PYTEST_TEMP_ROOT` when you want the wrapper to use an opt-in external pytest temp root instead of its default repo-local generated temp area.

The default new-agent validation sequence is sync, normal native build, doctor, then `uv run python run_pytest.py --confidence -q`.

For the standard Codex validation loops:

```powershell
uv run python run_pytest.py --generic -q
uv run python run_pytest.py --confidence -q
```

`--generic` runs the fast core runtime, parameter-template, equation-registry, and regression API slice. `--confidence` is the default runtime-confidence check; it runs that same slice plus the native runtime contract tests. To keep pytest temp files outside the repo for an opt-in run, set `EPCSAFT_PYTEST_TEMP_ROOT`, for example:

```powershell
$env:EPCSAFT_PYTEST_TEMP_ROOT = Join-Path $env:TEMP 'epcsaft-pytest'
uv run python run_pytest.py --confidence -q
```

Use `uv run python scripts\build_epcsaft.py --clean` only as a repair step for stale CMake state or stale/locked `_core` artifacts. If a `_core*.pyd` is locked, stop the importing Python/test/IDE process before running the clean repair.

`CMakePresets.json` is optional Windows MinGW convenience for IDEs and manual CMake use. The canonical local native build remains `uv run python scripts\build_epcsaft.py`.

Build distributable artifacts at the packaging boundary:

```powershell
uv build
```

## Example

```python
import numpy as np
from epcsaft import ePCSAFTMixture

mixture = ePCSAFTMixture.from_params(
    {"m": np.asarray([2.8149]), "s": np.asarray([3.7169]), "e": np.asarray([285.69])},
    species=["Toluene"],
)
state = mixture.state(T=320.0, x=np.asarray([1.0]), P=101325.0)
print(state.density())
print(state.ares())
print(state.ares(return_contribution_terms=True)["terms"]["hc"])
```

## Package Layout

- `src/epcsaft/`: Python package, pybind11 binding source, and native C++ equation code
- `data/epcsaft_parameters/`: source-checkout example parameter datasets for inspection, comparison, and tests
- `docs/`: user documentation and reference material
- `scripts/`: uv/CMake build, doctor, and validation helpers

## Public API

The main entry points are `ePCSAFTMixture`, `ePCSAFTState`, `create_parameter_template`, `fit_pure_neutral(...)`, and the structured result objects returned by solver-style methods.

`create_parameter_template(...)` creates a blank dataset folder with the expected `pure/`, `mixed/`, and `user_options.json` layout. After you fill in the files, `ePCSAFTMixture.from_dataset(...)` can load the folder path you created yourself.

Phase 1 of the regression workflow is intentionally narrow: `fit_pure_neutral(...)` fits only nonassociating neutral-component `m`, `s`, and `e` against liquid-density and vapor-pressure data. Use `write_fit_result(...)` when you want to persist a fit back into a user-owned dataset folder. Ion and binary regression are deferred for now.

## Documentation

- [Start here](docs/pages/README.rst)
- [Getting started](docs/pages/getting_started.rst)
- [Create your own parameter folder](docs/pages/user_parameter_templates.rst)
- [Parameter regression guide](docs/pages/parameter_regression.rst)
- [User options reference](docs/pages/user_options.rst)
- [Task-based package guide](docs/pages/package_guide.rst)
- [API reference](docs/pages/api_reference.rst)

## License

GNU General Public License v3.0
