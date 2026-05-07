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

`uv.toml` routes uv's cache into `build/uv-cache` so local and sandboxed runs do not touch the user-level uv cache.
The current development and CI smoke baseline is Python 3.13, while package metadata still allows Python `>=3.9`.

```powershell
uv sync --no-install-project
uv run python scripts\build_epcsaft.py
uv run python scripts\validate_project.py quick
```

Direct pytest also works, for example `uv run python -m pytest tests\api\test_runtime.py -q`, but `uv run python run_pytest.py ...` is preferred for source-checkout validation because it sets the source import path and manages pytest temporary directories more predictably. Set `EPCSAFT_PYTEST_TEMP_ROOT` when you want the wrapper to use an opt-in external pytest temp root instead of its default repo-local generated temp area.

The default source-checkout validation sequence is sync, normal native build, then `uv run python scripts\validate_project.py quick`. Use `uv run python scripts\validate_project.py confidence` before handoff when runtime confidence matters.

For maintainers, the [development workflow guide](docs/pages/development_workflows.rst) is the source-of-truth command matrix for setup, fast rebuilds, focused tests, profiling, packaging, and repair-only cleanup. The [project structure guide](docs/pages/project_structure.rst) explains where package code, reusable reference data, and analysis workflows belong.
For downstream projects that install this checkout as a uv path dependency, use the [downstream local install guide](docs/pages/downstream_local_installs.rst). For suspected package bugs, use the [downstream dependency protocol](docs/pages/downstream_dependency_protocol.rst) to file a complete GitHub issue and validation command.

GitHub Actions are intentionally lightweight while this package is still in local development. Pull requests run a fast workflow/package contract smoke; pushes to `main` run a Windows downstream path-install smoke on Python 3.13. Full Linux/macOS/Windows wheel and sdist packaging checks are available as a manual workflow dispatch when release confidence is needed.

For the standard validation loops:

```powershell
uv run python scripts\validate_project.py quick
uv run python scripts\validate_project.py confidence
uv run python scripts\validate_project.py docs
uv run python run_pytest.py --list-slices
uv run python run_pytest.py -q
uv run python run_pytest.py --runtime -q
uv run python run_pytest.py --generic -q
uv run python run_pytest.py --confidence -q
uv run python run_pytest.py --equilibrium-confidence -q -s
uv run python run_pytest.py --all -q
uv run python run_pytest.py --profile -q
uv run python run_pytest.py --profile-full -q -s
```

`validate_project.py` is the high-level validation orchestrator; `run_pytest.py` remains the lower-level test selector. `run_pytest.py -q` is the default fast contract suite: it runs representative API, native, regression, equilibrium, and workflow smoke checks without regenerating plots or full scientific reproduction suites. `--runtime` runs runtime API plus native contract tests. `--generic` runs the same fast contract target list as the default route. `--confidence` adds a few native runtime contracts and is the default runtime-confidence check before handoff. `--equilibrium-confidence` is a bounded Khudaida fixture plus cached fixed-phase residual contract; native confidence solving and full electrolyte reports remain explicit opt-ins. Generated plot output tests are manual, targeted pytest runs rather than named validation slices. `--all` is the explicit exhaustive historical suite and can take many minutes. `--profile` enables and runs the quick opt-in runtime-only profiling check. `--profile-full` runs the slower opt-in runtime, MIAC, and regression profile suite; it can take about a minute locally, so use a runner timeout of at least 120 seconds. To keep pytest temp files outside the repo for an opt-in run, set `EPCSAFT_PYTEST_TEMP_ROOT`, for example:

```powershell
$env:EPCSAFT_PYTEST_TEMP_ROOT = Join-Path $env:TEMP 'epcsaft-pytest'
uv run python run_pytest.py --confidence -q
```

For repeated calculations, create an `ePCSAFTMixture` and `ePCSAFTState` once and reuse them inside loops. The profile report flags full rebuilds when they dominate runtime.

Use `uv run python scripts\build_epcsaft.py --clean` only as a repair step for stale CMake state or stale/locked `_core` artifacts. If a `_core*.pyd` is locked, stop the importing Python/test/IDE process before running the clean repair.

`CMakePresets.json` is optional Windows MinGW convenience for IDEs and manual CMake use. The canonical local native build remains `uv run python scripts\build_epcsaft.py`.

Build distributable artifacts at the packaging boundary:

```powershell
uv run python scripts\build_dist.py
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
- `data/reference/`: source-checkout reference data, including example parameter datasets under `data/reference/epcsaft_parameters/`
- `analyses/`: source-owned paper validation, fit validation, and figure workflows with local `scripts/`, `data/`, and `results/`
- `docs/`: user documentation and reference material
- `scripts/`: uv/CMake build, doctor, and validation helpers

## Public API

The main entry points are `ePCSAFTMixture`, `ePCSAFTState`, `create_parameter_template`, `fit_pure_neutral(...)`, `fit_pure_ion(...)`, `fit_binary_pair(...)`, `fit_mea_co2_h2o_electrolyte(...)`, `solve_reactive_speciation(...)`, `runtime_build_info()`, `capabilities()`, and the structured result objects returned by solver-style methods.

`create_parameter_template(...)` creates a blank dataset folder with the expected `pure/`, `mixed/`, and `user_options.json` layout. After you fill in the files, `ePCSAFTMixture.from_dataset(...)` can load the folder path you created yourself.

The regression workflow is record-driven: `fit_pure_neutral(...)` fits nonassociating neutral-component `m`, `s`, and `e`; `fit_pure_ion(...)` fits ion `s`/`e` and optional `d_born`; `fit_binary_pair(...)` fits V1 VLE `k_ij` values; and `fit_mea_co2_h2o_electrolyte(...)` exposes the opt-in MEA-CO2-H2O electrolyte pure-parameter benchmark. Use `write_fit_result(...)` when you want to persist a fit back into a user-owned dataset folder.

Homogeneous reaction chemistry is available through `solve_reactive_speciation(...)`, where the caller supplies species, balances, reactions, equilibrium constants, and a mixture factory for native activity-coefficient evaluations. Electrolyte bubble-pressure and composed reactive electrolyte bubble-pressure entry points are native-backend placeholders and raise `InputError` until the corresponding C++ solvers exist; the package does not expose Python equilibrium solver fallbacks.

## Documentation

- [Start here](docs/pages/README.rst)
- [Getting started](docs/pages/getting_started.rst)
- [Downstream local installs](docs/pages/downstream_local_installs.rst)
- [Development workflow guide](docs/pages/development_workflows.rst)
- [Create your own parameter folder](docs/pages/user_parameter_templates.rst)
- [Parameter regression guide](docs/pages/parameter_regression.rst)
- [User options reference](docs/pages/user_options.rst)
- [Task-based package guide](docs/pages/package_guide.rst)
- [Downstream dependency protocol](docs/pages/downstream_dependency_protocol.rst)
- [Electrolyte VLE and reactive speciation](docs/pages/electrolyte_vle_reactive_workflow.rst)
- [API reference](docs/pages/api_reference.rst)

## License

GNU General Public License v3.0
