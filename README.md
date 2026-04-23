# ePC-SAFT

`epcsaft` is a Cython/C++ implementation of the ePC-SAFT equation of state with dipole, association, and electrolyte terms.

## Install

Most users should install from PyPI:

```bash
pip install epcsaft
```

`epcsaft` includes a compiled Cython/C++ extension. If a wheel is available for your Python version and platform, pip installs it automatically. If a wheel is not available, pip falls back to a source build, which requires a working native build toolchain.

For a source checkout of this repository:

```bash
pip install .
```

For editable development from this source tree:

```bash
python scripts/build_epcsaft.py
```

That command expects an editable checkout to already be importable. It rebuilds the in-place Cython/C++ extension when native sources are stale, and otherwise exits without touching the environment. For a first install or a broken editable link, use:

```bash
python scripts/build_epcsaft.py --reinstall-editable
```

If you want to call pip directly, use:

```bash
pip install -e . --no-build-isolation --config-settings editable_mode=compat
```

If you are iterating on the native code, rerun `python scripts/build_epcsaft.py` to refresh the editable install. On Windows, the compiled extension may appear as a `.pyd` during local builds, but it is a build artifact, not a file users normally copy into a project by hand.

## Example

```python
import numpy as np
from epcsaft import ePCSAFTMixture

mixture = ePCSAFTMixture.from_params(
    {"m": np.asarray([2.8149]), "s": np.asarray([3.7169]), "e": np.asarray([285.69])},
    species=["Toluene"],
)
state = mixture.state(T=320.0, x=np.asarray([1.0]), P=101325.0)
# Pressure-based states solve and cache density during construction.
print(state.density())               # mol/m^3 by default
print(state.density(units="mass"))   # kg/m^3 when MW is available
print(state.molar_density())         # explicit molar-density alias
print(state.ares())                  # short alias for residual_helmholtz()
print(state.ares(return_contribution_terms=True)["terms"]["hc"])
```

## Package Layout

- `src/epcsaft/`: installable runtime package and Cython/C++ sources
- `data/epcsaft_parameters/`: source-checkout example parameter datasets for inspection, comparison, and tests
- `data/`: other datasets and figures that are not required by the package
- `docs/`: user documentation and reference material

## Public API

The main entry points are `ePCSAFTMixture`, `ePCSAFTState`, `create_parameter_template`, `fit_pure_neutral(...)`, and the structured result objects returned by solver-style methods.

`create_parameter_template(...)` creates a blank dataset folder with the expected `pure/`, `mixed/`, and `user_options.json` layout. After you fill in the files, `ePCSAFTMixture.from_dataset(...)` can load the folder path you created yourself. If you are working from a checkout of this repository, the example folders under `data/epcsaft_parameters/` are available for inspection and comparison.

Phase 1 of the regression workflow is intentionally narrow: `fit_pure_neutral(...)` fits only nonassociating neutral-component `m`, `s`, and `e` against liquid-density and vapor-pressure data. The public Python surface stays the same, but the optimizer now runs natively with the package's least-squares workflow only. Use `write_fit_result(...)` when you want to persist a fit back into a user-owned dataset folder. Ion and binary regression are deferred for now.

## Documentation

- [Start here](docs/README.rst)
- [Getting started](docs/getting_started.rst)
- [Create your own parameter folder](docs/user_parameter_templates.rst)
- [Parameter regression guide](docs/parameter_regression.rst)
- [User options reference](docs/user_options.rst)
- [Task-based package guide](docs/package_guide.rst)
- [API reference](docs/api_reference.rst)

## Author

Tanner Polley

## License

GNU General Public License v3.0


