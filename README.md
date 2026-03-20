# PC-SAFT

`pcsaft` is a Cython/C++ implementation of the PC-SAFT equation of state with dipole, association, and electrolyte terms.

This repository now serves three roles at once:

- the installable `pcsaft` Python package
- the runtime parameter-data source for ePC-SAFT datasets
- the in-repo analysis and paper-validation workspace used to recreate figures and compare model behavior

## Install

Standard installs:

```bash
pip install .
```

Editable installs for local development and analysis:

```bash
pip install -e . --no-build-isolation
```

The build requires a C++ toolchain, Cython, NumPy, and the vendored Eigen headers included in this repository.

## Package Layout

- `src/pcsaft/`: installable runtime package, Cython/C++ sources, and packaged parameter datasets
- `data/`: workspace and analysis datasets that are not required for package runtime
- `scripts/`: validation, fitting, and paper-reproduction workflows
- `tests/`: core package regression tests
- `docs/`: user documentation and paper/reference material

## Runtime Parameter Data

The runtime dataset loader lives in `pcsaft.parameters`, and the package bundles the dataset files needed at runtime.

## Development Workflow

For local work inside this repo, use an editable install and rerun it when package metadata or compiled sources change:

```bash
python scripts/build_pcsaft.py
python run_pytest.py
```

`python scripts/build_pcsaft.py` now skips the reinstall when the editable build is already current. Use `python scripts/build_pcsaft.py --force` to force a fresh editable reinstall, or `python run_pytest.py --force-build` to force that rebuild before tests.

Analysis and validation scripts are expected to run from the active `PC-SAFT` environment with `pcsaft` installed editable. A source checkout by itself is not a supported package-import path.

The paper-validation and exploratory analysis scripts intentionally remain outside `pytest`; they are workspace tools, not package-unit tests.

## MIAC Data Workflow

Canonical MIAC experimental datasets live under `data/MIAC`.

- Canonical MIAC CSVs: `data/MIAC/<solvent_system>/*.csv`
- Canonical MIAC fit outputs: `data/MIAC/<solvent_system>/miac_fits/*.png`
- Canonical MIAC_m fit outputs: `data/MIAC/<solvent_system>/miac_m_fits/*.png`

Regenerate the MIAC figures with:

```bash
conda run -n PC-SAFT python scripts\\fits\\validate_miac_fits.py
```

The legacy `data/MIAC_m`, `data/MIAC_combined`, and `data/MIAC/*/plot_fits` layouts are no longer part of the supported workflow.
