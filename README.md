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
- `docs/latex/`: separate Git submodule for the Overleaf-linked LaTeX project

## Runtime Parameter Data

The runtime dataset loader lives in `pcsaft.parameters`, and the package bundles the dataset files needed at runtime.

The preferred object API is `PCSAFTMixture` plus `PCSAFTState`. Use `PCSAFTMixture.from_dataset(...)` for dataset-backed systems and `PCSAFTMixture.from_params(...)` when you already have a resolved parameter payload.

## Development Workflow

For local work inside this repo, use an editable install and rerun it when package metadata or compiled sources change:

```bash
python scripts/build_pcsaft.py
python run_pytest.py
```

`python scripts/build_pcsaft.py` now skips the reinstall when the editable build is already current. Use `python scripts/build_pcsaft.py --force` to force a fresh editable reinstall, or `python run_pytest.py --force-build` to force that rebuild before tests.

Analysis and validation scripts are expected to run from the repo-named `PC-SAFT` Conda environment with `pcsaft` installed editable. A source checkout by itself is not a supported package-import path.

The paper-validation and exploratory analysis scripts intentionally remain outside `pytest`; they are workspace tools, not package-unit tests.

## LaTeX / Overleaf Workflow

The LaTeX project lives in the `docs/latex` submodule and keeps its own Overleaf remote and branch history.

Typical workflow:

```bash
git -C docs/latex status
git -C docs/latex add -A
git -C docs/latex commit -m "Update paper draft"
git -C docs/latex push origin master
git add docs/latex
git commit -m "Bump latex submodule pointer"
```

Clone and initialize a fresh checkout with:

```bash
git clone --recurse-submodules <repo-url>
git submodule update --init --recursive
```

If the submodule is already cloned but needs to be refreshed after a parent checkout:

```bash
git submodule update --init --recursive
```

## Codex / Agent Workflow

For Codex agents and other automated local tooling on Windows, this repo overrides the generic "default to `conda run`" policy: prefer the repo-local environment wrapper first.

```bash
pwsh -NoProfile -ExecutionPolicy Bypass -File scripts/run_in_repo_env.ps1 scripts/codex_doctor.py
pwsh -NoProfile -ExecutionPolicy Bypass -File scripts/run_in_repo_env.ps1 run_pytest.py
```

The wrapper resolves the repo-named Conda environment's `python.exe` directly and only falls back to `conda run` if needed, which avoids some Windows temp-file issues seen with `conda run`.

Use raw `conda run -n PC-SAFT ...` only when you already know it is healthy in the current session. In that case, these are equivalent:

```bash
conda run -n PC-SAFT python scripts/codex_doctor.py
conda run -n PC-SAFT python run_pytest.py
```

`scripts/codex_doctor.py` reports which checkout `pcsaft` is importing from, whether that checkout is compatible with the current repo, and whether a rebuild is needed before testing.

`scripts/build_pcsaft.py` now treats a same-commit editable install from another worktree as compatible for testing, so agents do not need to reinstall the package just because the active environment still points at a sibling checkout with the same git `HEAD`.

When an editable reinstall is required, the build helper uses repo-local pip temp/cache directories under `build/` instead of machine-global temp paths.

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
