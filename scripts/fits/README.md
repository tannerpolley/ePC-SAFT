# Fit And Validation Scripts

This folder contains research and validation scripts, not the core package build workflow.

Generated CSV/PNG/SVG outputs belong under `scripts/fits/out/**`, not under `data/**`. These outputs are ignored. The external plot browser should discover source-local `out/` folders directly.

## Active Validation

- `validate_miac_fits.py`: main dataset-driven MIAC validation and plot generator.
- `profile_miac_runtime.py`: opt-in runtime profiler used by `tests/profile/test_miac_profile.py`.

## Retained Analysis Helpers

- `plot_2025_presentation_miac.py`: presentation-oriented MIAC plots built from the validation workflow.
- `plot_2025_salt_solvent_comparison.py`: comparison plots for 2025 salt/solvent analysis.

## Historical Plot Generators

- `dielc_fit.py` and `dielc_diff.py`: dielectric fit and derivative visualizations retained for provenance.
- `validation_2014_repro.py`: Held 2014-style osmotic reproduction plot retained for provenance.

Keep this folder stable unless the corresponding datasets, presentation references, or runtime-profile tests are updated.
