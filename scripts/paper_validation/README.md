# Paper Validation Scripts

This folder contains paper-reproduction and diagnostic analysis assets. It is separate from the core package build and test workflow.

## Gallery Workflow

Regenerate the local PNG browser with:

```powershell
uv run python scripts/paper_validation/tools/build_analysis_galleries.py
```

The generator writes the master gallery to `docs/plots/index.html` and nested galleries under `docs/plots/**/index.html`. Open `docs/plots/index.html` locally to browse generated PNG outputs.

The master page and parent folders show only navigation. A page renders PNGs only for files directly inside the current folder; use the subfolder dropdown or folder links to drill down.

## CSV-Backed Figure Data

`2025_Figiel_analysis` is the first validation slice where model/literature figure payloads are stored as tracked CSV artifacts under `docs/plots/paper_validation/2025_Figiel/<figure>/data/`.

Regenerate the CSV baseline with:

```powershell
uv run python scripts/paper_validation/2025_Figiel_analysis/generate_figure_data.py
```

Validate the committed CSV baseline, regenerate the PNGs from CSV, and rebuild the gallery with:

```powershell
uv run python scripts/paper_validation/2025_Figiel_analysis/validate_figure_data.py
```

The 2025 Figiel plot scripts should read these generated CSV payloads and should not call ePC-SAFT model evaluation directly. Keep source/literature CSV inputs under the analysis folder; generated numeric figure payloads belong under `docs/plots`.

All generated PNGs under `docs/plots/**` should also have a companion `data/<png-stem>_plot_data.csv`. The shared `scripts.plot_outputs.save_plot_figure(...)` helper writes this CSV when a script saves a Matplotlib figure. Legacy or unavailable outputs may use an `existing_png_backfill` row until their source workflow is made reproducible again.

## Folder Roles

- `tools/`: shared tooling for paper-validation maintenance, starting with the gallery generator.
- `*_analysis/`: per-paper figure scripts and digitized/source data. Generated PNGs now live under `docs/plots/paper_validation/<paper>/...`.
- `2020_Bulow_analysis/**/diagnostics`: exploratory and audit material. Review and tag these before pruning; do not delete them as routine cleanup.

Keep figure scripts and diagnostics in place unless a cleanup explicitly archives or removes a complete, reviewed analysis slice.
