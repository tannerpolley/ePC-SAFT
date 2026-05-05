# Paper Validation Scripts

This folder contains paper-reproduction and diagnostic analysis assets. It is separate from the core package build and test workflow.

## Gallery Workflow

Regenerate the local PNG browser with:

```powershell
uv run python scripts/paper_validation/tools/build_analysis_galleries.py
```

The generator writes the master gallery to `docs/plots/index.html` from the tracked `docs/plots/manifest.json`. Open `docs/plots/index.html` locally to browse generated PNG outputs when they exist.

The gallery is a single root page. It lists source-owned plot outputs under `scripts/**/out`, `tests/**/out`, and `src/**/out` and links the matching SVG and CSV assets when they exist. These generated assets are ignored; refresh the manifest with `uv run python scripts/build_plot_manifest.py --refresh` only after intentionally regenerating local outputs.

## CSV-Backed Figure Data

`2025_Figiel_analysis` stores generated model/literature figure payloads under each figure folder's local `out/` directory.

Regenerate the CSV baseline with:

```powershell
uv run python scripts/paper_validation/2025_Figiel_analysis/generate_figure_data.py
```

Validate the committed CSV baseline, regenerate the PNGs from CSV, and rebuild the gallery with:

```powershell
uv run python scripts/paper_validation/2025_Figiel_analysis/validate_figure_data.py
```

The 2025 Figiel plot scripts should read generated CSV payloads and should not call ePC-SAFT model evaluation directly during rendering. Keep source/literature CSV inputs under the analysis folder; generated numeric figure payloads belong under the owning figure folder's ignored `out/` directory.

All generated PNGs should also have a companion `<png-stem>_plot_data.csv` in the same `out/` folder. The shared `scripts.plot_outputs.save_plot_figure(...)` helper writes this CSV when a script saves a Matplotlib figure.

## Folder Roles

- `tools/`: shared tooling for paper-validation maintenance, starting with the gallery generator.
- `*_analysis/`: per-paper figure scripts and digitized/source data. Generated PNGs and generated figure CSVs live under each owning folder's ignored `out/`.
- `2020_Bulow_analysis/**/diagnostics`: exploratory and audit material. Review and tag these before pruning; do not delete them as routine cleanup.

Keep figure scripts and diagnostics in place unless a cleanup explicitly archives or removes a complete, reviewed analysis slice.
