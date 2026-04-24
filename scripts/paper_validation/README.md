# Paper Validation Scripts

This folder contains paper-reproduction and diagnostic analysis assets. It is separate from the core package build and test workflow.

## Gallery Workflow

Regenerate the local PNG browser with:

```powershell
uv run python scripts/paper_validation/tools/build_analysis_galleries.py
```

The generator writes the master gallery to `docs/plots/index.html` and nested galleries under `docs/plots/**/index.html`. Open `docs/plots/index.html` locally to browse generated PNG outputs.

## Folder Roles

- `tools/`: shared tooling for paper-validation maintenance, starting with the gallery generator.
- `*_analysis/`: per-paper figure scripts and digitized/source data. Generated PNGs now live under `docs/plots/paper_validation/<paper>/...`.
- `2020_Bulow_analysis/**/diagnostics`: exploratory and audit material. Review and tag these before pruning; do not delete them as routine cleanup.

Keep figure scripts and diagnostics in place unless a cleanup explicitly archives or removes a complete, reviewed analysis slice.
