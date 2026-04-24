# Paper Validation Scripts

This folder contains paper-reproduction and diagnostic analysis assets. It is separate from the core package build and test workflow.

## Gallery Workflow

Regenerate the local PNG browser with:

```powershell
uv run python scripts/paper_validation/tools/build_analysis_galleries.py
```

The generator writes `index.html` in this folder and in each `*_analysis` folder. Open `scripts/paper_validation/index.html` locally to browse generated PNG outputs.

## Folder Roles

- `tools/`: shared tooling for paper-validation maintenance, starting with the gallery generator.
- `*_analysis/`: per-paper figure scripts, digitized data, generated PNGs, and local `index.html` gallery pages.
- `2020_Bulow_analysis/**/diagnostics`: exploratory and audit material. Review and tag these before pruning; do not delete them as routine cleanup.

Keep figure scripts and diagnostics in place unless a cleanup explicitly archives or removes a complete, reviewed analysis slice.
