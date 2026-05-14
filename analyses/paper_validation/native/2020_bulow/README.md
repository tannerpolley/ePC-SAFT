# 2020 Bulow Analysis

Paper-validation workflow for the Bulow 2020 figures.

Canonical Gibbs-energy and MIAC source tables are staged into each figure's `input/`
folder by the figure-local `generate_data.py` entrypoints. Plot and diagnostic scripts
should read figure-local inputs or figure-owned outputs only.

Entry point:

```powershell
uv run python analyses\paper_validation\native\2020_bulow\scripts\run_all.py
```
