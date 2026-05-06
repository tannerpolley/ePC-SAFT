# MIAC Fits

Mean ionic activity coefficient fitting and validation workflows.

Reusable input datasets are loaded from `data/reference/MIAC/**`. Generated figures belong under `analyses/miac_fits/results/final/figures/`; durable generated tables or reports belong under `results/final/tables/` and `results/final/reports/`. Run-specific or exploratory outputs should be written under `results/runs/`.

The main validation entrypoint is:

```powershell
uv run python analyses\miac_fits\scripts\validate_miac_fits.py
```
