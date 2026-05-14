# MIAC Fits

Mean ionic activity coefficient fitting and validation workflows.

Reusable input datasets are loaded from `data/reference/MIAC/**`. For new or reorganized figure workflows, keep each retained figure under its own `figures/<figure_id>/` folder with local `input/`, `output/`, and `scripts/` subfolders. Historical MIAC fit outputs now belong under the figure-owned layout instead of the old analysis-level plot-set roots.

The main validation entrypoint is:

```powershell
uv run python analyses\miac_fits\scripts\validate_miac_fits.py
```
