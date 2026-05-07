# MIAC Fits

Mean ionic activity coefficient fitting and validation workflows.

Reusable input datasets are loaded from `data/reference/MIAC/**`. Curated generated plot sets belong under `analyses/miac_fits/results/<plot_set>/`, with figure files, exact plotted data snapshots, and `<plot_set>.mpl.yaml` sidecars kept together for `mplgallery` discovery. Run-specific or exploratory outputs should be written under `results/runs/`.

The main validation entrypoint is:

```powershell
uv run python analyses\miac_fits\scripts\validate_miac_fits.py
```
