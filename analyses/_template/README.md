# Analysis Template

Use this folder shape for new modeling, validation, and figure workflows. Copy the structure, rename the analysis folder with a short lowercase ID, and delete optional folders that are not needed.

Put each retained figure under `figures/<figure_id>/`. Hand-curated or digitized inputs go under `figures/<figure_id>/input/`. Generated model tables, plotted CSV snapshots, rendered figures, and `.mpl.yaml` sidecars go under `figures/<figure_id>/output/`. Disposable figure-local runs go under `figures/<figure_id>/output/runs/`.
