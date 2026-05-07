# Analysis Workflows

This directory contains source-controlled scientific analysis, validation, and figure workflows that are useful for developing and checking `epcsaft` but are not package runtime code.

Each analysis should be self-contained:

```text
analyses/<short_id>/
  README.md
  analysis.yaml
  config/
  scripts/
  data/
    input/
    raw/
    processed/
  results/
    runs/
    <plot_set>/
      <plot_set>.csv
      <plot_set>.svg
      <plot_set>.png
      <plot_set>.mpl.yaml
  notebooks/
  tests/
```

Only create optional folders when the analysis needs them. Stable literature inputs that are reused by multiple analyses belong under `data/reference/`; analysis-specific inputs stay local under `analyses/<short_id>/data/input/` or beside the figure script when the current workflow has not yet been normalized.

Generated run payloads belong under `results/runs/` and are ignored. Curated plot-producing outputs belong under `results/<plot_set>/`, with the figure, exact plotted data snapshot, and `<plot_set>.mpl.yaml` sidecar in the same folder so `mplgallery` can discover plot sets directly. The package-level `scripts/` directory is reserved for repo tooling such as builds, validation, data curation, packaging, and docs.
