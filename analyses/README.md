# Analysis Workflows

This directory contains source-controlled scientific analysis, validation, and figure workflows that are useful for developing and checking `epcsaft` but are not package runtime code.

The analysis taxonomy is:

```text
analyses/
  _template/
  paper_validation/
    native/
    application/
  data_validation/
  package_validation/
```

Each real analysis should be self-contained inside one of those category folders:

```text
analyses/<category>/<short_id>/
  README.md
  analysis.yaml
  config/
  scripts/
  figures/
    <figure_id>/
      input/
      output/
      scripts/
  notebooks/
  tests/
```

Only create optional folders when the analysis needs them. Stable literature inputs that are reused by multiple analyses belong under `data/reference/`; analysis-specific inputs stay local under the owning `figures/<figure_id>/input/` folder.

Generated figure-local run payloads belong under `figures/<figure_id>/output/runs/` and are ignored. Retained model CSVs, plotted data snapshots, rendered figures, and `.mpl.yaml` sidecars belong together under the owning `figures/<figure_id>/output/` folder. The package-level `scripts/` directory is reserved for repo tooling such as builds, validation, data curation, packaging, and docs.
