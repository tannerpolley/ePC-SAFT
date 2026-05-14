Project Structure
=================

This repository is organized as a Python package with separate source-owned analysis workflows.

Package Surfaces
----------------

``src/epcsaft/``
    Public Python package code, pure-Python helpers, and native-extension wrappers.

``tests/``
    Package/API/native/workflow contracts. Top-level domains stay broad, with second-level subgroup folders such as ``tests/api/runtime/``, ``tests/equilibrium/electrolyte/``, ``tests/native/cppad/``, and ``tests/workflows/repo/`` when a domain gets large. Default tests should stay fast and should not reproduce full scientific studies, regenerate plot galleries, or run long fitting/equilibrium sweeps.

``scripts/``
    Repository tooling only: native builds, doctor checks, validation orchestration, packaging, docs, reference-data curation, LaTeX sync, and issue triage. Analysis-specific coordinators may live under ``analyses/<category>/<short_id>/scripts/``, while figure-local generation and rendering scripts belong under ``analyses/<category>/<short_id>/figures/<figure_id>/scripts/``.

``data/reference/``
    Reusable checkout data: parameter datasets, benchmark records, literature tables, and curation sources that may be shared by several analyses or package tests.

``docs/``
    Package documentation and archival paper material. Archival PDFs and paper notes remain under ``docs/papers/``.

Analysis Workflows
------------------

Scientific reproductions, validation plots, fits, and paper figure workflows live under categorized analysis roots such as ``analyses/paper_validation/native/<short_id>/`` or ``analyses/data_validation/<short_id>/``.

Use this layout for new analyses:

.. code-block:: text

   analyses/<category>/<short_id>/
     README.md
     analysis.yaml
     references.bib
     config/
     scripts/
     figures/
       <figure_id>/
         input/
         output/
         scripts/
     notebooks/
     tests/

Only create optional folders when they are useful. Keep each analysis self-contained: local scripts should read local analysis inputs or ``data/reference/`` and write outputs under the same analysis folder.

Data Ownership
--------------

Use ``data/reference/`` for stable reusable inputs:

- package example parameter datasets under ``data/reference/epcsaft_parameters/``
- equilibrium benchmark fixtures under ``data/reference/equilibrium_benchmarks/``
- reusable literature data such as ``MIAC``, ``osmotic``, ``pure_component``, and regression tables

Use ``data/reference/`` for stable reusable inputs shared by multiple analyses. Use ``analyses/<category>/<short_id>/figures/<figure_id>/input/`` for hand-curated, digitized, or figure-owned parameter snapshots. Use ``analyses/<category>/<short_id>/figures/<figure_id>/output/`` for generated model tables and the exact plotted data retained with that figure.

Output Policy
-------------

Generated figure outputs use figure-owned ``output/`` folders:

- ``figures/<figure_id>/output/runs/`` for ignored run-specific payloads, logs, sweeps, and exploratory output
- ``figures/<figure_id>/output/`` for curated generated tables, exact plotted data snapshots, rendered figures, and editable sidecars that are intentionally retained

Each figure-owned output folder keeps the figure, exact plotted data snapshot, and editable Matplotlib sidecar together:

.. code-block:: text

   figures/response_curve/output/
     response_curve.csv
     response_curve.svg
     response_curve.png
     response_curve.mpl.yaml

Use one meaningful figure folder per figure or figure family. This layout is designed so external tools such as ``mplgallery`` can discover analysis-owned plots directly without a package-owned gallery index, server, or manifest.

Do not add new gallery index, server, or manifest workflows to this package. External visualization tools should discover analysis-owned assets directly.
