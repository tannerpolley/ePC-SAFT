Project Structure
=================

This repository is organized as a Python package with separate source-owned analysis workflows.

Package Surfaces
----------------

``src/epcsaft/``
    Public Python package code, pure-Python helpers, and native-extension wrappers.

``tests/``
    Package/API/native/workflow contracts. Default tests should stay fast and should not reproduce full scientific studies, regenerate plot galleries, or run long fitting/equilibrium sweeps.

``scripts/``
    Repository tooling only: native builds, doctor checks, validation orchestration, packaging, docs, reference-data curation, LaTeX sync, and issue triage. Analysis-specific scripts belong under ``analyses/<short_id>/scripts/``.

``data/reference/``
    Reusable checkout data: parameter datasets, benchmark records, literature tables, and curation sources that may be shared by several analyses or package tests.

``docs/``
    Package documentation and archival paper material. Archival PDFs and paper notes remain under ``docs/papers/``.

Analysis Workflows
------------------

Scientific reproductions, validation plots, fits, and paper figure workflows live under ``analyses/<short_id>/``.

Use this layout for new analyses:

.. code-block:: text

   analyses/<short_id>/
     README.md
     analysis.yaml
     references.bib
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

Only create optional folders when they are useful. Keep each analysis self-contained: local scripts should read local analysis inputs or ``data/reference/`` and write outputs under the same analysis folder.

Data Ownership
--------------

Use ``data/reference/`` for stable reusable inputs:

- package example parameter datasets under ``data/reference/epcsaft_parameters/``
- equilibrium benchmark fixtures under ``data/reference/equilibrium_benchmarks/``
- reusable literature data such as ``MIAC``, ``osmotic``, ``pure_component``, and regression tables

Use ``analyses/<short_id>/data/input/`` for hand-curated or digitized inputs that are specific to one analysis. Use ``analyses/<short_id>/data/processed/`` for canonical plot-ready intermediate tables when those tables should be retained as analysis inputs.

Output Policy
-------------

Generated outputs use ``results/``:

- ``results/runs/`` for ignored run-specific payloads, logs, sweeps, and exploratory output
- ``results/<plot_set>/`` for curated plot sets that are intentionally retained

Each curated plot set keeps the figure, exact plotted data snapshot, and editable Matplotlib sidecar together:

.. code-block:: text

   results/response_curve/
     response_curve.csv
     response_curve.svg
     response_curve.png
     response_curve.mpl.yaml

Use one meaningful plot-set folder per figure or figure family. This layout is designed so external tools such as ``mplgallery`` can discover analysis-owned plots directly without a package-owned gallery index, server, or manifest.

Do not add new gallery index, server, or manifest workflows to this package. External visualization tools should discover analysis-owned assets directly.
