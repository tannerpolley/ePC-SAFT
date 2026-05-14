# Scripts Layout

The root `scripts/` package is intentionally split between shared helper modules and runnable entrypoints.

Canonical rules:

- Keep shared helper modules at the root when analysis code imports them directly.
- Put runnable repo workflows under purpose-specific subfolders such as `dev/`, `benchmarks/`, `validation/`, `docs/`, `data/`, and `support/`.
- Move study-specific or figure-specific tooling into `analyses/<study>/scripts/` instead of keeping it under the repo-root `scripts/` tree.

Current layout:

- `scripts/_env.py`, `scripts/_epcsaft_oop.py`, `scripts/plot_outputs.py`: shared analysis helpers.
- `scripts/sandbox_tempfile_site/`: shared support for sandbox-safe temporary directories during build/package workflows.
- `scripts/dev/`: environment setup, native build, repair, validation, and packaging entrypoints.
- `scripts/benchmarks/`: benchmark and profiling entrypoints.
- `scripts/validation/`: opt-in specialty validation entrypoints.
- `scripts/docs/`: documentation, equation-registry, and LaTeX-mirror maintenance entrypoints.
- `scripts/data/`: reference-data curation helpers.
- `scripts/support/`: issue triage and other support tooling.
