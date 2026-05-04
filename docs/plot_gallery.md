# Generated Plot Gallery

`docs/plots/index.html` is the single tracked gallery entrypoint. Generated PNG,
SVG, CSV, and related plot assets live in source-local `out/` folders beside the
scripts, tests, or package modules that create them.

Use:

```powershell
uv run python scripts/paper_validation/tools/build_analysis_galleries.py
```

Open `docs/plots/index.html` to browse all generated plot groups. The gallery is
a single-page file explorer: expand folders in the left sidebar, check one or
more folders, and the page displays every PNG in those selected subtrees without
navigating to separate folder pages.

To rebuild and serve the gallery as a localhost site:

```powershell
uv run python scripts/paper_validation/tools/serve_plot_gallery.py
```

The server serves the repository root so `docs/plots/index.html` can load assets
from `scripts/**/out`, `tests/**/out`, and `src/**/out`.

Folder roles:

- `scripts/paper_validation/**/out`: paper-reproduction and diagnostic figures.
- `scripts/fits/out`: MIAC, dielectric, osmotic, and presentation fit figures.
- `tests/plots/out`: opt-in regression/test visualizations.
- `src/**/out`: package-owned report figures, such as equilibrium confidence.

Plot notation conventions:

- Use proper thermodynamic/math notation whenever a standard symbol exists.
- Keep units outside the math expression where practical, for example
  `r"$\rho$ / mol m$^{-3}$"` or `r"$P$ / Pa"`.
- Apply this to plot titles, axes, legends, hover labels, and colorbars.

To rebuild the gallery bundle for one source test with an explicit plot recipe:

```powershell
uv run python scripts/build_test_plot_gallery.py tests/equilibrium/test_electrolyte_lle.py
```

The command resolves the source test through `tests/plots/plot_registry.py`,
runs the registered plot-producing pytest target, validates that every affected
PNG under `tests/plots/out/**` has a sibling `<stem>_plot_data.csv`, creates
missing SVG companions, refreshes `docs/plots/index.html`, and writes the asset
report under `tests/plots/out/plot_asset_report.csv`.
