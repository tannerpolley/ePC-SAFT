# 2012 Paper Data Gaps

The current folder currently reproduces these 2012 validation figures with machine-readable data:

- Figure 2d (NaI in methanol/ethanol): `figure_2d/plot_figure_2d.py`
- Figure 3 (NaI solvent activity coefficients): `figure_3/plot_figure_3.py`
- Figure 5 (NaCl in water/ethanol): `figure_5/plot_figure_5.py` using `figure_5/esteso_1989_table1_selected_series.csv` from Esteso et al. (1989) Table I.
- Figure 6 (NaBr in water/methanol)
- Figure 7 (NaCl in water/methanol/ethanol)

Resolved by local figure-input mirroring in this folder:

- Figure 2d now stages `ethanol-NaI.csv` and `methanol-NaI.csv` into `figures/figure_2d/input/`.
- Figure 6 now stages `water-methanol-NaBr.csv` into `figures/figure_6/input/`.
- Figure 7 already uses figure-local copies of the Hernandez 2007 and dielectric-mixture CSVs under `figures/figure_7/input/`.

Still missing machine-readable data or dedicated recreation work for full-paper reproduction:

- Figure 1: NH4Br/ethanol osmotic + MIAC pair (osmotic curve not available as CSV).
- Figure 2a/2b/2c: density, vapor pressure, and osmotic-coefficient panels.
- Figure 4: water/ethanol/NaI density panel.

Deferred repo-local digitization source:

- `docs/papers/pdf/Held, Prinz, Wallmeyer, Sadowski - 2012 - Measuring and modeling alcohol-salt systems.pdf`

