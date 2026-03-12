# 2014 Paper Data Gaps

Current reproducible outputs in this folder:

- Figure 1: NaCl/KBr osmotic coefficients at 298.15 K (`figure_1/plot_figure_1.py`).
- Figure 2: model curves for NaCl/KCl at 273.15 K and 298.15 K with digitized panel data for KCl (273.15 K, 298.15 K) and NaCl (273.15 K), plus exact NaCl 298.15 K data (`figure_2/plot_figure_2.py`, `data/osmotic/water/held_2014_figure_2_digitized.csv` mirrored to `figure_2/figure_2_digitized.csv`).
- Figure 3: LiAc/NaAc/KAc osmotic coefficients with strategy-2 model and a classical ePC-SAFT reference using the new `held_2009` dataset (`figure_3/plot_figure_3.py`) and the richer acetate data file `data/osmotic/water/LiAc-NaAc-KAc.csv` mirrored to `figure_3/LiAc-NaAc-KAc.csv`.
- Figure 4a: amino-acid/KCl osmotic coefficients (alanine, glycine) with strategy-2 model and strategy-1-like reference model (`figure_4/plot_figure_4a.py`) using `figure_4a_digitized.csv` data points.
- Figure 4b (digitized): amino-acid solubility panel reproduced from digitized paper data/model traces (`figure_4/plot_figure_4b.py`, `figure_4b_digitized.csv`).
- Figure 5 (digitized): benzene solubility panel using digitized data/model points from paper panel (`figure_5/plot_figure_5.py`, `figure_5_digitized.csv`).
- Figure 6a/6b: both panels now read the provided `figure_6/1-butanol-NH4Cl-water-LLE.csv`; the data points are reused in the ternary panel and the model tie lines are generated with the built-in LLE solver using `held_2014` Tables 1-4 assumptions plus the constant-rule, salt-free solvent-only dielectric blend from Held 2014 Eq. (23) (`figure_6/_shared.py`, `figure_6/plot_figure_6a.py`, `figure_6/plot_figure_6b.py`).
- Figure 6b (digitized): ternary phase diagram reproduced from digitized tie-line endpoints and envelope points (`figure_6/plot_figure_6b.py`, `figure_6b_digitized.csv`).

Remaining gap:

- Figures 4b and 6b are currently digitized reproductions, not fresh ePC-SAFT equilibrium solves from raw thermodynamic constraints.

Optional digitized-data hook for Figure 2:

- `figure_2/figure_2_digitized.csv` with columns:
  `salt,temperature_K,molality,osmotic`
- Supported `salt` values: `NaCl`, `KCl`.
- Supported `temperature_K` values: `273.15`, `298.15`.

