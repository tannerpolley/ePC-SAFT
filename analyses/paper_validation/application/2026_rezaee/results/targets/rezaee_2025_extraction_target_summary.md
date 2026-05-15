# Rezaee 2025 Extraction Target Summary

## Scope

This report summarizes source-backed extraction-percentage targets that can be used to compare a Rezaee-style reactive extraction calculation against the paper.

The SI equilibrium-composition data are now tracked separately in `data/input/rezaee_2025_extraction_equilibrium_mole_fractions.csv` and summarized in `data/processed/rezaee_2025_extraction_equilibrium_summary.csv`.

## Target Counts

- Table 5 DOE Li/Na response rows: `26`.
- Table 7 optimum-neighborhood rows: `5`.
- Table 8 real-brine Li/Na comparison rows: `1` combined Li/Na target.
- SI Tables S1/S2 equilibrium mole-fraction rows: `26`.

## Selected Operating Point

- `T_C`: `23.0`.
- `pH`: `10.4`.
- `topo_wt_pct`: `10.0`.
- `na_li_mass_ratio`: `5.0`.
- `li_extraction_pct_exp`: `48.57`.
- `na_extraction_pct_exp_inferred`: `11.013605442176871`.
- `li_na_selectivity_exp`: `4.41`.

## Iran-Source Model Brine Check

- `li_extraction_pct_exp`: `51.63`.
- `na_extraction_pct_exp`: `9.97`.
- `li_na_selectivity_exp`: `5.1785356068204615`.

## Remaining Model Gate

The response and equilibrium-composition rows are enough to start a source-backed Li/Na fitting target. The remaining gap is model closure: connect the aqueous ions, organic DES/TOPO pseudo-components, and RLi/RNa complex species through a calibrated reactive/electrolyte LLE objective before claiming extraction-percent reproduction.
