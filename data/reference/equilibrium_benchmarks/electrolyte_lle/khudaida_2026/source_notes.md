# Khudaida 2026 Electrolyte LLE Benchmark Fixture

This fixture promotes the already curated Khudaida 2026 water + ethanol + isobutanol + NaCl tie-line data into a stable confidence-validation input.

- Experimental tie-line rows come from the paper-table data retained in `analyses/paper_validation/application/2026_khudaida/data/input/table_3_4_experimental_tielines.csv`.
- Feed rows come from the digitized feed-composition CSVs under `analyses/paper_validation/application/2026_khudaida/figures/figure_2` through `figure_7`.
- Paper Tables 5-7 and Supporting Information Tables S1-S4 are retained under `analyses/paper_validation/application/2026_khudaida/data/input/`.
- CSV rows use formula-basis NaCl mole fractions. The confidence runner expands NaCl to explicit Na+ and Cl- mole fractions before calling the public native-backed equilibrium API.
- Thresholds are report bands for V5 confidence validation, not strict scientific acceptance criteria.
