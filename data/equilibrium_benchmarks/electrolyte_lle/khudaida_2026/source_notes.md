# Khudaida 2026 Electrolyte LLE Benchmark Fixture

This fixture promotes the already curated Khudaida 2026 water + ethanol + isobutanol + NaCl tie-line data into a stable confidence-validation input.

- Experimental tie-line rows come from the paper-table data curated in `scripts/paper_validation/2026_Khudaida_analysis/_common.py`.
- Feed rows come from the digitized feed-composition CSVs under `scripts/paper_validation/2026_Khudaida_analysis/figure_2` through `figure_7`.
- CSV rows use formula-basis NaCl mole fractions. The confidence runner expands NaCl to explicit Na+ and Cl- mole fractions before calling the public native-backed equilibrium API.
- Thresholds are report bands for V5 confidence validation, not strict scientific acceptance criteria.
