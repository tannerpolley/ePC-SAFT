# Neutral Pure-Component Regression Basis

This folder is the phase-1 data foundation for neutral-component regression of `m`, `s`, and `e`.

Current scope:

- species: `Methane`, `Ethane`, `Propane`
- target parameters from the workbook: `m`, `s`, `e`
- property data basis: saturation pressure and saturated liquid density

Files:

- `hydrocarbon_basis_workbook_reference.csv`
  - workbook reference targets for `m`, `s`, and `e`
- `methane_nist_saturation.csv`
- `ethane_nist_saturation.csv`
- `propane_nist_saturation.csv`

Notes:

- The saturation CSV files were collected from the NIST Chemistry WebBook fluid-property tables using saturation-property queries over explicit temperature grids.
- These are trusted NIST thermophysical-property values and are suitable as an initial neutral-only fitting basis.
- They are not a substitute for a literature-curated experimental regression dataset if the goal is paper-grade parameter validation.
- The workbook reference file is taken from:
  - `workbooks/PC-SAFT Calculations - Hydrocarbon Basis.xlsm`
  - sheet: `PC-SAFT Liquid`
  - species row: `B13:D13`
  - `m`: `B16:D16`
  - `s`: `B17:D17`
  - `e`: `B18:D18`
