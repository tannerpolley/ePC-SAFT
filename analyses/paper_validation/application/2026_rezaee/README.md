# Rezaee 2025/2026 Paper Validation

This analysis validates the package-facing Rezaee lithium-extraction paper workflow from the local source transcriptions. It is not a surrogate-generation, PrOMMiS/IDAES, slide-deck, or downstream process-design workflow.

## Boundary

- Active objective: run the Rezaee 2025/2026 source-data checks and Section 3.2 reaction-coordinate replication with the local `epcsaft` package.
- Not claimed: a direct published-constant closure of the reported extraction percentages, a calibrated distribution surrogate, PrOMMiS/IDAES transfer variables, costing rows, or deck-ready process outputs.
- Source rows: the local machine-readable supporting-information table contains 26 designed-experiment equilibrium-composition rows. The 2026 paper text also says 36 equilibrium data points; this workflow treats that as a source-text mismatch unless a new source-backed table is supplied.

## Source Assets

Paper markdown used for source checking is stored under `docs/papers/md/`:

- `Rezaee et al. - 2025 - Application of Response Surface Methodology for Selective Extraction of Lithium Using a Hydrophobic DES.md`
- `Rezaee et al. - 2025 - Supporting information - Application of Response Surface Methodology for Selective Extraction of Lithium.md`
- `Rezaee et al. - 2026 - Thermodynamic modeling of lithium extraction from synthetic brine using deep eutectic solvents A PC.md`
- `Rezaee et al. - 2026 - Supplementary material - Thermodynamic modeling of lithium extraction from synthetic brine using deep eutectic solvents.md`

Machine-readable source inputs are in `analyses/paper_validation/application/2026_rezaee/data/input/`.

## Paper-Validation Commands

Run these from the repository root:

```powershell
uv run python analyses\paper_validation\application\2026_rezaee\scripts\run_all.py
```

Surrogate, downstream bridge, PrOMMiS/IDAES, and deck scripts are intentionally not included in this package-validation lane.

## Generated Evidence

The strict pre-surrogate evidence bundle is:

- `data/processed/rezaee_2025_extraction_target_summary.csv`
- `data/processed/rezaee_2025_extraction_equilibrium_summary.csv`
- `results/targets/rezaee_2025_extraction_target_summary.md`
- `data/processed/rezaee_2026_section32_basis_inference_rows.csv`
- `results/reaction_equilibrium/rezaee_2026_section32_basis_inference_summary.json`
- `results/reaction_equilibrium/rezaee_2026_section32_basis_inference.md`
- `data/processed/rezaee_2026_section32_equilibrium_replication_rows.csv`
- `results/reaction_equilibrium/rezaee_2026_section32_equilibrium_replication_summary.json`
- `results/reaction_equilibrium/rezaee_2026_section32_equilibrium_replication.md`

Additional guardrail outputs in `results/reaction_equilibrium/` record replay, convention-scan, ePC-SAFT option-scan, and paper-basis reaction-coordinate diagnostics.
The lane-level gate summary is `results/reaction_equilibrium/summary.json`.

## Current Result

The package can run the paper-validation workflow from this copied application folder. The validation does not support claiming a direct published-constant Rezaee Section 3.2 closure:

- The source-supported Eq. 14/15 activity variant has a large combined median absolute log residual.
- No simple sign, reciprocal constant, activity on/off, water/OH, H+/NH4+, TOPO, Born, or dielectric option scan closes the gap.
- The 26-row Section 3.2 replication runs and writes the pre-surrogate result table, but the direct source-aligned case predicts essentially zero Li extraction and about 100% Li-extraction AARD against the paper's reported 7.89% post-Table-9 benchmark.

The honest completion state is a reproducible paper-validation result with a documented source/reference-state gap, not a promoted direct reactive-LLE surrogate basis.
