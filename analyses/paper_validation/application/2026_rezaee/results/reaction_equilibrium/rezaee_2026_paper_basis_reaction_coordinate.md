# Rezaee 2026 Paper-Basis Reaction-Coordinate Benchmark

## Purpose

This benchmark follows the Section 3.2 reaction-coordinate idea directly from the 2025 DOE basis: equal aqueous/organic masses, 1000 ppm Li in the aqueous feed, Na/Li mass ratio from Table 5, pH-derived H+/OH-, TOPO wt% in the organic phase, and published Rezaee reaction constants and organic PC-SAFT parameters. It is deliberately separate from the package phase-equilibrium smoke.

## Result

- Modes tested: `3`.
- Rows tested: `26`.
- Best mode by Li extraction absolute error: `2025_born_no_ssm_empirical`.
- Mean Li extraction absolute error: `35.8154` percentage points.
- Mean Na extraction absolute error: `11.7562` percentage points.
- Median calculated Li extraction: `3.47205e-10%`.
- Median calculated Na extraction: `7.66336e-11%`.
- Median initial OH/Li mole ratio from pH basis: `0.000694`.

## Interpretation

Literal pH-derived OH- stoichiometry makes Reaction A/B OH-limited by orders of magnitude. This benchmark therefore does not reproduce the reported extraction percentages with the published K values before any parameter refit.

The pH-derived hydroxide amount is far too small to support 30-50% lithium extraction if Eq. 14/15 is interpreted as a literal one-OH-minus stoichiometric reaction extent. This points to an apparent-activity/reference-state interpretation in the paper workflow, not a simple ePC-SAFT Born/dielectric option issue.

## Generated Files

- `data\processed\rezaee_2026_paper_basis_reaction_coordinate_rows.csv`
- `results\reaction_equilibrium\rezaee_2026_paper_basis_reaction_coordinate_summary.json`
