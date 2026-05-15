# Rezaee Section 3.2 Initial-Basis Inference

## Purpose

This diagnostic tests whether the 2025 Table 5 extraction percentages and 2025 SI aqueous/organic equilibrium mole fractions define a self-consistent mole basis for the Rezaee 2026 Eq. 17 reaction-coordinate calculation.

The check does not fit parameters and does not use package equilibrium solvers. It only applies metal conservation to infer phase totals from Li, Na, RLi, and RNa rows.

## Result

- Rows evaluated: `26`.
- Median aqueous total consistency, Li-derived / Na-derived: `0.884441445349643`.
- Median organic total consistency, RLi-derived / RNa-derived: `3.6375026056361`.
- Max absolute aqueous charge residual: `8.000000000008868e-07`.
- Median SI xOH / pH-derived xOH estimate: `0.9999690361888166`.
- Median SI H moles / extracted Li+Na moles: `0.7407749435988276`.
- Median SI OH moles / extracted Li+Na moles: `0.0007005621825199862`.

## Interpretation

The SI aqueous rows are nearly charge balanced and the OH- mole fractions track the reported pH. However, combining Table 5 extraction percentages with SI organic RLi/RNa mole fractions does not produce one consistent organic phase total for one-metal RLi/RNa complexes. This is a source-basis blocker for exact Section 3.2 reproduction until the original initial-mole/reaction-coordinate worksheet or a documented phase-amount convention is available.

The direct replication should therefore keep using the 26 rows as the benchmark set, but it should not claim the paper's post-Table-9 AARD values until the hidden initial-mole/phase-amount convention is supplied or reconstructed from a source-backed worksheet.

## Generated Files

- `data\processed\rezaee_2026_section32_basis_inference_rows.csv`
- `results\reaction_equilibrium\rezaee_2026_section32_basis_inference_summary.json`
