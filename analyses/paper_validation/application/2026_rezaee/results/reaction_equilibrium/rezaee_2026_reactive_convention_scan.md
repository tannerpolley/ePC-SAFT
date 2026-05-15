# Rezaee 2026 Reactive Convention Scan

## Purpose

This diagnostic scans simple activity, stoichiometry, and reference-constant conventions for the Rezaee Li/Na reactive equilibrium lane. It is meant to prevent future agents from repeatedly trying the same one-line fixes without evidence.

## Source Boundary

- Machine-readable SI equilibrium rows available here: `26`.
- Benchmark equilibrium rows used here: `26`.
- Rezaee 2026 text mentions 36 equilibrium data points, but this workflow treats that as a clerical source-text mismatch because the 2025 SI Tables S1/S2 provide 26 designed-experiment equilibrium rows.
- The source-supported equation remains Rezaee Eq. 14/15 with activity coefficients.

## Source-Supported Result

- Variant: `paper_eq14_with_activity_vs_paper_k`.
- Li median ln residual: `32.5168`.
- Na median ln residual: `37.9726`.
- Li median absolute ln residual: `32.5168`.
- Na median absolute ln residual: `37.9726`.

## Best Numerical Simple Variant

- Variant: `paper_eq14_no_activity_vs_inverse_k`.
- Source supported: `False`.
- Combined median absolute ln residual: `9.50658`.
- Description: Mole-fraction quotient compared against reciprocal reported constants.

## Interpretation

No scanned simple convention closes the published-constant reactive equilibrium to a defensible tolerance for both Li and Na. The best numerical variant is not the equation stated in the paper, so it cannot be used as a published Rezaee reproduction claim.

The practical package-validation workflow remains: keep the ePC-SAFT activity and stability diagnostics in the validity layer, and treat direct published-constant closure as blocked until the missing source rows or reference-state convention are resolved.

## Generated Files

- `data\processed\rezaee_2026_reactive_convention_scan_rows.csv`
- `data\processed\rezaee_2026_reactive_convention_scan_summary.csv`
- `results\reaction_equilibrium\rezaee_2026_reactive_convention_scan_summary.json`
