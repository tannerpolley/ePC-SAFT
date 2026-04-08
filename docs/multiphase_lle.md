# Multiphase Electrolyte LLE (V1)

Archived reference note. The Python `PCSAFTState.multiphase_lle(...)` method has been removed from the supported object API and will be rewritten later in native code. This page is retained only for historical context.

## Notes

- Scope was liquid-liquid only in v1.
- Feed must be species-level and electroneutral.
- Ions were treated as fully dissociated species.
- The archived solver used:
  - Ascani-style E-matrix preprocessing for independent mean-ionic constraints.
  - TPDF random search as phase-stability seed.
  - Nonlinear least-squares solve for two-phase equilibrium.

## Output

The solver returns a structured dictionary with:

- `n_phases` (`1` or `2`)
- `phases` list (each phase has `beta`, `x`, `rho`, `lnfugcoef`, `lnfug`)
- `tpdf_min`, `tpdf_seed_x`
- `converged`, `status`, `message`, `residual_norm`
- `e_matrix`, `ion_pair_rows`

