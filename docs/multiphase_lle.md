# Multiphase Electrolyte LLE (V1)

`pcsaft_multiphase_lle` solves a two-liquid-phase electrolyte split at fixed `T`, `P`, and overall feed composition.

## API

```python
from pcsaft import pcsaft_multiphase_lle

result = pcsaft_multiphase_lle(
    t=298.15,
    p=1.0e5,
    z_feed=z_feed,          # species-level mole fractions (sum = 1)
    params=params,          # PC-SAFT/ePC-SAFT parameter dict
    species=species,        # aligned species labels
    options={
        "tpdf_global_trials": 4000,
        "tpdf_local_trials": 2000,
        "tpdf_tol": -1e-8,
        "solver_tol": 1e-9,
        "max_nfev": 200,
    },
)
```

## Notes

- Scope is liquid-liquid only in v1.
- Feed must be species-level and electroneutral.
- Ions are treated as fully dissociated species.
- The method uses:
  - Ascani-style E-matrix preprocessing for independent mean-ionic constraints.
  - TPDF random search as phase-stability seed.
  - Nonlinear least-squares solve for two-phase equilibrium.

## Experimental IPOPT Prototype

- A standalone `cyipopt` experiment lives at `scripts/multiphase_model_analysis/cyipopt_two_phase_experiment.py`.
- It reuses the same TPDF seed and EOS callbacks as `pcsaft_multiphase_lle`, but replaces the final `least_squares` stage with a `cyipopt.minimize_ipopt(...)` solve on the reduced two-phase residual system.
- The production package API is unchanged. This path is for side-by-side solver testing only.
- Manual benchmark entrypoint:

```bash
conda run -n PC-SAFT python scripts/multiphase_model_analysis/cyipopt_two_phase_experiment.py --bench-ascani-case2
```

## Output

The solver returns a structured dictionary with:

- `n_phases` (`1` or `2`)
- `phases` list (each phase has `beta`, `x`, `rho`, `lnfugcoef`, `lnfug`)
- `tpdf_min`, `tpdf_seed_x`
- `converged`, `status`, `message`, `residual_norm`
- `e_matrix`, `ion_pair_rows`

