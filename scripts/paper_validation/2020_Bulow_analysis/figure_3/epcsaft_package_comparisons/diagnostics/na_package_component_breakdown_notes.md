# Na+ Package Component Breakdown Notes

These CSV files compare the current `ePC-SAFT` package, `feos`, and `Clapeyron.jl` for the Na+/Cl-/water infinite-dilution state used in Figure 3 package comparisons.

## Column meanings

- `a_kj_mol`: contribution Helmholtz term $RT a^\alpha$.
- `z_kj_mol`: contribution compressibility term $RT Z^\alpha$.
- `target_mu_kj_mol`: package-exposed Na+ contribution on the $\mu^\alpha$ basis when available.
- `target_lnfug_kj_mol`: Na+ fugacity-style contribution reconstructed as $RT\left(\mu^\alpha - \frac{Z^\alpha}{Z-1}\ln Z\right)$.
- `simplex_fd_*`: finite-difference proxy for $RT\,\partial a^\alpha/\partial x_k$ at fixed $T,\rho$ or fixed $T,V$ depending on the package API.

## Exposure differences

- `epcsaft` exposes explicit analytical derivative pieces for `hc`, `disp`, and `assoc`, including the same structures documented in `equations_v2.tex` such as $\partial a^{hs}/\partial x_k$, $\partial g^{hs}/\partial x_k$, and the association-site derivatives solved from the linear system.
- `feos` exposes contribution-resolved residual Helmholtz energies, pressures, and chemical-potential contributions, but not the sub-derivatives behind hard-sphere, hard-chain, or dispersion. The CSV therefore records raw contribution labels plus finite-difference proxies.
- `Clapeyron.jl` exposes the per-contribution Helmholtz functions and obtains $\mu$-style quantities through `VT_molar_gradient`. It exposes real association site fractions $X$, association strengths $\Delta$, and the site matrix, but not named derivative subpieces like $dX/dx_k$ or $d\Delta/dx_k` through a public API.

## Files

- `hc`: rows for epcsaft, feos, Clapeyron.jl.
- `disp`: rows for epcsaft, feos, Clapeyron.jl.
- `assoc`: rows for epcsaft, feos, Clapeyron.jl.
- `dh`: rows for epcsaft, feos, Clapeyron.jl.
- `born`: rows for epcsaft, feos, Clapeyron.jl.

