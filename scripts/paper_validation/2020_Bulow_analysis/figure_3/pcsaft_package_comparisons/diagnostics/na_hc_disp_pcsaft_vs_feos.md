# Na+ hard-chain / dispersion comparison: current repo vs feos

This note compares only the `Na+` Figure 3 water infinite-dilution state after the local `feos` branch-bookkeeping patch.

## Main outcome

- Current repo dispersion still has the known residual: $\mu^{disp} - \partial a^{disp}/\partial x_k = -2.890093761\ \mathrm{kJ\ mol^{-1}}$.
- Current repo hard-chain still satisfies $\mu^{hc} = \partial a^{hc}/\partial x_k$ for this state: residual `-1.421085472e-14` kJ/mol.
- Patched `feos` hard-chain now matches the current repo full $\tilde{\mu}^{hc}$ branch to roundoff: `feos - pcsaft mu = 4.914806340e-09` kJ/mol.
- Patched `feos` dispersion now matches the current repo full $\tilde{\mu}^{disp}$ branch to roundoff: `feos - pcsaft mu = -5.444082518e-09` kJ/mol.
- The patched `feos` public branches are still not the same as a fixed-density finite difference of the exposed Helmholtz branch values: hc `mu_public = 24.398547739` vs `fd(h=1e-6) = -11.756351741`, disp `mu_public = -28.733396753` vs `fd(h=1e-6) = 12.960404913` kJ/mol.

## Interpretation

- The non-differential branch pieces already matched closely across packages, so the actual fix target was the public branch-bookkeeping helper, not the underlying $a^\alpha$ or $Z^\alpha$ values.
- The local patch in `state/properties.rs::chemical_potential_contributions(...)` now returns the repo-style branch quantity $\tilde{\mu}^\alpha = a^\alpha + Z^\alpha + dadx^\alpha - \sum_j x_j dadx_j^\alpha$.
- That makes the patched `feos` branch values align with the current repo branch values while preserving the already-correct total hydration energy.
- The finite-difference mismatch is expected because the public helper is no longer returning a raw fixed-density derivative of `molar_helmholtz_energy_contributions(...)`; it is returning the assembled branch chemical-potential expression.

## Files

- CSV: `C:\Users\Tanner\Documents\git\PC-SAFT\scripts\paper_validation\2020_Bulow_analysis\figure_3\pcsaft_package_comparisons\diagnostics\na_hc_disp_pcsaft_vs_feos.csv`
- feos patched Python contribution path: `C:\Users\Tanner\Documents\git\feos\crates\feos-core\src\state\properties.rs`
- feos residual helper not used by Python here: `C:\Users\Tanner\Documents\git\feos\crates\feos-core\src\state\residual_properties.rs`
- feos hard-chain contribution source: `C:\Users\Tanner\Documents\git\feos\crates\feos\src\epcsaft\eos\hard_chain.rs`
- feos dispersion contribution source: `C:\Users\Tanner\Documents\git\feos\crates\feos\src\epcsaft\eos\dispersion.rs`
