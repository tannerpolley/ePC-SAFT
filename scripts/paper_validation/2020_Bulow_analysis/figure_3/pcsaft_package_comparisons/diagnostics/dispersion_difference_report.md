# Dispersion and feos contribution audit

## 1. Fluoride is a parameter-loading problem in the current Clapeyron comparison path

- In the current shipped comparison load (`mixed_current`), `F-` dispersion is `-21.953601 kJ/mol`, versus `-17.630768 kJ/mol` for the current repo.
- The same mixed path loads water/sodium and water/fluoride unlike values from the revised table: `k(H2O,Na+) = 2.379990` and `k(H2O,F-) = -0.250000`.
- The current repo uses `k(H2O,Na+) = 0.004500` and `k(H2O,F-) = 0.000000` at the Figure 3 state.
- Loading only the advanced unlike table moves fluoride to `-17.489011 kJ/mol`, already close to the current repo.
- Loading the repo override file on top leaves fluoride unchanged at `-17.489011 kJ/mol` even though it restores `k(Na+,F-)` from `0.000000` to `0.665000`.
- That means the `Na/F` unlike pair is not what drives the Figure 3 fluoride discrepancy here. The large shift comes from fixing the water-ion rows, and the remaining `0.141757 kJ/mol` gap is small enough to attribute mainly to the remaining pure-model differences, not the unlike-parameter path.

## 2. Why fluoride stands out while the other ions move together

- `Cl-`, `Br-`, and `I-` are unchanged by the alternate Clapeyron unlike loads used in this audit, so their remaining average offset of `0.738694 kJ/mol` versus the current repo is not a k-parameter load-order issue. It is a baseline model or pure-parameter difference between the two frameworks.
- `Na+` moves only from `-28.306436` to `-28.184333 kJ/mol` when the water/sodium row is fixed, a change of `0.122103 kJ/mol`.
- `F-` moves from `-21.953601` to `-17.489011 kJ/mol`, a much larger change of `4.464590 kJ/mol`.
- That contrast isolates the fluoride anomaly to the water/fluoride unlike row. The other ions may still differ from the current repo, but not because the current comparison path is loading the wrong water-ion `k_ij` row for them.

## 3. feos does not expose the same branch definition as the current repo

- `feos` dispersion matches the current repo `dadx_disp` branch to roundoff for every audited ion (`max |feos - dadx_disp| = 9.111154e-09 kJ/mol`).
- `feos` hard-chain differs from the current repo `dadx_hc` by an ion-independent `38.071805 kJ/mol` offset across all seven ions.
- That pattern is what you would expect from an internal state-scalar mismatch, not from different physical parameters: the ion ranking is preserved, but the branch baseline is shifted.
- Source inspection of `feos` shows the public contribution path is not computing the same quantity as the repo's `mu^alpha = a^alpha + Z^alpha + dadx^alpha - \sum x_j dadx_j^alpha` formula. In `feos-core/src/state/residual_properties.rs`, `residual_chemical_potential_contributions(...)` differentiates `molar_helmholtz_energy_contributions(t, v, x)` directly with respect to composition, and the local variable for `v` is incorrectly initialized from reduced temperature instead of reduced inverse density.
- That explains the observed mix of behaviors: dispersion coming out as the repo `dadx` branch, association nearly matching because its state scalar is zero here, and hard-chain carrying a large constant offset because the density-sensitive part is evaluated at the wrong reduced variable.

## 4. Practical conclusion

- The current `Clapeyron` fluoride dispersion discrepancy is mostly a load-order issue, not evidence that its dispersion expression is fundamentally different.
- The `advanced_only` config is not fully available for every Figure 3 ion in Clapeyron's shipped database, so the diagnostic CSV records unsupported cases as missing rather than pretending that advanced-only is a complete water-ion set.
- The current `feos` `hc`/`disp` discrepancy is not a legitimate model-vs-model difference; it is a contribution-API mismatch, and the exposed branch values should not be treated as comparable `mu` contributions until that upstream path is fixed.
