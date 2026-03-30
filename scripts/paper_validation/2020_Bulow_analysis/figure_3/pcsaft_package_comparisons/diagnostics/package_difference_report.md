# Package Difference Report

This report compares the Figure 3 contribution bookkeeping for `Na+` and `Cl-` across the current repo `PC-SAFT`, `feos`, and `Clapeyron.jl` after correcting the package-side total check to use the fugacity-basis sum where it is reconstructible.

## Basis Checks

- `Na+`: `pcsaft total - adjusted_sum = +0.000000` kJ/mol, `feos total - adjusted_sum = -0.000000` kJ/mol, `clapeyron total - adjusted_sum = -0.000000` kJ/mol.
- `Cl-`: `pcsaft total - adjusted_sum = +0.000000` kJ/mol, `feos total - adjusted_sum = -0.000000` kJ/mol, `clapeyron total - adjusted_sum = +0.000000` kJ/mol.

The `pcsaft` and `Clapeyron` corrected sums close their totals to roundoff. `feos` does not: even after applying the same $-\frac{Z^\alpha}{Z-1}\ln Z$ correction using its pressure contributions, the reconstructed sum still misses the package total by about `-40.962` kJ/mol for both `Na+` and `Cl-`. That means the exposed `feos` contribution labels are not equivalent to the current repo's per-term $\mu^\alpha/Z^\alpha$ bookkeeping.

## Term-Level Observations

- `Na+` `hc` mu term: `pcsaft = +24.398548`, `feos = +24.398548`, `clapeyron = +30.719180` kJ/mol.
- `Na+` `disp` mu term: `pcsaft = -28.733397`, `feos = -28.733397`, `clapeyron = -28.306436` kJ/mol.
- `Na+` `assoc` mu term: `pcsaft = -3.215655`, `feos = -3.215655`, `clapeyron = -4.370860` kJ/mol.
- `Na+` `born` mu term: `pcsaft = -552.066713`, `feos = -552.095649`, `clapeyron = -408.223145` kJ/mol.
- `Cl-` `hc` mu term: `pcsaft = +23.287466`, `feos = +23.287466`, `clapeyron = +28.759551` kJ/mol.
- `Cl-` `disp` mu term: `pcsaft = -30.085203`, `feos = -30.085203`, `clapeyron = -29.684460` kJ/mol.
- `Cl-` `assoc` mu term: `pcsaft = -3.018089`, `feos = -3.018089`, `clapeyron = -4.007542` kJ/mol.
- `Cl-` `born` mu term: `pcsaft = -565.527846`, `feos = -565.557488`, `clapeyron = -354.060342` kJ/mol.

For `feos`, the association and Born `mu` terms stay very close to the current repo, while the hard-chain term is much more positive and the dispersion term shifts enough that the package-side branch sum no longer matches either `mu_total` or the corrected fugacity sum. That pattern points to a decomposition mismatch more than a total-EOS mismatch, because the overall hydration totals remain very close to the current repo.

For `Clapeyron`, the corrected fugacity sum closes exactly, so its remaining differences are genuine model/state differences rather than a missing basis conversion. Its Born term is much less negative than the current repo for both `Na+` and `Cl-`, and its association term is more negative than both the current repo and `feos`.

## Likely Causes

### 1. `feos` branch bookkeeping is not equivalent to the repo's `mu` decomposition

- `feos` uses `chemical_potential_contributions(...)` for branch terms, but the exposed branch set does not close to either its `mu_total` or its fugacity-basis total after the standard $Z$ correction.
- The reconstructed `feos` fugacity-basis sum still misses the package total by the same offset as the package `mu_total - mu_sum` gap. That makes it unlikely that the remaining discrepancy is just a missing $Z$ correction.
- `feos` Born uses the hard-sphere diameter directly in [born.rs](/Users/Tanner/Documents/git/feos/crates/feos/src/epcsaft/eos/born.rs), not an explicit `d_Born` table like the current repo or Clapeyron.

### 2. `Clapeyron.jl` is a genuinely different advanced-like model setup

- The extractor uses `ESElectrolyte + pharmaPCSAFT + DHBorn + LinMixRSP`, which is the closest composable analogue in-tree, but it is not the same implementation as the current repo's native `2020_Bulow` path.
- Current repo explicit Born diameters: `Na+ = 3.445` A, `Cl- = 4.100` A.
- Clapeyron explicit Born diameters from `born_like.csv`: `Na+ = 3.360` A, `Cl- = 3.874` A.
- Those Born-radius differences directly explain part of the Born-term gap. The remaining association/hc/disp differences are consistent with the different neutralmodel implementation and state point that follow from the composable Clapeyron setup.

### 3. Binary interaction parameters are not the main remaining cause for `feos`

- Current repo `2020_Bulow` uses `k_ij(H2O,Na+) = 0.004500`, `k_ij(H2O,Cl-) = -0.250000`, `k_ij(Na+,Cl-) = 0.317000`.
- `feos` Held-2014 binary JSON has `water/sodium ion k_ij = [0.0045, 0.0, 0.0, 0.0]`, `water/chloride ion k_ij = [-0.25, 0.0, 0.0, 0.0]`, `sodium ion/chloride ion k_ij = [0.317, 0.0, 0.0, 0.0]`.
- Clapeyron advanced-like unlike CSV has `water/sodium k = 0.004500`, `water/chloride k = -0.250000`, `sodium/chloride k = 0.317000`.
- So after fixing the missing binary files, the large residual `feos` branch discrepancy is still present even though the key aqueous `k` values now line up with the current repo. That reinforces that the remaining issue is the exposed contribution split, not simply missing unlike parameters.

## Output Artifacts

- Detailed cross-package CSV: `C:\Users\Tanner\Documents\git\PC-SAFT\scripts\paper_validation\2020_Bulow_analysis\figure_3\pcsaft_package_comparisons\diagnostics\figure3_package_mu_breakdown_comparison.csv`
- Main comparison totals/plots: `C:\Users\Tanner\Documents\git\PC-SAFT\scripts\paper_validation\2020_Bulow_analysis\figure_3\pcsaft_package_comparisons`
