# Ascani Case 2 Dataset Comparison

## Model presets

- `ascani_2022`: dielc rule 3 (combined); rel_perm diff=0; Born=on; d_Born_mode=0; shell=False; sat=False
- Coverage: Full ascani_2022 dataset coverage for water/butanol/Na+/K+/Cl-.
- `figiel_2025_overlay`: dielc rule 4 (empirical); rel_perm diff=1; Born=on; d_Born_mode=3; shell=True; sat=True
- Coverage: figiel_2025 has no Butanol entry, so case 2 uses ascani_2022 as the butanol base and overlays figiel_2025 shared H2O/ion pure parameters and available k_ij values.

## Feed composition

| Symbol | Value |
|---|---:|
| $w_{water}$ | 0.8094 |
| $w_{butanol}$ | 0.1728 |
| $w_{NaCl}$ | 0.0054 |
| $w_{KCl}$ | 0.0124 |
| $z_{water}$ | 0.940373 |
| $z_{butanol}$ | 0.0487962 |
| $z_{Na^+}$ | 0.00193393 |
| $z_{K^+}$ | 0.00348132 |
| $z_{Cl^-}$ | 0.00541526 |

## Paper vs model results

| Quantity | Paper (Ascani 2022) | ascani_2022 | figiel_2025_overlay |
|---|---:|---:|---:|
| $x_{water}^{(org)}$ | 0.4426 | 0.526293 | 0.522535 |
| $x_{butanol}^{(org)}$ | 0.557 | 0.473084 | 0.473206 |
| $x_{NaCl}^{(org)}$ | 4.15e-05 | 4.7683e-05 | 0.000131934 |
| $x_{KCl}^{(org)}$ | 0.00042 | 0.00026372 | 0.00199727 |
| $x_{water}^{(aq)}$ | 0.9627 | 0.964137 | 0.966944 |
| $x_{butanol}^{(aq)}$ | 0.0122 | 0.0244471 | 0.0218072 |
| $x_{NaCl}^{(aq)}$ | 0.0076 | 0.00204218 | 0.00204852 |
| $x_{KCl}^{(aq)}$ | 0.0174 | 0.00366598 | 0.0035757 |
| $\ln(f_{water}/bar)$ | -3.521 | -3.51343 | -3.51538 |
| $\ln(f_{butanol}/bar)$ | -5.088 | -5.24549 | -5.24741 |
| $\ln(f_{\pm,NaCl}/bar)$ | -224.891 | -199.363 | -146.319 |
| $\ln(f_{\pm,KCl}/bar)$ | -206.733 | -183.419 | -131.481 |
| $\hat g_{feed}$ (J/mol) | -27361.3 | -13832.2 | -12424.8 |
| $\hat g_{eq}$ (J/mol) | -27479.9 | -13842.3 | -12437.8 |
| $\Delta\hat g$ (J/mol) | -118.543 | -10.1106 | -12.9344 |
| $\beta_{org}$ |  | 0.0542737 | 0.0597898 |
| $\beta_{aq}$ |  | 0.945726 | 0.94021 |
| $TPDF_{min}$ |  | -0.0960286 | -0.10934 |
| Residual norm |  | 0.0547122 | 0.0590082 |
| Charge residual org |  | 2.46141e-07 | 2.51843e-08 |
| Charge residual aq |  | -1.41256e-08 | -1.60152e-09 |
| Mass balance max error |  | 0 | 0 |
| $\eta_{water\to org}$ (%) |  | 3.0375 | 3.32233 |
| $\eta_{butanol\to org}$ (%) |  | 52.6188 | 57.9817 |
| $\eta_{Na^+\to aq}$ (%) |  | 99.8662 | 99.5921 |
| $\eta_{K^+\to aq}$ (%) |  | 99.5889 | 96.5698 |
| $\eta_{Cl^-\to aq}$ (%) |  | 99.6881 | 97.6492 |

## Phase-resolved notes

- The PNG figure compares paper and model bars for organic-phase composition, aqueous-phase composition, and equilibrium fugacity metrics.
- A fourth panel shows model-only ion partitioning to the aqueous phase because the paper does not report ion-resolved phase-share data.
