# Ascani 2022 Case-2 Comparison

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

| Quantity | Paper (Ascani 2022) | Model 2020 | Model 2025 numeric |
|---|---:|---:|---:|
| $x_{water}^{(org)}$ | 0.4426 | 0.526473 | 0.407464 |
| $x_{butanol}^{(org)}$ | 0.557 | 0.472904 | 0.591114 |
| $x_{NaCl}^{(org)}$ | 4.15e-05 | 4.77021e-05 | 0.000115184 |
| $x_{KCl}^{(org)}$ | 0.00042 | 0.00026379 | 0.00059569 |
| $x_{water}^{(aq)}$ | 0.9627 | 0.964145 | 0.940404 |
| $x_{butanol}^{(aq)}$ | 0.0122 | 0.0244383 | 0.0487649 |
| $x_{NaCl}^{(aq)}$ | 0.0076 | 0.00204226 | 0.00193404 |
| $x_{KCl}^{(aq)}$ | 0.0174 | 0.00366612 | 0.00348149 |
| $\ln(f_{water}/bar)$ | -3.521 | -3.51334 | -3.58712 |
| $\ln(f_{butanol}/bar)$ | -5.088 | -5.24565 | -5.15283 |
| $\ln(f_{\pm,NaCl}/bar)$ | -224.891 | -199.363 | -137.68 |
| $\ln(f_{\pm,KCl}/bar)$ | -206.733 | -183.419 | -124.387 |
| $\hat g_{feed}$ (J/mol) | -27361.3 | -13832.2 | -12220.6 |
| $\hat g_{eq}$ (J/mol) | -27479.9 | -13842.3 | -12220.6 |
| $\Delta\hat g=\hat g_{eq}-\hat g_{feed}$ (J/mol) | -118.543 | -10.1148 | -0.0108222 |
| $\beta_{org}$ |  | 0.054314 | 5.7877e-05 |
| $\beta_{aq}$ |  | 0.945686 | 0.999942 |
| $TPDF_{min}$ |  | -0.0960286 | -35.897 |
| Residual norm |  | 0.0545518 | 0.212724 |
| Phase split favored ($\Delta\hat g<0$) |  | True | True |
| Water prefers organic ($x_{water}^{org}/x_{water}^{aq}>1$) |  | False | False |
| $x_{water}^{org}/x_{water}^{aq}$ |  | 0.546052 | 0.433286 |
| $\eta_{NaCl\to aq}$ (%) |  | 99.866 | 99.9997 |
| $\eta_{KCl\to aq}$ (%) |  | 99.5884 | 99.999 |
| $\eta_{water\to org}$ (%) |  | 3.0408 | 0.00250781 |
| Charge residual org |  | -5.14633e-07 | 4.09762e-09 |
| Charge residual aq |  | 2.95571e-08 | -2.37172e-13 |
| Mass-balance max error |  | 2.1684e-19 | 0 |
| $k_{ij}(water,butanol)$ used |  | -0.0143439 | -0.0143439 |
| $l_{ij}(water,butanol)$ used |  | -0.0044 | -0.0044 |
| $k_{ij}^{hb}(water,butanol)$ used |  | 0.026 | 0.026 |

## Transfer interpretation

- Both models have $\Delta\hat g<0$, so phase split is thermodynamically favored at the specified feed.
- Water partition ratio $x_{water}^{org}/x_{water}^{aq}<1$ for both models, so water does not preferentially transfer to the organic-rich phase.
- At equilibrium, chemical potentials are equal across phases; no net transfer occurs after convergence.
