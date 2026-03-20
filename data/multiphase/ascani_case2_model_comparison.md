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
| $x_{water}^{(org)}$ | 0.4426 | 0.526628 | 0.538042 |
| $x_{butanol}^{(org)}$ | 0.557 | 0.472748 | 0.437422 |
| $x_{NaCl}^{(org)}$ | 4.15e-05 | 4.77755e-05 | 0.00301444 |
| $x_{KCl}^{(org)}$ | 0.00042 | 0.000264154 | 0.00925323 |
| $x_{water}^{(aq)}$ | 0.9627 | 0.964159 | 0.964554 |
| $x_{butanol}^{(aq)}$ | 0.0122 | 0.0244241 | 0.0254395 |
| $x_{NaCl}^{(aq)}$ | 0.0076 | 0.00204236 | 0.00186899 |
| $x_{KCl}^{(aq)}$ | 0.0174 | 0.00366627 | 0.00313443 |
| $\ln(f_{water}/bar)$ | -3.521 | -3.51327 | -3.5068 |
| $\ln(f_{butanol}/bar)$ | -5.088 | -5.24583 | -5.30531 |
| $\ln(f_{\pm,NaCl}/bar)$ | -224.891 | -199.363 | -143.468 |
| $\ln(f_{\pm,KCl}/bar)$ | -206.733 | -183.419 | -129.628 |
| $\hat g_{feed}$ (J/mol) | -27361.3 | -13832.2 | -12369.2 |
| $\hat g_{eq}$ (J/mol) | -27479.9 | -13842.3 | -12378 |
| $\Delta\hat g=\hat g_{eq}-\hat g_{feed}$ (J/mol) | -118.543 | -10.1195 | -8.84882 |
| $\beta_{org}$ |  | 0.0543628 | 0.0566935 |
| $\beta_{aq}$ |  | 0.945637 | 0.943307 |
| $TPDF_{min}$ |  | -0.0960286 | -0.0800311 |
| Residual norm |  | 0.0544131 | 0.0473408 |
| Phase split favored ($\Delta\hat g<0$) |  | True | True |
| Water prefers organic ($x_{water}^{org}/x_{water}^{aq}>1$) |  | False | False |
| $x_{water}^{org}/x_{water}^{aq}$ |  | 0.546205 | 0.557815 |
| $\eta_{NaCl\to aq}$ (%) |  | 99.8657 | 91.1631 |
| $\eta_{KCl\to aq}$ (%) |  | 99.5875 | 84.9311 |
| $\eta_{water\to org}$ (%) |  | 3.04443 | 3.24376 |
| Charge residual org |  | -4.46147e-07 | 1.71846e-08 |
| Charge residual aq |  | 2.56481e-08 | -1.03281e-09 |
| Mass-balance max error |  | 0 | 0 |
| $k_{ij}(water,butanol)$ used |  | -0.0143439 | -0.0143439 |
| $l_{ij}(water,butanol)$ used |  | -0.0044 | -0.0044 |
| $k_{ij}^{hb}(water,butanol)$ used |  | 0.026 | 0.026 |

## Transfer interpretation

- Both models have $\Delta\hat g<0$, so phase split is thermodynamically favored at the specified feed.
- Water partition ratio $x_{water}^{org}/x_{water}^{aq}<1$ for both models, so water does not preferentially transfer to the organic-rich phase.
- At equilibrium, chemical potentials are equal across phases; no net transfer occurs after convergence.
