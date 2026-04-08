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
| $x_{water}^{(org)}$ | 0.4426 | 0.526384 | 0.537917 |
| $x_{butanol}^{(org)}$ | 0.557 | 0.472993 | 0.437547 |
| $x_{NaCl}^{(org)}$ | 4.15e-05 | 4.77089e-05 | 0.00301459 |
| $x_{KCl}^{(org)}$ | 0.00042 | 0.00026381 | 0.00925329 |
| $x_{water}^{(aq)}$ | 0.9627 | 0.964142 | 0.964539 |
| $x_{butanol}^{(aq)}$ | 0.0122 | 0.0244414 | 0.0254533 |
| $x_{NaCl}^{(aq)}$ | 0.0076 | 0.00204223 | 0.00186904 |
| $x_{KCl}^{(aq)}$ | 0.0174 | 0.00366606 | 0.00313474 |
| $\ln(f_{water}/bar)$ | -3.521 | -3.51339 | -3.50686 |
| $\ln(f_{butanol}/bar)$ | -5.088 | -5.24558 | -5.30515 |
| $\ln(f_{\pm,NaCl}/bar)$ | -224.891 | -199.363 | -143.468 |
| $\ln(f_{\pm,KCl}/bar)$ | -206.733 | -183.419 | -129.628 |
| $\hat g_{feed}$ (J/mol) | -27361.3 | -13832.2 | -12369.2 |
| $\hat g_{eq}$ (J/mol) | -27479.9 | -13842.3 | -12378 |
| $\Delta\hat g=\hat g_{eq}-\hat g_{feed}$ (J/mol) | -118.543 | -10.1129 | -8.84487 |
| $\beta_{org}$ |  | 0.0542967 | 0.0566448 |
| $\beta_{aq}$ |  | 0.945703 | 0.943355 |
| $TPDF_{min}$ |  | -0.0960286 | -0.0800311 |
| Residual norm |  | 0.0546316 | 0.0474452 |
| Phase split favored ($\Delta\hat g<0$) |  | True | True |
| Water prefers organic ($x_{water}^{org}/x_{water}^{aq}>1$) |  | False | False |
| $x_{water}^{org}/x_{water}^{aq}$ |  | 0.545961 | 0.557694 |
| $\eta_{NaCl\to aq}$ (%) |  | 99.8661 | 91.1703 |
| $\eta_{KCl\to aq}$ (%) |  | 99.5885 | 84.9439 |
| $\eta_{water\to org}$ (%) |  | 3.03931 | 3.24023 |
| Charge residual org |  | 2.42232e-08 | 4.928e-08 |
| Charge residual aq |  | -1.39075e-09 | -2.95907e-09 |
| Mass-balance max error |  | 0 | 1.11022e-16 |
| $k_{ij}(water,butanol)$ used |  | -0.0143439 | -0.0143439 |
| $l_{ij}(water,butanol)$ used |  | -0.0044 | -0.0044 |
| $k_{ij}^{hb}(water,butanol)$ used |  | 0.026 | 0.026 |

## Transfer interpretation

- Both models have $\Delta\hat g<0$, so phase split is thermodynamically favored at the specified feed.
- Water partition ratio $x_{water}^{org}/x_{water}^{aq}<1$ for both models, so water does not preferentially transfer to the organic-rich phase.
- At equilibrium, chemical potentials are equal across phases; no net transfer occurs after convergence.
