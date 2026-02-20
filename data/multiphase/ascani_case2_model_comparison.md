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
| $x_{water}^{(org)}$ | 0.4426 | 0.42103 | 0.92169 |
| $x_{butanol}^{(org)}$ | 0.557 | 0.573482 | 0.0694691 |
| $x_{NaCl}^{(org)}$ | 4.15e-05 | 5.83494e-05 | 0.0013582 |
| $x_{KCl}^{(org)}$ | 0.00042 | 0.00268541 | 0.00306204 |
| $x_{water}^{(aq)}$ | 0.9627 | 0.97921 | 0.940456 |
| $x_{butanol}^{(aq)}$ | 0.0122 | 0.00956013 | 0.0487043 |
| $x_{NaCl}^{(aq)}$ | 0.0076 | 0.00207419 | 0.00193649 |
| $x_{KCl}^{(aq)}$ | 0.0174 | 0.00354084 | 0.00348319 |
| $\ln(f_{water}/bar)$ | -3.521 | -3.52256 | -3.48421 |
| $\ln(f_{butanol}/bar)$ | -5.088 | -5.11895 | -4.62049 |
| $\ln(f_{\pm,NaCl}/bar)$ | -224.891 | -200.637 | -140.214 |
| $\ln(f_{\pm,KCl}/bar)$ | -206.733 | -183.315 | -126.07 |
| $\hat g_{feed}$ (J/mol) | -27361.3 | -13772.3 | -12203.3 |
| $\hat g_{eq}$ (J/mol) | -27479.9 | -13818.1 | -12203.3 |
| $\Delta\hat g=\hat g_{eq}-\hat g_{feed}$ (J/mol) | -118.543 | -45.8509 | -0.0151041 |
| $\beta_{org}$ |  | 0.0695772 | 0.00442784 |
| $\beta_{aq}$ |  | 0.930423 | 0.995572 |
| $TPDF_{min}$ |  | -0.320445 | -26.0346 |
| Residual norm |  | 0.0919713 | 0.129762 |
| Phase split favored ($\Delta\hat g<0$) |  | True | True |
| Water prefers organic ($x_{water}^{org}/x_{water}^{aq}>1$) |  | False | False |
| $x_{water}^{org}/x_{water}^{aq}$ |  | 0.42997 | 0.980046 |
| $\eta_{NaCl\to aq}$ (%) |  | 99.7901 | 99.689 |
| $\eta_{KCl\to aq}$ (%) |  | 94.633 | 99.6105 |
| $\eta_{water\to org}$ (%) |  | 3.11516 | 0.433986 |
| Charge residual org |  | 5.81176e-08 | -3.90936e-07 |
| Charge residual aq |  | -4.34605e-09 | 1.7387e-09 |
| Mass-balance max error |  | 0 | 0 |
| $k_{ij}(water,butanol)$ used |  | 0 | 0 |
| $l_{ij}(water,butanol)$ used |  | -0.0044 | -0.0044 |
| $k_{ij}^{hb}(water,butanol)$ used |  | 0.026 | 0.026 |

## Transfer interpretation

- Both models have $\Delta\hat g<0$, so phase split is thermodynamically favored at the specified feed.
- Water partition ratio $x_{water}^{org}/x_{water}^{aq}<1$ for both models, so water does not preferentially transfer to the organic-rich phase.
- At equilibrium, chemical potentials are equal across phases; no net transfer occurs after convergence.
