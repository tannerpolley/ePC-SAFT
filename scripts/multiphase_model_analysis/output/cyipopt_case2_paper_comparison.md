# Ascani Case 2: cyipopt vs Current Solver vs Paper

## Validation basis

- Fixed pure-component and binary-interaction parameter basis: `2022_Ascani`.
- Two backends are compared for each option set: the current package `least_squares` solver and the experimental `cyipopt` IPOPT backend.
- Paper values are the case-2 targets already encoded in the repo and correspond only to the quantities tabulated in the paper basis.
- Phase-specific fugacity values below are implementation outputs; the paper comparison for fugacity uses the same phase-averaged basis as the existing case-study script.

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

## ascani2022_params_bulow2020_opts

- Option dataset: `2020_Bulow`
- Electrolyte settings: dielc rule 3 (combined); rel_perm diff=0; Born=on; d_Born_mode=1; shell=False; sat=False
- Coverage note: Pure-component and binary-interaction parameters are fixed to 2022_Ascani. Runtime/electrolyte options come from the current 2020_Bulow user_options.json, with solvent-specific ion sigma/dispersion precomputes disabled because 2022_Ascani only provides pure/any_solvent.csv.

### Paper-comparable summary

| Quantity | Paper | Current | |Error| current | cyipopt | |Error| cyipopt |
|---|---:|---:|---:|---:|---:|
| $x_{water}^{(org)}$ | 0.4426 | 0.523297 | 0.0806972 | 0.940046 | 0.497446 |
| $x_{butanol}^{(org)}$ | 0.557 | 0.476355 | 0.0806448 | 0.0491423 | 0.507858 |
| $x_{NaCl}^{(org)}$ | 4.15e-05 | 2.46402e-05 | 1.68598e-05 | 0.00192867 | 0.00188717 |
| $x_{KCl}^{(org)}$ | 0.00042 | 0.00014921 | 0.00027079 | 0.00347738 | 0.00305738 |
| $x_{water}^{(aq)}$ | 0.9627 | 0.964707 | 0.00200729 | 0.941023 | 0.0216767 |
| $x_{butanol}^{(aq)}$ | 0.0122 | 0.0238506 | 0.0116506 | 0.0481096 | 0.0359096 |
| $x_{NaCl}^{(aq)}$ | 0.0076 | 0.00204533 | 0.00555467 | 0.00194437 | 0.00565563 |
| $x_{KCl}^{(aq)}$ | 0.0174 | 0.00367573 | 0.0137243 | 0.00348916 | 0.0139108 |
| $\ln(f_{water}/bar)$ | -3.521 | -3.51449 | 0.00650858 | -3.49172 | 0.029282 |
| $\ln(f_{butanol}/bar)$ | -5.088 | -5.24277 | 0.154771 | -5.09185 | 0.00385214 |
| $\ln(f_{\pm,NaCl}/bar)$ | -224.891 | -226.379 | 1.48777 | -226.289 | 1.39799 |
| $\ln(f_{\pm,KCl}/bar)$ | -206.733 | -208.363 | 1.6296 | -208.331 | 1.59835 |
| $\hat g_{feed}$ (J/mol) | -27361.3 | -14521.2 | 12840.1 | -14521.2 | 12840.1 |
| $\hat g_{eq}$ (J/mol) | -27479.9 | -14532 | 12947.9 | -14521.2 | 12958.7 |
| $\Delta\hat g$ (J/mol) | -118.543 | -10.8376 | 107.705 | 0.000144413 | 118.543 |
| $\beta_{org}$ |  | 0.055128 |  | 0.664888 |  |
| $\beta_{aq}$ |  | 0.944872 |  | 0.335112 |  |
| $TPDF_{min}$ |  | -0.10144 |  | -0.10144 |  |
| Residual norm |  | 0.0566087 |  | 0.000516475 |  |
| Charge residual org |  | 1.04903e-07 |  | -2.17781e-11 |  |
| Charge residual aq |  | -6.12052e-09 |  | 4.32094e-11 |  |
| Mass balance max error |  | 0 |  | 1.11022e-16 |  |
| Neutral fugacity gap max |  | 0.0565902 |  | 0.000494359 |  |
| Mean-ionic gap max |  | 7.27746e-05 |  | 0.000147297 |  |
| $\eta_{Na^+\to aq}$ (%) |  | 99.9298 |  | 33.6922 |  |
| $\eta_{K^+\to aq}$ (%) |  | 99.7637 |  | 33.5866 |  |
| $\eta_{Cl^-\to aq}$ (%) |  | 99.8231 |  | 33.6243 |  |

### Solver diagnostics

| Backend | Converged | Status | Message | Nit | Objective | Residual norm |
|---|---:|---:|---|---:|---:|---:|
| Current least_squares solver | True | 3 | `xtol` termination condition is satisfied. |  |  | 0.0566087 |
| cyipopt IPOPT experiment | True | 0 | b'Algorithm terminated successfully at a locally optimal point, satisfying the convergence tolerances (can be specified by options).' | 37 | 1.33373e-07 | 0.000516475 |

### Phase-resolved implementation values

#### Current least_squares solver

**Organic-rich phase**

- Phase fraction $\beta$: 0.055128
| Species | Mole fraction | Share of feed (%) | $\ln(f/bar)$ |
|---|---:|---:|---:|
| H2O | 0.523297 | 3.06775 | -3.54279 |
| Butanol | 0.476355 | 53.8166 | -5.24349 |
| Na+ | 2.46402e-05 | 0.0702385 | -225.233 |
| K+ | 0.00014921 | 0.23628 | -189.201 |
| Cl- | 0.000173746 | 0.176875 | -227.524 |
| NaCl_pm |  |  | -226.379 |
| KCl_pm |  |  | -208.363 |

**Aqueous-rich phase**

- Phase fraction $\beta$: 0.944872
| Species | Mole fraction | Share of feed (%) | $\ln(f/bar)$ |
|---|---:|---:|---:|
| H2O | 0.964707 | 96.9322 | -3.4862 |
| Butanol | 0.0238506 | 46.1834 | -5.24205 |
| Na+ | 0.00204533 | 99.9298 | -224.961 |
| K+ | 0.00367573 | 99.7637 | -188.929 |
| Cl- | 0.00572107 | 99.8231 | -227.796 |
| NaCl_pm |  |  | -226.379 |
| KCl_pm |  |  | -208.363 |

#### cyipopt IPOPT experiment

**Organic-rich phase**

- Phase fraction $\beta$: 0.664888
| Species | Mole fraction | Share of feed (%) | $\ln(f/bar)$ |
|---|---:|---:|---:|
| H2O | 0.940046 | 66.4656 | -3.49173 |
| Butanol | 0.0491423 | 66.9603 | -5.0916 |
| Na+ | 0.00192867 | 66.3078 | -224.844 |
| K+ | 0.00347738 | 66.4134 | -188.928 |
| Cl- | 0.00540605 | 66.3757 | -227.734 |
| NaCl_pm |  |  | -226.289 |
| KCl_pm |  |  | -208.331 |

**Aqueous-rich phase**

- Phase fraction $\beta$: 0.335112
| Species | Mole fraction | Share of feed (%) | $\ln(f/bar)$ |
|---|---:|---:|---:|
| H2O | 0.941023 | 33.5344 | -3.49171 |
| Butanol | 0.0481096 | 33.0397 | -5.0921 |
| Na+ | 0.00194437 | 33.6922 | -224.843 |
| K+ | 0.00348916 | 33.5866 | -188.928 |
| Cl- | 0.00543353 | 33.6243 | -227.735 |
| NaCl_pm |  |  | -226.289 |
| KCl_pm |  |  | -208.331 |

### Notes

- The paper-comparable fugacity entries are the phase-averaged neutral and mean-ionic values already used in the repo’s original Ascani case-2 comparison.
- The per-phase $\ln(f/bar)$ rows are implementation diagnostics and are not directly tabulated in the paper target set.

## ascani2022_params_figiel2025_opts

- Option dataset: `2025_Figiel`
- Electrolyte settings: dielc rule 4 (empirical); rel_perm diff=1; Born=on; d_Born_mode=3; shell=True; sat=True
- Coverage note: Pure-component and binary-interaction parameters are fixed to 2022_Ascani. Runtime/electrolyte options come from the current 2025_Figiel user_options.json, with solvent-specific ion sigma/dispersion precomputes disabled because 2022_Ascani only provides pure/any_solvent.csv.

### Paper-comparable summary

| Quantity | Paper | Current | |Error| current | cyipopt | |Error| cyipopt |
|---|---:|---:|---:|---:|---:|
| $x_{water}^{(org)}$ | 0.4426 | 0.538042 | 0.0954423 | 0.608962 | 0.166362 |
| $x_{butanol}^{(org)}$ | 0.557 | 0.437422 | 0.119578 | 0.367839 | 0.189161 |
| $x_{NaCl}^{(org)}$ | 4.15e-05 | 0.00301444 | 0.00297294 | 0.00297206 | 0.00293056 |
| $x_{KCl}^{(org)}$ | 0.00042 | 0.00925323 | 0.00883323 | 0.00862742 | 0.00820742 |
| $x_{water}^{(aq)}$ | 0.9627 | 0.964554 | 0.00185364 | 0.968933 | 0.00623264 |
| $x_{butanol}^{(aq)}$ | 0.0122 | 0.0254395 | 0.0132395 | 0.0213027 | 0.0091027 |
| $x_{NaCl}^{(aq)}$ | 0.0076 | 0.00186899 | 0.00573101 | 0.00184447 | 0.00575553 |
| $x_{KCl}^{(aq)}$ | 0.0174 | 0.00313443 | 0.0142656 | 0.00303786 | 0.0143621 |
| $\ln(f_{water}/bar)$ | -3.521 | -3.5068 | 0.0141974 | -3.48157 | 0.0394313 |
| $\ln(f_{butanol}/bar)$ | -5.088 | -5.30531 | 0.21731 | -5.37455 | 0.286547 |
| $\ln(f_{\pm,NaCl}/bar)$ | -224.891 | -143.468 | 81.4228 | -143.471 | 81.4204 |
| $\ln(f_{\pm,KCl}/bar)$ | -206.733 | -129.628 | 77.1049 | -129.632 | 77.1014 |
| $\hat g_{feed}$ (J/mol) | -27361.3 | -12369.2 | 14992.1 | -12369.2 | 14992.1 |
| $\hat g_{eq}$ (J/mol) | -27479.9 | -12378 | 15101.8 | -12379 | 15100.9 |
| $\Delta\hat g$ (J/mol) | -118.543 | -8.84882 | 109.694 | -9.77901 | 108.764 |
| $\beta_{org}$ |  | 0.0566935 |  | 0.0793382 |  |
| $\beta_{aq}$ |  | 0.943307 |  | 0.920662 |  |
| $TPDF_{min}$ |  | -0.0800311 |  | -0.0800311 |  |
| Residual norm |  | 0.0473408 |  | 0.000315652 |  |
| Charge residual org |  | 1.71846e-08 |  | 3.7498e-10 |  |
| Charge residual aq |  | -1.03281e-09 |  | -3.2314e-11 |  |
| Mass balance max error |  | 0 |  | 0 |  |
| Neutral fugacity gap max |  | 0.0473257 |  | 0.000275446 |  |
| Mean-ionic gap max |  | 8.66139e-05 |  | 9.76282e-05 |  |
| $\eta_{Na^+\to aq}$ (%) |  | 91.1631 |  | 87.8073 |  |
| $\eta_{K^+\to aq}$ (%) |  | 84.9311 |  | 80.3384 |  |
| $\eta_{Cl^-\to aq}$ (%) |  | 87.1567 |  | 83.0058 |  |

### Solver diagnostics

| Backend | Converged | Status | Message | Nit | Objective | Residual norm |
|---|---:|---:|---|---:|---:|---:|
| Current least_squares solver | True | 3 | `xtol` termination condition is satisfied. |  |  | 0.0473408 |
| cyipopt IPOPT experiment | True | 0 | b'Algorithm terminated successfully at a locally optimal point, satisfying the convergence tolerances (can be specified by options).' | 33 | 4.98181e-08 | 0.000315652 |

### Phase-resolved implementation values

#### Current least_squares solver

**Organic-rich phase**

- Phase fraction $\beta$: 0.0566935
| Species | Mole fraction | Share of feed (%) | $\ln(f/bar)$ |
|---|---:|---:|---:|
| H2O | 0.538042 | 3.24376 | -3.53047 |
| Butanol | 0.437422 | 50.8215 | -5.3059 |
| Na+ | 0.00301444 | 8.83688 | -156.968 |
| K+ | 0.00925323 | 15.0689 | -129.288 |
| Cl- | 0.0122677 | 12.8433 | -129.968 |
| NaCl_pm |  |  | -143.468 |
| KCl_pm |  |  | -129.628 |

**Aqueous-rich phase**

- Phase fraction $\beta$: 0.943307
| Species | Mole fraction | Share of feed (%) | $\ln(f/bar)$ |
|---|---:|---:|---:|
| H2O | 0.964554 | 96.7562 | -3.48314 |
| Butanol | 0.0254395 | 49.1785 | -5.30471 |
| Na+ | 0.00186899 | 91.1631 | -157.077 |
| K+ | 0.00313443 | 84.9311 | -129.397 |
| Cl- | 0.00500342 | 87.1567 | -129.859 |
| NaCl_pm |  |  | -143.468 |
| KCl_pm |  |  | -129.628 |

#### cyipopt IPOPT experiment

**Organic-rich phase**

- Phase fraction $\beta$: 0.0793382
| Species | Mole fraction | Share of feed (%) | $\ln(f/bar)$ |
|---|---:|---:|---:|
| H2O | 0.608962 | 5.13774 | -3.48171 |
| Butanol | 0.367839 | 59.8072 | -5.37461 |
| Na+ | 0.00297206 | 12.1927 | -156.984 |
| K+ | 0.00862742 | 19.6616 | -129.306 |
| Cl- | 0.0115995 | 16.9942 | -129.957 |
| NaCl_pm |  |  | -143.471 |
| KCl_pm |  |  | -129.632 |

**Aqueous-rich phase**

- Phase fraction $\beta$: 0.920662
| Species | Mole fraction | Share of feed (%) | $\ln(f/bar)$ |
|---|---:|---:|---:|
| H2O | 0.968933 | 94.8623 | -3.48143 |
| Butanol | 0.0213027 | 40.1928 | -5.37449 |
| Na+ | 0.00184447 | 87.8073 | -157.08 |
| K+ | 0.00303786 | 80.3384 | -129.402 |
| Cl- | 0.00488233 | 83.0058 | -129.861 |
| NaCl_pm |  |  | -143.471 |
| KCl_pm |  |  | -129.632 |

### Notes

- The paper-comparable fugacity entries are the phase-averaged neutral and mean-ionic values already used in the repo’s original Ascani case-2 comparison.
- The per-phase $\ln(f/bar)$ rows are implementation diagnostics and are not directly tabulated in the paper target set.
