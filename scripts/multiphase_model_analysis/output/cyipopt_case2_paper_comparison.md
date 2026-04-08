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
| $x_{water}^{(org)}$ | 0.4426 | 0.523349 | 0.0807492 | 0.935859 | 0.493259 |
| $x_{butanol}^{(org)}$ | 0.557 | 0.476303 | 0.080697 | 0.0535634 | 0.503437 |
| $x_{NaCl}^{(org)}$ | 4.15e-05 | 2.46274e-05 | 1.68726e-05 | 0.00186347 | 0.00182197 |
| $x_{KCl}^{(org)}$ | 0.00042 | 0.000149121 | 0.000270879 | 0.00342552 | 0.00300552 |
| $x_{water}^{(aq)}$ | 0.9627 | 0.964705 | 0.00200469 | 0.940396 | 0.0223039 |
| $x_{butanol}^{(aq)}$ | 0.0122 | 0.0238532 | 0.0116532 | 0.0487721 | 0.0365721 |
| $x_{NaCl}^{(aq)}$ | 0.0076 | 0.00204533 | 0.00555467 | 0.00193429 | 0.00566571 |
| $x_{KCl}^{(aq)}$ | 0.0174 | 0.00367574 | 0.0137243 | 0.00348161 | 0.0139184 |
| $\ln(f_{water}/bar)$ | -3.521 | -3.51447 | 0.00653239 | -3.49172 | 0.0292827 |
| $\ln(f_{butanol}/bar)$ | -5.088 | -5.24278 | 0.154777 | -5.09184 | 0.003845 |
| $\ln(f_{\pm,NaCl}/bar)$ | -224.891 | -226.379 | 1.48776 | -226.289 | 1.398 |
| $\ln(f_{\pm,KCl}/bar)$ | -206.733 | -208.363 | 1.62959 | -208.331 | 1.59832 |
| $\hat g_{feed}$ (J/mol) | -27361.3 | -14521.2 | 12840.1 | -14521.2 | 12840.1 |
| $\hat g_{eq}$ (J/mol) | -27479.9 | -14532 | 12947.9 | -14521.2 | 12958.7 |
| $\Delta\hat g$ (J/mol) | -118.543 | -10.8381 | 107.705 | 1.70353e-05 | 118.543 |
| $\beta_{org}$ |  | 0.0551289 |  | 0.00503598 |  |
| $\beta_{aq}$ |  | 0.944871 |  | 0.994964 |  |
| $TPDF_{min}$ |  | -0.10144 |  | -0.10144 |  |
| Residual norm |  | 0.0565613 |  | 0.000218797 |  |
| Charge residual org |  | -2.72731e-07 |  | 1.12409e-10 |  |
| Charge residual aq |  | 1.59126e-08 |  | -5.68957e-13 |  |
| Mass balance max error |  | 2.1684e-19 |  | 0 |  |
| Neutral fugacity gap max |  | 0.0565402 |  | 0.000180885 |  |
| Mean-ionic gap max |  | 5.00744e-05 |  | 0.000122181 |  |
| $\eta_{Na^+\to aq}$ (%) |  | 99.9298 |  | 99.5148 |  |
| $\eta_{K^+\to aq}$ (%) |  | 99.7639 |  | 99.5045 |  |
| $\eta_{Cl^-\to aq}$ (%) |  | 99.8228 |  | 99.5081 |  |

### Solver diagnostics

| Backend | Converged | Status | Message | Nit | Objective | Residual norm |
|---|---:|---:|---|---:|---:|---:|
| Current least_squares solver | True | 4 | Both `ftol` and `xtol` termination conditions are satisfied. |  |  | 0.0565613 |
| cyipopt IPOPT experiment | True | 0 | b'Algorithm terminated successfully at a locally optimal point, satisfying the convergence tolerances (can be specified by options).' | 37 | 2.39361e-08 | 0.000218797 |

### Phase-resolved implementation values

#### Current least_squares solver

**Organic-rich phase**

- Phase fraction $\beta$: 0.0551289
| Species | Mole fraction | Share of feed (%) | $\ln(f/bar)$ |
|---|---:|---:|---:|
| H2O | 0.523349 | 3.06811 | -3.54274 |
| Butanol | 0.476303 | 53.8117 | -5.24355 |
| Na+ | 2.46274e-05 | 0.0702032 | -225.234 |
| K+ | 0.000149121 | 0.236143 | -189.202 |
| Cl- | 0.000174021 | 0.177159 | -227.523 |
| NaCl_pm |  |  | -226.379 |
| KCl_pm |  |  | -208.363 |

**Aqueous-rich phase**

- Phase fraction $\beta$: 0.944871
| Species | Mole fraction | Share of feed (%) | $\ln(f/bar)$ |
|---|---:|---:|---:|
| H2O | 0.964705 | 96.9319 | -3.4862 |
| Butanol | 0.0238532 | 46.1883 | -5.24201 |
| Na+ | 0.00204533 | 99.9298 | -224.961 |
| K+ | 0.00367574 | 99.7639 | -188.929 |
| Cl- | 0.00572106 | 99.8228 | -227.796 |
| NaCl_pm |  |  | -226.379 |
| KCl_pm |  |  | -208.363 |

#### cyipopt IPOPT experiment

**Organic-rich phase**

- Phase fraction $\beta$: 0.00503598
| Species | Mole fraction | Share of feed (%) | $\ln(f/bar)$ |
|---|---:|---:|---:|
| H2O | 0.935859 | 0.50118 | -3.49171 |
| Butanol | 0.0535634 | 0.552797 | -5.09194 |
| Na+ | 0.00186347 | 0.485248 | -224.846 |
| K+ | 0.00342552 | 0.495525 | -188.93 |
| Cl- | 0.00528898 | 0.491855 | -227.732 |
| NaCl_pm |  |  | -226.289 |
| KCl_pm |  |  | -208.331 |

**Aqueous-rich phase**

- Phase fraction $\beta$: 0.994964
| Species | Mole fraction | Share of feed (%) | $\ln(f/bar)$ |
|---|---:|---:|---:|
| H2O | 0.940396 | 99.4988 | -3.49172 |
| Butanol | 0.0487721 | 99.4472 | -5.09175 |
| Na+ | 0.00193429 | 99.5148 | -224.844 |
| K+ | 0.00348161 | 99.5045 | -188.928 |
| Cl- | 0.0054159 | 99.5081 | -227.734 |
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
| $x_{water}^{(org)}$ | 0.4426 | 0.537917 | 0.0953175 | 0.935982 | 0.493382 |
| $x_{butanol}^{(org)}$ | 0.557 | 0.437547 | 0.119453 | 0.0529901 | 0.50401 |
| $x_{NaCl}^{(org)}$ | 4.15e-05 | 0.00301459 | 0.00297309 | 0.00195325 | 0.00191175 |
| $x_{KCl}^{(org)}$ | 0.00042 | 0.00925329 | 0.00883329 | 0.00356095 | 0.00314095 |
| $x_{water}^{(aq)}$ | 0.9627 | 0.964539 | 0.00183913 | 0.940396 | 0.0223042 |
| $x_{butanol}^{(aq)}$ | 0.0122 | 0.0254533 | 0.0132533 | 0.0487747 | 0.0365747 |
| $x_{NaCl}^{(aq)}$ | 0.0076 | 0.00186904 | 0.00573096 | 0.00193383 | 0.00566617 |
| $x_{KCl}^{(aq)}$ | 0.0174 | 0.00313474 | 0.0142653 | 0.00348092 | 0.0139191 |
| $\ln(f_{water}/bar)$ | -3.521 | -3.50686 | 0.0141411 | -3.48714 | 0.033857 |
| $\ln(f_{butanol}/bar)$ | -5.088 | -5.30515 | 0.217146 | -5.17868 | 0.0906806 |
| $\ln(f_{\pm,NaCl}/bar)$ | -224.891 | -143.468 | 81.4228 | -143.492 | 81.3991 |
| $\ln(f_{\pm,KCl}/bar)$ | -206.733 | -129.628 | 77.1049 | -129.656 | 77.0766 |
| $\hat g_{feed}$ (J/mol) | -27361.3 | -12369.2 | 14992.1 | -12369.2 | 14992.1 |
| $\hat g_{eq}$ (J/mol) | -27479.9 | -12378 | 15101.8 | -12369.2 | 15110.7 |
| $\Delta\hat g$ (J/mol) | -118.543 | -8.84487 | 109.698 | 2.58876e-05 | 118.543 |
| $\beta_{org}$ |  | 0.0566448 |  | 0.00510471 |  |
| $\beta_{aq}$ |  | 0.943355 |  | 0.994895 |  |
| $TPDF_{min}$ |  | -0.0800311 |  | -0.0800311 |  |
| Residual norm |  | 0.0474452 |  | 0.00038465 |  |
| Charge residual org |  | 4.928e-08 |  | 4.0274e-10 |  |
| Charge residual aq |  | -2.95907e-09 |  | -2.06642e-12 |  |
| Mass balance max error |  | 1.11022e-16 |  | 0 |  |
| Neutral fugacity gap max |  | 0.0474282 |  | 0.000358964 |  |
| Mean-ionic gap max |  | 0.000237244 |  | 0.000126174 |  |
| $\eta_{Na^+\to aq}$ (%) |  | 91.1703 |  | 99.4844 |  |
| $\eta_{K^+\to aq}$ (%) |  | 84.9439 |  | 99.4779 |  |
| $\eta_{Cl^-\to aq}$ (%) |  | 87.1676 |  | 99.4802 |  |

### Solver diagnostics

| Backend | Converged | Status | Message | Nit | Objective | Residual norm |
|---|---:|---:|---|---:|---:|---:|
| Current least_squares solver | True | 3 | `xtol` termination condition is satisfied. |  |  | 0.0474452 |
| cyipopt IPOPT experiment | True | 0 | b'Algorithm terminated successfully at a locally optimal point, satisfying the convergence tolerances (can be specified by options).' | 19 | 7.39779e-08 | 0.00038465 |

### Phase-resolved implementation values

#### Current least_squares solver

**Organic-rich phase**

- Phase fraction $\beta$: 0.0566448
| Species | Mole fraction | Share of feed (%) | $\ln(f/bar)$ |
|---|---:|---:|---:|
| H2O | 0.537917 | 3.24023 | -3.53057 |
| Butanol | 0.437547 | 50.7923 | -5.30577 |
| Na+ | 0.00301459 | 8.82971 | -156.968 |
| K+ | 0.00925329 | 15.0561 | -129.288 |
| Cl- | 0.0122678 | 12.8324 | -129.968 |
| NaCl_pm |  |  | -143.468 |
| KCl_pm |  |  | -129.628 |

**Aqueous-rich phase**

- Phase fraction $\beta$: 0.943355
| Species | Mole fraction | Share of feed (%) | $\ln(f/bar)$ |
|---|---:|---:|---:|
| H2O | 0.964539 | 96.7598 | -3.48314 |
| Butanol | 0.0254533 | 49.2077 | -5.30452 |
| Na+ | 0.00186904 | 91.1703 | -157.077 |
| K+ | 0.00313474 | 84.9439 | -129.397 |
| Cl- | 0.00500379 | 87.1676 | -129.859 |
| NaCl_pm |  |  | -143.468 |
| KCl_pm |  |  | -129.628 |

#### cyipopt IPOPT experiment

**Organic-rich phase**

- Phase fraction $\beta$: 0.00510471
| Species | Mole fraction | Share of feed (%) | $\ln(f/bar)$ |
|---|---:|---:|---:|
| H2O | 0.935982 | 0.508087 | -3.48715 |
| Butanol | 0.0529901 | 0.554344 | -5.1785 |
| Na+ | 0.00195325 | 0.51557 | -157.089 |
| K+ | 0.00356095 | 0.522146 | -129.418 |
| Cl- | 0.0055142 | 0.519798 | -129.895 |
| NaCl_pm |  |  | -143.492 |
| KCl_pm |  |  | -129.656 |

**Aqueous-rich phase**

- Phase fraction $\beta$: 0.994895
| Species | Mole fraction | Share of feed (%) | $\ln(f/bar)$ |
|---|---:|---:|---:|
| H2O | 0.940396 | 99.4919 | -3.48713 |
| Butanol | 0.0487747 | 99.4457 | -5.17886 |
| Na+ | 0.00193383 | 99.4844 | -157.09 |
| K+ | 0.00348092 | 99.4779 | -129.419 |
| Cl- | 0.00541475 | 99.4802 | -129.894 |
| NaCl_pm |  |  | -143.492 |
| KCl_pm |  |  | -129.656 |

### Notes

- The paper-comparable fugacity entries are the phase-averaged neutral and mean-ionic values already used in the repo’s original Ascani case-2 comparison.
- The per-phase $\ln(f/bar)$ rows are implementation diagnostics and are not directly tabulated in the paper target set.
