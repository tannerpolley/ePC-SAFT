# Ascani Case 2 Dataset Comparison

## Validation basis

- Fixed pure-component and binary-interaction parameter basis: `2022_Ascani`.
- Only the electrolyte/runtime user options are swapped between the two runs.

## Model presets

- `ascani2022_params_bulow2020_opts`: params=`2022_Ascani`, options=`2020_Bulow`, dielc rule 1 (linear-mole); rel_perm diff=0; Born=on; d_Born_mode=1; shell=False; sat=False
- Coverage: Pure-component and binary-interaction parameters are fixed to 2022_Ascani. Runtime/electrolyte options come from the current 2020_Bulow user_options.json, with solvent-specific ion sigma/dispersion precomputes disabled because 2022_Ascani only provides pure/any_solvent.csv.
- `ascani2022_params_figiel2025_opts`: params=`2022_Ascani`, options=`2025_Figiel`, dielc rule 4 (empirical); rel_perm diff=1; Born=on; d_Born_mode=3; shell=True; sat=True
- Coverage: Pure-component and binary-interaction parameters are fixed to 2022_Ascani. Runtime/electrolyte options come from the current 2025_Figiel user_options.json, with solvent-specific ion sigma/dispersion precomputes disabled because 2022_Ascani only provides pure/any_solvent.csv.

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

| Quantity | Paper (Ascani 2022) | ascani2022_params_bulow2020_opts | ascani2022_params_figiel2025_opts |
|---|---:|---:|---:|
| $x_{water}^{(org)}$ | 0.4426 | 0.546018 | 0.538042 |
| $x_{butanol}^{(org)}$ | 0.557 | 0.448321 | 0.437422 |
| $x_{NaCl}^{(org)}$ | 4.15e-05 | 0.000597464 | 0.00301444 |
| $x_{KCl}^{(org)}$ | 0.00042 | 0.00223299 | 0.00925323 |
| $x_{water}^{(aq)}$ | 0.9627 | 0.962068 | 0.964554 |
| $x_{butanol}^{(aq)}$ | 0.0122 | 0.0268173 | 0.0254395 |
| $x_{NaCl}^{(aq)}$ | 0.0076 | 0.00200745 | 0.00186899 |
| $x_{KCl}^{(aq)}$ | 0.0174 | 0.00355 | 0.00313443 |
| $\ln(f_{water}/bar)$ | -3.521 | -3.50931 | -3.5068 |
| $\ln(f_{butanol}/bar)$ | -5.088 | -5.26352 | -5.30531 |
| $\ln(f_{\pm,NaCl}/bar)$ | -224.891 | -226.539 | -143.468 |
| $\ln(f_{\pm,KCl}/bar)$ | -206.733 | -208.524 | -129.628 |
| $\hat g_{feed}$ (J/mol) | -27361.3 | -14528.7 | -12369.2 |
| $\hat g_{eq}$ (J/mol) | -27479.9 | -14536.4 | -12378 |
| $\Delta\hat g$ (J/mol) | -118.543 | -7.75405 | -8.84882 |
| $\beta_{org}$ |  | 0.0521442 | 0.0566935 |
| $\beta_{aq}$ |  | 0.947856 | 0.943307 |
| $TPDF_{min}$ |  | -0.0761361 | -0.0800311 |
| Residual norm |  | 0.0477113 | 0.0473408 |
| Charge residual org |  | -3.52227e-08 | 1.71846e-08 |
| Charge residual aq |  | 1.9377e-09 | -1.03281e-09 |
| Mass balance max error |  | 0 | 0 |
| Neutral fugacity gap max |  | 0.0476896 | 0.0473257 |
| Mean-ionic gap max |  | 0.000292027 | 8.66139e-05 |
| E-matrix rank |  | 2 | 2 |
| Independent ionic-pair count |  | 2 | 2 |
| $\eta_{water\to org}$ (%) |  | 3.0277 | 3.24376 |
| $\eta_{butanol\to org}$ (%) |  | 47.9081 | 50.8215 |
| $\eta_{Na^+\to aq}$ (%) |  | 98.3891 | 91.1631 |
| $\eta_{K^+\to aq}$ (%) |  | 96.6554 | 84.9311 |
| $\eta_{Cl^-\to aq}$ (%) |  | 97.2745 | 87.1567 |

## Algorithm checks

- `Neutral fugacity gap max` is the maximum absolute phase-to-phase difference in $\ln(f_i/bar)$ for neutral species.
- `Mean-ionic gap max` is the maximum absolute residual in $E(\ln f^{org}-\ln f^{aq})$ for the charged-species system.
- `E-matrix rank` should match the number of independent ionic-pair equations for the Ascani 2022 construction.
