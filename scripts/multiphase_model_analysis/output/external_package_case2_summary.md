# External Package Ascani Case-2 Comparison

## Worked example

- System: Water + 1-butanol + NaCl + KCl worked example
- Feed mass fractions: water=0.8094, 1-butanol=0.1728, NaCl=0.0054, KCl=0.0124
- Clapeyron no-charge build success: False
- Clapeyron explicit-charge build success: True
- Clapeyron flash method: `MichelsenTPFlash lle K0_A`
- Clapeyron phase betas: organic=0.050865, aqueous=0.949135
- Clapeyron max neutral ln fugacity gap: 1.729e-11
- Clapeyron max mean-ionic ln fugacity gap: 1.728e-11
- feos advanced build success: False
- feos revised build success: False

## Worked example vs paper

| Quantity | Paper | Package | |Δ| |
|---|---:|---:|---:|
| x_water_org | 0.4426 | 0.00750175 | 0.435098 |
| x_butanol_org | 0.557 | 0.9754 | 0.4184 |
| x_nacl_org | 4.15e-05 | 0.000317229 | 0.000275729 |
| x_kcl_org | 0.00042 | 0.016781 | 0.016361 |
| x_water_aq | 0.9627 | 0.995156 | 0.0324562 |
| x_butanol_aq | 0.0122 | 1.59597e-05 | 0.012184 |
| x_nacl_aq | 0.0076 | 0.00203071 | 0.00556929 |
| x_kcl_aq | 0.0174 | 0.00279715 | 0.0146028 |
| lnf_water_bar | -3.521 | -3.4669 | 0.0541045 |
| lnf_butanol_bar | -5.088 | -4.71403 | 0.373974 |

## Capability audit

| System | Clapeyron build no charge | Clapeyron build explicit charge | feos advanced | feos revised |
|---|---:|---:|---:|---:|
| Water + 1-butanol + NaCl + KCl worked example | False | True | False | False |
| Water + 1-butanol + KCl + NH4Cl | False | True | False | False |
| Water + 1-propanol + NaCl + KCl | True | True | False | False |
| Water + 1-propanol + KCl + NH4Cl | True | True | False | False |

## Notes

- Clapeyron alcohol-ion unlike rows for worked example: 0.
- feos alcohol permittivity present for worked example alcohol: False.
- feos water-alcohol binary present in shipped PC-SAFT data: True.
- feos alcohol-ion binary rows for worked example: False.
- Clapeyron mean-ionic ln fugacities are captured in the raw JSON, but the package-accessible ionic reference basis does not appear to match the paper's tabulated mean-ionic fugacity basis closely enough for a direct numeric comparison.
