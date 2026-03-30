# Ascani 2022 Case-2 ghat Audit

## Inputs

| Quantity | Value |
|---|---:|
| feed water mass fraction | 0.8094 |
| feed butanol mass fraction | 0.1728 |
| feed NaCl mass fraction | 0.0054 |
| feed KCl mass fraction | 0.0124 |
| paper ghat_feed_j_per_mol | -27361.317 |
| paper ghat_eq_j_per_mol | -27479.860 |

## Basis check

| Quantity | Value |
|---|---:|
| beta from water mass balance | 0.294751 |
| beta from butanol mass balance | 0.294787 |
| beta from NaCl mass balance | 0.291063 |
| beta from KCl mass balance | 0.294464 |
| least-squares beta if reported values are treated as formula-unit mole fractions | 0.0516339 |
| max residual under literal formula-unit mole-fraction interpretation | 0.013023 |

## ghat Reconstructions

| Reconstruction | ghat_eq (J/mol) | Delta vs paper (J/mol) |
|---|---:|---:|
| literal_stated_formula | -11743.153 | 15736.707 |
| mass_fraction_sum_once | -18591.256 | 8888.603859 |
| mass_fraction_sum_double_count_salts | -27938.472 | -458.612047 |
| feed_reference_formula_basis | -11749.325 |  |
| feed_reference_mass_basis_double_count_salts | -27974.795 |  |

## General Ionic Convention

| Quantity | Value |
|---|---:|
| ion-basis ghat_eq | -14535.781 |
| ion-basis ghat_feed | -14547.974 |
| ion-basis delta (eq-feed) | 12.192525 |
| ion-basis delta vs paper eq | 12944.079 |

## Conclusion

- The printed case-2 phase compositions behave like mass fractions, not mole fractions, and the printed ghat cannot be recovered from the printed numbers by the literal stated formula.
- Closest simple reconstruction: `mass_fraction_sum_double_count_salts` with delta -458.612047 J/mol.
