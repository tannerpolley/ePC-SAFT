# Figure 2-3 Accounting Audit

This audit compares the paper Figure 3 hydration contributions against several model-side intermediate quantities at the same infinite-dilution state used for Figure 2.

## Main Findings

- The paper Figure 3 contributions do **not** sum to the paper Figure 2 advanced totals. The missing amount `paper total - (hc + disp + assoc + born)` ranges from `31.309` to `55.680` kJ/mol across the ions.
- The paper quantity `Figure 2 total - Figure 3 Born` matches the model short-range sum `hc + disp + assoc` on the correct $\ln\varphi$ basis very closely: RMSE = `0.734` kJ/mol.
- Individually, the paper `hc`, `disp`, and `assoc` bars match the model **raw** $\tilde{\mu}^\alpha$ terms far better than the correct $\ln\varphi^\alpha$ contributions. This is strongest evidence that the paper Figure 3 short-range bars were not plotted on the same bookkeeping basis as the Figure 2 totals.
- The Born bars are insensitive to that distinction here: Born RMSE is `0.494` kJ/mol for both `mu` and `lnfug`, which explains why Born appears consistent while the other contributions do not.

## Candidate Ranking By Term

| term | best candidate | RMSE | MAE | corr |
| --- | --- | --- | --- | --- |
| hc | mu | 15.776 | 14.014 | 0.775 |
| disp | mu | 5.587 | 5.538 | 0.999 |
| assoc | mu | 4.480 | 4.438 | 0.955 |
| born | mu | 0.494 | 0.362 | 1.000 |

## Mu vs Lnfug Comparison

| term | RMSE(mu) | RMSE(lnfug) | fit slope mu->paper | fit slope lnfug->paper |
| --- | --- | --- | --- | --- |
| hc | 15.776 | 165.966 | 0.206 | 0.206 |
| disp | 5.587 | 148.121 | 1.048 | 1.048 |
| assoc | 4.480 | 59.589 | 1.136 | 1.137 |
| born | 0.494 | 0.494 | 1.000 | 1.000 |

## Ion-Level Total Check

| ion | paper Fig2 total | paper Fig3 sum | missing to total | paper total - born | model short sum lnfug |
| --- | --- | --- | --- | --- | --- |
| Li+ | -555.716 | -595.661 | 39.944 | -7.113 | -7.170 |
| Na+ | -541.823 | -580.388 | 38.564 | 11.249 | 10.357 |
| K+ | -445.731 | -493.392 | 47.661 | 20.750 | 20.039 |
| F- | -869.465 | -900.774 | 31.309 | 10.664 | 10.067 |
| Cl- | -558.032 | -594.846 | 36.814 | 7.331 | 8.091 |
| Br- | -502.460 | -545.969 | 43.509 | 5.361 | 5.421 |
| I- | -423.734 | -479.414 | 55.680 | 1.406 | 2.644 |

## Interpretation

- The most plausible explanation is that Figure 2 uses the correct Gibbs-energy quantity $RT\ln\varphi_i^\infty$, while the non-Born Figure 3 bars were generated from a pre-compressibility intermediate closer to the raw $RT\tilde{\mu}_i^\alpha$ terms.
- Converting from $\tilde{\mu}_i^\alpha$ to $\ln\varphi_i^\alpha$ adds the explicit $Z^\alpha$ share correction from Eq. `lnphi_alpha`. For hydration in water, those short-range $Z$-share corrections are very large and compensating, so the individual `hc`, `disp`, and `assoc` bars swing dramatically while the sum still closes to the correct total.
- The paper Figure 3 short-range bars do not close even on the raw-`mu` basis, so there is likely a second issue as well: either digitization error in the small bars, rounding from the publication, or inconsistent internal plotting/export bookkeeping in the paper.
