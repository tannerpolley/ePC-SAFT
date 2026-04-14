# Figure 4-5 Accounting Audit

This audit compares the paper Figure 5 transfer contributions against model-side transfer values formed from exposed `mu`, `lnfug`, and `Z`-share terms at the same infinite-dilution states used for Figure 4.

## C++ Contribution Construction

- In [epcsaft_electrolyte.cpp](C:/Users/Tanner/Documents/git/ePC-SAFT/epcsaft_electrolyte.cpp), the hard-chain and dispersion branches are built as
  `mu_hc[i] = ares_hc + Zhc + dahc_dx[i] - sum_j x[j]*dahc_dx[j]`
  `mu_disp[i] = ares_disp + Zdisp + dadisp_dx[i] - sum_j x[j]*dadisp_dx[j]`
- The Debye-Huckel and Born branches use the same pattern with `a_DH` / `a_born`, `Z_DH` / `Zborn`, and the corresponding `dadx` vectors.
- The current public API already exposes the downstream `mu_*`, `lnfugcoef_*`, and `z_*` terms, but it does **not** yet expose the internal `a^alpha`, `dadx^alpha`, or `sum_x dadx^alpha` pieces separately. This audit therefore uses the fully exposed downstream pieces first.

## Main Findings

### Methanol

- `hc`: best exposed match is `mu` with RMSE `2.213` kJ/mol; `mu` RMSE is `2.213` and `lnfug` RMSE is `67.173`.
- `disp`: best exposed match is `mu` with RMSE `3.424` kJ/mol; `mu` RMSE is `3.424` and `lnfug` RMSE is `52.981`.
- `assoc`: best exposed match is `mu` with RMSE `1.133` kJ/mol; `mu` RMSE is `1.133` and `lnfug` RMSE is `16.198`.
- `born`: best exposed match is `mu` with RMSE `4.854` kJ/mol; `mu` RMSE is `4.854` and `lnfug` RMSE is `4.854`.

### Ethanol

- `hc`: best exposed match is `mu` with RMSE `0.319` kJ/mol; `mu` RMSE is `0.319` and `lnfug` RMSE is `18.675`.
- `disp`: best exposed match is `mu` with RMSE `2.063` kJ/mol; `mu` RMSE is `2.063` and `lnfug` RMSE is `6.405`.
- `assoc`: best exposed match is `mu` with RMSE `0.331` kJ/mol; `mu` RMSE is `0.331` and `lnfug` RMSE is `17.746`.
- `born`: best exposed match is `mu` with RMSE `5.087` kJ/mol; `mu` RMSE is `5.087` and `lnfug` RMSE is `5.087`.

## Candidate Ranking By Solvent and Term

| solvent | term | best candidate | RMSE | MAE | corr | sign matches |
| --- | --- | --- | --- | --- | --- | --- |
| methanol | hc | mu | 2.213 | 1.474 | 1.000 | 3 |
| methanol | disp | mu | 3.424 | 2.613 | 0.992 | 3 |
| methanol | assoc | mu | 1.133 | 1.118 | 0.999 | 3 |
| methanol | born | mu | 4.854 | 4.632 | 1.000 | 3 |
| ethanol | hc | mu | 0.319 | 0.278 | 1.000 | 3 |
| ethanol | disp | mu | 2.063 | 2.052 | 1.000 | 3 |
| ethanol | assoc | mu | 0.331 | 0.312 | 1.000 | 3 |
| ethanol | born | mu | 5.087 | 5.055 | 1.000 | 3 |

## Figure 4 vs Figure 5 Cross-Check

| solvent | ion | paper total (Fig 4) | paper total - born | model short sum lnfug | model short sum mu | model born lnfug |
| --- | --- | --- | --- | --- | --- | --- |
| methanol | Na+ | 10.473 | -4.876 | 0.484 | 2.486 | 9.759 |
| methanol | Cl- | 15.122 | -0.599 | 4.918 | 6.920 | 9.997 |
| methanol | I- | 15.563 | 5.465 | 8.002 | 10.005 | 7.513 |
| ethanol | Na+ | 15.099 | 5.188 | -0.391 | 2.545 | 15.316 |
| ethanol | Cl- | 19.775 | 9.596 | 4.142 | 7.078 | 15.689 |
| ethanol | I- | 19.042 | 11.500 | 7.236 | 10.172 | 11.791 |

## Interpretation

- For both alcohols, the paper `hc`, `disp`, and `assoc` bars are much closer to the transfer in raw $\Delta\tilde{\mu}^\alpha$ than to the transfer in $\Delta\ln\varphi^\alpha$.
- The Born transfer bar is different: it does **not** match the model well even on the raw-`mu` route, which means Figure 5 is not explained by a single bookkeeping switch the way Figure 3 largely is.
- Ethanol is the clearest case: the paper short-range bars are close to the raw `mu` transfer values, while the model totals only close when the correct $\ln\varphi$ transfer contributions are used. That is why Figure 4 totals can agree while Figure 5 component bars disagree strongly, including sign flips.
- The detailed CSV contains the water-side and organic-side state values for each branch, so you can inspect exactly where each sign change enters.

