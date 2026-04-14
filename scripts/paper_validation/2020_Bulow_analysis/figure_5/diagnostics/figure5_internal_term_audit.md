# Figure 5 Internal Term Audit

This audit compares the paper Figure 5 transfer bars against the model-side transfer deltas of the internal contribution pieces exposed from `epcsaft_lnfugcoef_terms(...)`.

For each contribution $\alpha$, the report tracks the state-level pieces used in the code path:

- $a^\alpha$
- $Z^\alpha_{\mathrm{raw}}$
- $(\partial a^\alpha/\partial x_i)$ at the ion index
- $\sum_j x_j (\partial a^\alpha/\partial x_j)$
- the final code-level $\mu_i^\alpha$
- the final $\ln\varphi_i^\alpha$

Transfer deltas are formed as organic minus water at the same infinite-dilution reference states used for Figure 4 and Figure 5.

## Best-Matching Candidate By Solvent and Contribution

### Ethanol

- `hc`: best candidate is `reconstructed_mu` with RMSE `0.319 kJ/mol`, MAE `0.278 kJ/mol`, sign matches `3/3`.
- `disp`: best candidate is `dadx_minus_sum` with RMSE `0.805 kJ/mol`, MAE `0.776 kJ/mol`, sign matches `3/3`.
- `assoc`: best candidate is `reconstructed_mu` with RMSE `0.331 kJ/mol`, MAE `0.312 kJ/mol`, sign matches `3/3`.
- `dh`: best candidate is `neg_sum_x_dadx` with RMSE `0.000 kJ/mol`, MAE `0.000 kJ/mol`, sign matches `0/3`.
- `born`: best candidate is `dadx` with RMSE `5.087 kJ/mol`, MAE `5.055 kJ/mol`, sign matches `3/3`.

### Methanol

- `hc`: best candidate is `reconstructed_mu` with RMSE `2.213 kJ/mol`, MAE `1.474 kJ/mol`, sign matches `3/3`.
- `disp`: best candidate is `dadx` with RMSE `2.420 kJ/mol`, MAE `1.760 kJ/mol`, sign matches `3/3`.
- `assoc`: best candidate is `dadx_minus_sum` with RMSE `1.049 kJ/mol`, MAE `1.032 kJ/mol`, sign matches `3/3`.
- `dh`: best candidate is `neg_sum_x_dadx` with RMSE `0.000 kJ/mol`, MAE `0.000 kJ/mol`, sign matches `0/3`.
- `born`: best candidate is `mu` with RMSE `4.854 kJ/mol`, MAE `4.632 kJ/mol`, sign matches `3/3`.

## Notes

- `reconstructed_mu` is $\Delta a^\alpha + \Delta Z^\alpha_{raw} + \Delta(\partial a^\alpha/\partial x_i) - \Delta\sum_j x_j (\partial a^\alpha/\partial x_j)$.
- For `hc`, `disp`, `polar`, `dh`, and `born`, `reconstructed_mu` should match the exposed `mu` branch to numerical precision.
- The association branch is dumped exactly as the code currently assembles it. In the current implementation, the exposed `dadx_assoc` column is the same direct code-level quantity as `mu_assoc`, not a separate post-processed correction path.
- The detailed per-ion dump is in `figure5_internal_term_detail.csv`, and the candidate ranking table is in `figure5_internal_term_summary.csv`.

