# Figure 6 feos Analysis Notes

## Source provenance

- This Figure 6 `feos` run is not independent of this repo's parameter set.
- The custom `feos` ePC-SAFT pure/binary JSON files are built directly from `get_prop_dict("2020_Bulow", ["Li+", "Br-", "Ethanol"], ...)` in this repo.
- That means the parameter source for the current Figure 6 `feos` comparison is the `PC-SAFT` repo's `2020_Bulow` values, not stock `feos` parameter files.
- The comparison is still useful for checking implementation behavior and bookkeeping across molality, but it is not an external parameter-set validation in the same sense as the water-only Figure 3 work.

## Outputs

- `figure_6a_feos_vs_pcsaft.png` compares Bulow-2020 experimental points, the current repo total curve, and the `feos` total curve.
- `figure_6b_feos_vs_pcsaft_mu.png` compares the current repo Figure 6b-style $\mu$-basis contributions against `feos` $\mu$-basis contributions.
- `figure_6b_feos_bookkeeping.png` shows how the `feos` total compares with the summed `feos` $\mu$ and reconstructed per-term $\ln\varphi$ contributions across molality.
- `figure6_feos_vs_pcsaft_stats.csv` reports per-series deltas between the `PC-SAFT` and `feos` curves.

## Important caveat

- `feos` Born uses the hard-sphere diameter path in its ePC-SAFT implementation rather than the repo's explicit `d_born` values, so the Born comparison is not a strict apples-to-apples port even though the other pure/binary values are ported from `2020_Bulow`.

## Fit summary

- Repo total RMSE vs Bulow-2020 data: 1.0530
- feos total RMSE vs Bulow-2020 data: 1.0530
- Total mean |feos - pcsaft|: 0.000027
- Total max |feos - pcsaft|: 0.000028
- Max feos closure gap |total - mu_sum|: 0.2368
- Max feos closure gap |total - lnfug_sum|: 0.0000
