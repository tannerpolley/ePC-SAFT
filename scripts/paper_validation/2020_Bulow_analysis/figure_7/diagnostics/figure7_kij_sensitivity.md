# Figure 7 `k_ij` sensitivity

Metric basis: RMSE/bias are computed against the visible paper points in each panel. For Figure 7a, this excludes the off-scale Bulow point above the plotted `y_max = 4`, because that point dominates any full-range metric while not affecting the displayed curve comparison.

## Main findings

- A literal `Li-H2O` unlike parameter does not enter Figure 7a or Figure 7c directly, because those runtime species sets are `Li+/Br-/Ethanol` and `Li+/Cl-/Methanol`.
- The practical analog is the active Li-solvent pair:
  - Figure 7a: `Li-Ethanol`
  - Figure 7c: `Li-Methanol`
- Figure 7a is only mildly sensitive within the tested range. Both increasing `k_ij(Li,Br)` and making `k_ij(Li,Ethanol)` more negative reduce the low-bias a little, but neither change removes the mismatch cleanly.
- Figure 7c is much more sensitive to `k_ij(Li,Methanol)`. Making it more negative moves the advanced curve upward and largely removes the small low-bias over the visible range.

## Summary table

| Panel | Requested pair | Runtime pair | Baseline k | Best k in scan | Baseline RMSE | Best RMSE | Baseline bias | Best bias |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 7a | Li-Br | Li-Br | 0.591 | 1.100 | 0.8713 | 0.8180 | -0.4697 | -0.4372 |
| 7a | Li-H2O | Li-Ethanol | 0.000 | -0.500 | 0.8713 | 0.8035 | -0.4697 | -0.4295 |
| 7c | Li-H2O | Li-Methanol | 0.000 | -0.500 | 0.0795 | 0.0167 | -0.0480 | +0.0057 |

## Interpretation

- Figure 7a: the residual is not especially sensitive to `k_ij(Li,Br)` alone. That pair can help a little, but not enough to explain the whole gap by itself.
- Figure 7a: the Li-solvent pair is the stronger lever than Li-Br, but the improvement is still modest over the tested range.
- Figure 7c: the Li-solvent pair is a strong lever. If the methanol curve is the one you think is systematically a bit low, this pair is the first parameter I would suspect.
- Because the best points for all three scans sit at the edge of the tested range, the trend is monotonic over the scanned interval rather than showing a clean interior optimum.
