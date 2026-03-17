# Figure 2 Quantity Audit

## Key Findings

- The original Figure 2 replication bug was that it zeroed all $k_{ij}$ values, but the 2019 paper uses the Table 6 water-IL unlike parameters.
- For panel 2a, the raw API cation derivative is not the closest salt-basis match. The strongest current candidate is the mean-ionic salt-basis quantity $0.25\,(d a^{ion}/d x_{cat} + d a^{ion}/d x_{anion})$.
- For panel 2b, the paper scale around $14$ at low $x_{IL}$ is incompatible with the raw API $\ln\varphi_{cat}$. The closest current candidate is the concentration-referenced ionic quantity $\ln\varphi^{ion}_{cat} - \ln x_{IL} + \ln(\rho RT/P)$.
- Numerical vs analytical DH derivatives are essentially identical for these curves; the main mismatch is the plotted quantity definition, not numerical derivative error.

## Best Current Matches

- Panel 2a `epc`: best current quantity is `dadx_mean_salt` (`numerical`), RMSE `0.4782`, shifted RMSE `0.3667`.
- Panel 2a `orig_il`: best current quantity is `dadx_mean_salt` (`numerical`), RMSE `0.9483`, shifted RMSE `0.4726`.
- Panel 2a `orig_water`: best current quantity is `dadx_mean_salt` (`analytical`), RMSE `0.1270`, shifted RMSE `0.0567`.
- Panel 2b `epc`: best current quantity is `lnphi_mean_salt_stdil` (`analytical`), RMSE `1.7389`, shifted RMSE `1.2854`.
- Panel 2b `orig_il`: best current quantity is `lnphi_mean_salt_stdil` (`analytical`), RMSE `2.5275`, shifted RMSE `2.5034`.
- Panel 2b `orig_water`: best current quantity is `lnphi_mean_salt_stdil` (`analytical`), RMSE `1.5489`, shifted RMSE `1.2345`.
