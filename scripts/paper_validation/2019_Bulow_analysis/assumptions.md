# Bulow 2019 Assumptions

- The `bulow_2019` parameter set is built from Table 3 pure IL-ion parameters and Table 6 water-ion `k_ij` relationships in the local markdown paper copy.
- The current 2019 scripts use the general electrolyte LLE solver with each ionic liquid represented as `H2O + cation + anion`.
- `ePC-SAFT` in this folder means concentration-dependent dielectric constant in the Debye-Huckel term with `rel_perm.rule=1` and Born disabled.
- `original ePC-SAFT` curves are approximated with constant dielectric constants by forcing all species dielectric constants to either water-like (`78.09`) or IL-like (`11.0`) values and turning off concentration dependence.
- The paper's phase-dependent dielectric treatment is approximated by stitching the water-rich branch from the `78.09` run to the IL-rich branch from the `11.0` run.
- Experimental data points and bars are not yet fully digitized in this first pass, so several figures are model-only reproductions of the paper layout.
