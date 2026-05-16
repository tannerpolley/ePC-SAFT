# Rezaee 2026 ePC-SAFT Parameter Smoke Report

Last updated: 2026-05-15

## Boundary

This is a diagnostic package/regression smoke test. Rezaee 2026 uses PC-SAFT for the organic DES/TOPO phase and ePC-SAFT for the aqueous electrolyte phase. It is not the flagship Shan/Gando HBTA/TOPO/sulfonated-kerosene parameterization.

The 2026 supporting information supplies the reaction constants, Gibbs-energy basis, density tables, organic PC-SAFT parameters, and organic binary interactions used here. The 2025 supporting information supplies the designed-experiment aqueous/organic equilibrium mole fractions used to start the Li/Na reactive extraction target.

## Density Fit

- Density source: `Rezaee2026_SI_TableS5` from `analyses/rezaee_2026_pcsaft_epcsaft/data/input/rezaee_2026_si_density_tables.csv`.
- DES molecular-weight basis: `0.20748` kg/mol.
- MW note: Average formula-unit molecular weight for TBAC + 2 decanoic acid divided by 3. This is a diagnostic pseudo-component basis for the density fit, not a general DES molecular definition.
- Fit success: `False`.
- Fitted nonassociating `m`: `11.319494`.
- Fitted nonassociating `sigma`: `3.09386` A.
- Fitted nonassociating `epsilon/k`: `311.1918` K.
- Density metric: `0.010303518856323288`.

## Equilibrium Smoke

- Electrolyte stability status: `success`.
- Stable flag: `False`.
- Minimum TPD: `-0.5346836582364658`.
- Electrolyte LLE status: `not_accepted`.
- Electrolyte LLE accepted: `False`.
- Electrolyte LLE split detected: `False`.
- Electrolyte LLE residual norm: `2.359750191692001`.
- Electrolyte LLE acceptance gate: `predictive_budget_exhausted`.
- Electrolyte LLE best-effort phases returned: `True`.

## Source-Gated Extraction Comparison

- The 2025 SI equilibrium mole fractions are tracked at `analyses/rezaee_2026_pcsaft_epcsaft/data/input/rezaee_2025_extraction_equilibrium_mole_fractions.csv`.
- Rezaee 2026 closes extraction by phase-specific reaction equilibrium, not by a conventional same-species LLE fugacity equality.
- Run `scripts/rezaee_reactive_equilibrium_replay.py` for the current chemical-equilibrium replay and source-convention diagnostics.
- Legacy package-local fitting was removed; future calibration must route through the native regression gate.

## Interpretation

The density regression and electrolyte-stability calls exercise the current ePC-SAFT package successfully. The direct LLE call is kept as a diagnostic: if it returns a collapsed/non-predictive candidate for this pseudo-DES system, that is recorded as model-support evidence rather than hidden behind downstream calibration.
