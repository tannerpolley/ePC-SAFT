# Theory Overview

## Core Model

ePC-SAFT is a Helmholtz-residual equation of state. Most thermodynamic properties are derived from the reduced residual Helmholtz energy and its density, temperature, and composition derivatives.

The implementation splits the residual Helmholtz contribution into named terms:

- `hc`: hard-chain reference contribution.
- `disp`: dispersion contribution.
- `assoc`: association contribution.
- `ion`: electrolyte / Debye-Huckel contribution.
- `born`: Born solvation contribution.

The Python wrappers expose these terms through scalar and vector contribution payloads backed by `epcsaft._core`, with `hc`, `disp`, `assoc`, `ion`, `born`, and `total` keys.

## Contribution Families

- Hard-chain: segment packing, segment diameters, packing fractions, hard-sphere contact values, and chain-reference terms.
- Dispersion: segment-energy and size mixing rules, dispersion polynomial state, and density/composition derivatives.
- Association: association sites, cross-association energy/volume mixing, site fractions, and association contribution derivatives.
- Ion / Debye-Huckel: ionic diameter choices, charge-weighted terms, dielectric state, screening parameter behavior, and electrolyte long-range contribution.
- Born: ion solvation contribution driven by charge/radius and relative-permittivity behavior, including SSM/DS paths when enabled.
- Relative permittivity and electrolyte helpers: solvent/ion index sets, molecular-weight averages, solvent composition helpers, and dielectric model support.

## Derived Properties

- `residual_helmholtz(...)`: primary residual Helmholtz value, with optional contribution terms.
- `compressibility_factor(...)`: depends on density derivative of residual Helmholtz and returns `Z`, with optional contribution terms.
- Density closure: pressure-based states solve `T, P, x -> rho` during state construction.
- `residual_chemical_potential(...)`: composition derivative based residual chemical potential terms.
- `fugacity_coefficient(...)`: defaults to natural-log fugacity coefficients; `natural_log=False` returns coefficient form.
- `activity_coefficient(...)`: returns component activity coefficients or mean ionic values depending on `mean_ionic_form` and `basis`.

## Practical Review Rule

When a numerical behavior changes, check the contribution accounting first: confirm the affected contribution term, then trace how it enters residual Helmholtz, `Z`, residual chemical potential, fugacity, activity, or density closure.
