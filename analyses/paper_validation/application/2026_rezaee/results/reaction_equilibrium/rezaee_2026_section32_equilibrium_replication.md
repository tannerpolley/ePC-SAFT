# Rezaee 2026 Section 3.2 Equilibrium Replication

## Scope

This script starts at the paper text immediately after Table 8: after all PC-SAFT parameters are determined, reaction equilibrium is calculated to predict equilibrium composition, lithium extraction percentage, and Li/Na selectivity.

The validation uses the paper equations directly rather than package-owned equilibrium solvers:

- Eq. 17 supplies mole fractions from reaction coordinates.
- Eqs. 14 and 15 update RLi and RNa mole fractions from activity coefficients.
- Eqs. 18 and 19 compute Li extraction and selectivity.
- Eq. 20 computes AARD over the 26 benchmark rows.

The package is used only to compute activity coefficients. The first source-aligned package option is Held-2014-S2-like ePC-SAFT with no Born term and constant dielectric behavior for the aqueous phase, plus Table 9 organic binary interactions.

## Result

- Rows evaluated: `26`.
- Paper reference AARD after Table 9 kij: `7.89%` for lithium extraction and `8.63%` for selectivity.
- Direct Held-2014/Table-9/pH-stoichiometric Li extraction AARD: `99.99999999872658`.
- Direct Held-2014/Table-9/pH-stoichiometric selectivity AARD: `56.35551722412048`.
- Best diagnostic case: `held_2014_s2_no_born_table9_kij_inverseK_diagnostic`.
- Best diagnostic Li extraction AARD: `99.68882523044951`.
- Best diagnostic selectivity AARD: `185341.88896748304`.
- Best diagnostic median calculated Li extraction: `0.06900526897921457`.
- Best diagnostic median calculated selectivity: `875.5396522401461`.

## Interpretation

The script follows Rezaee Section 3.2 after Table 8: Eq. 17 updates reaction-coordinate mole fractions, ePC-SAFT/PC-SAFT supplies activity coefficients for Eqs. 14 and 15, and Eqs. 18-20 evaluate Li extraction, Li/Na selectivity, and AARD over the 26 benchmark rows. The direct Held-2014-like no-Born pH-stoichiometric basis is the first source-aligned case.

If the direct pH-stoichiometric case does not approach the paper AARD, the next source task is not a package phase-equilibrium solve. It is identifying the paper's exact initial mole basis for OH-/NH4OH/base inventory behind Eq. 17, because Section 3.2 says those initial moles are inputs but the main-text DOE table only exposes pH.

## Generated Files

- `data\processed\rezaee_2026_section32_equilibrium_replication_rows.csv`
- `results\reaction_equilibrium\rezaee_2026_section32_equilibrium_replication_summary.json`
