# Yu 2024 Li Mg IL Extraction Validation

This lane is a source-backed scaffold for Yu 2024, "Highly efficient lithium extraction from magnesium-rich brines with ionic liquid based collaborative extractants."

The validation target is Li+/Mg2+ extraction from simulated brine using [HOEMIM][Tf2N] plus TOP. It is a later-stage reactive-electrolyte LLE benchmark because it combines phase equilibrium, cation exchange, Li/Mg selectivity, and source-reported ePC-SAFT modeling.

## Source Summary

- Introduction: motivates high-Mg/Li brine separation and identifies IL plus organophosphorus ligand systems as a route away from FeCl3/TBP volatile-solvent extraction.
- Experimental section: defines simulated brine with Mg/Li molar ratio 40:1, extraction protocol, ICP-OES measurement, pH and phase-ratio studies, washing, stripping, and recycle tests.
- Thermodynamic modeling: uses ePC-SAFT with hard-chain, dispersion, association, and Debye-Huckel terms. The paper uses composition-dependent dielectric behavior and reports pure-component and binary parameters for TOP, IL, water, Li+, Mg2+, and chloride-containing systems.
- Quantum chemical calculations: DFT, ESP, and IGM calculations are used to rationalize Li+ transfer, Mg2+ rejection, and the Li+-TOP-[Tf2N]- complex.
- Phase-equilibrium method: the paper explicitly states that phase equilibria are calculated using the isofugacity criterion and phi-phi fugacity equality, with fugacity coefficients derived from residual Helmholtz energy through chemical potentials.
- Extraction chemistry: the work describes cation exchange and a 1:1:1 Li+-TOP-[Tf2N]- complex. The IL cation enters the aqueous phase while Li+ transfers into the organic phase as a stabilized neutral complex.
- Results: useful benchmark figures include ligand screening, IL cation screening, pH/IL concentration/OA optimization, ePC-SAFT extraction-efficiency prediction, Van't Hoff thermodynamic analysis, QC mechanism figures, washing/stripping, and recycling.

## Package Benchmark Role

This is a strong but demanding benchmark candidate for public `reactive_electrolyte_lle` behavior because it exercises:

- Li+ and Mg2+ competition in a high-Mg/Li brine;
- phase-tagged cation exchange;
- selective Li+ organic transfer without forcing fugacity equality between chemically distinct phase species;
- TOP/Li+ association and IL-anion stabilization;
- negligible Mg2+ extraction as a qualitative and quantitative diagnostic;
- source comparison for extraction efficiency and selectivity.

It should follow Hubach 2024 because the route is chemically richer and the currently supplied markdown references supplementary tables that still need to be available locally before this lane can be executable.

## Current Status

Status is `scaffolded_blocked`. The main-paper markdown is currently observed in a sibling checkout, not in this worktree, and the supplementary tables referenced by the paper were not available in the supplied local markdown set. This lane should not be promoted until source tables are present, normalized, and a public native Ipopt reactive-electrolyte LLE solve is accepted.

