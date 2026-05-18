# Hubach 2024 Li TCB Extraction Validation

This lane is a source-backed scaffold for Hubach et al. 2024, "Li+ Extraction from Aqueous Medium Using Tetracyanoborate Ionic Liquids - Experiments and ePC-SAFT Modeling."

The validation target is lithium extraction from aqueous LiCl using tetracyanoborate ionic liquids plus TBP, with emphasis on public package electrolyte/reactive LLE outputs rather than application-specific package APIs.

## Source Summary

- Introduction: motivates direct lithium extraction and positions IL/TBP systems as a lower-acidity alternative to traditional TBP/FeCl3 extraction systems.
- Experimental section: defines the extraction procedure, Li+ feed ranges, organic/aqueous phase ratio, TBP/IL ratio, and analytical measurements used to compute extraction efficiency.
- Thermodynamic modeling: uses ePC-SAFT advanced with hard-chain, dispersion, association, Debye-Huckel, and Born terms. The paper explicitly notes composition-dependent permittivity and temperature-corrected ion diameters.
- Phase-equilibrium method: the paper states that phase equilibria are calculated with an isofugacity criterion and a phi-phi fugacity relation, with fugacity coefficients derived from residual Helmholtz energy through chemical potentials.
- Mechanism: the extraction is treated as IL-assisted Li+ transfer with TBP/Li+ coordination represented through an additional Li+ association site. Cross-association is allowed for Li+/TBP and disabled for Li+/water and Li+/IL.
- Results: useful benchmark figures include IL cation solubility, TBP solubility, IL-cation screening, Li feed-concentration dependence, TBP-content dependence, organic/aqueous ratio dependence, and ePC-SAFT extraction-efficiency comparison.

## Package Benchmark Role

This is a strong benchmark candidate for electrolyte LLE and reactive-electrolyte transfer because it exercises:

- distributed ions across aqueous and organic liquid phases;
- liquid-root phase equilibrium at fixed temperature and pressure;
- TBP/Li+ association or reaction-like coordination;
- Born plus Debye-Huckel electrolyte terms;
- composition-dependent dielectric behavior;
- source-reported experimental and modeled rows in the supporting information.

The first executable gate should reproduce current package model outputs against Figure 7 / Table S11 because that is the clearest full extraction-efficiency comparison with paper ePC-SAFT values.

## Current Status

Status is `scaffolded_blocked`. The paper markdown is currently observed in a sibling checkout, not in this worktree. This lane should not be promoted until source assets are copied or normalized into this repo and a public package route produces accepted native Ipopt results.

