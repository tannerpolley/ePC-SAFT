# Rezaee 2026 Reactive Equilibrium Replay

## Source Basis

- 2025 source: `papers/pdf/Rezaee et al. - 2025 - Application of Response Surface Methodology for Selective Extraction of Lithium Using a Hydrophobic DES.pdf`.
- 2025 supporting information: `papers/pdf/Rezaee et al. - 2025 - Supporting information - Application of Response Surface Methodology for Selective Extraction of Lithium.pdf`.
- Main source: `papers/pdf/Rezaee et al. - 2026 - Thermodynamic modeling of lithium extraction from synthetic brine using deep eutectic solvents A PC.pdf`.
- 2026 supporting information: `papers/pdf/Rezaee et al. - 2026 - Supplementary material - Thermodynamic modeling of lithium extraction from synthetic brine using deep eutectic solvents.pdf`.
- Searchable source text: `papers/md/Rezaee et al. - 2026 - Thermodynamic modeling of lithium extraction from synthetic brine using deep eutectic solvents A PC.md`.
- The replay follows the paper's phase-specific reaction-equilibrium formulation rather than a conventional same-species LLE fugacity equality.
- The package check uses `ReactionDefinition.phase_stoichiometry` and the native phase-tagged cross-phase residual block.
- Aqueous phase: H2O, Li+, Na+, Cl-, H+, OH-, NH4+ with ePC-SAFT component activity coefficients.
- Organic phase: DES, TOPO, RLi, RNa with PC-SAFT activity coefficients calculated as mixture fugacity coefficient over pure-component fugacity coefficient.

## Result

- Rows replayed: `26`.
- Status: `source_mismatch`.
- Paper lnK Li/Na: `-19.52985865272526`, `-27.493348310803476`.
- Median lnQ at experimental rows Li/Na: `12.986988440452834`, `10.479242476903497`.
- Median lnQ-lnK Li/Na: `32.516847093178086`, `37.97259078770698`.
- Package phase-tagged cross-phase rows evaluated: `26`.
- Package native median cross-phase residual Li/Na: `408.03479797400666`, `396.5589577644374`.
- Median absolute RLi/RNa error from paper K replay: `0.004442850658281856`, `0.020811015739249972`.
- Mean relative RLi/RNa error from paper K replay: `0.9999999980654698`, `0.9999999980656575`.

## Interpretation

The package can evaluate the phase-tagged cross-phase reaction residual required by Rezaee's formulation. However, using the paper-reported Table 2 equilibrium constants together with the paper/SI composition rows and the paper-reported organic parameters does not reproduce the reported RLi/RNa complex mole fractions under this activity-reference convention.

This is not the same as the old four-species fixed-composition LLE smoke. This replay includes the chemical-equilibrium equations that control lithium and sodium extraction. The current blocker is the source/model convention gap exposed by lnQ-lnK, not an omitted call to `electrolyte_lle`.

Next implementation step: resolve the convention gap by checking the 2026 supporting information/group-contribution worksheet and the 2025 phase-amount basis, then either correct the stored constants/reference-state convention or add a calibrated Rezaee parameter-refit lane that uses these EOS activity calls directly.
