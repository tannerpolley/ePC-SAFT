# ePC-SAFT for Lithium Extraction from Produced Water (14-Slide Draft Text)

Use the bullet points under each slide as the *actual on-slide text*.  
Under **Insert (optional)** are the repo assets you can drop into PowerPoint.

---

## Slide 1 — ePC-SAFT for Lithium Extraction from Produced Water
- Predictive electrolyte thermodynamics to support lithium solvent-extraction design
- Status: electrolyte model is calibrated and validated; produced-water + LLE coupling is close
- Today: review model evolution, what was implemented, and what the fits show
- Decision needed: align on the recommended baseline (2020) and the next validation milestones

---

## Slide 2 — Why Produced Water Needs a Predictive Model
- Produced water chemistry is variable: salinity, competing ions, temperature, and organics all change performance
- Lithium is a trace species; selectivity depends on the full electrolyte matrix, not just LiCl-in-water
- A physics-based EOS reduces trial-and-error and makes extrapolation defensible
- Accurate activities and phase splits are prerequisites for any extraction unit model

---

## Slide 3 — Model Lineage: 2001 → 2005 → 2008 → 2014 → 2020
- **2001 (PC-SAFT):** neutral-fluid backbone (hard-chain + dispersion + association)
- **2005 (ePC-SAFT):** adds Debye–Hückel electrostatics for ion–ion interactions in a dielectric continuum
- **2008:** expands aqueous salt coverage and MIAC-focused parameter fitting/validation
- **2014 (revised):** adds explicit cation–anion dispersion to improve concentrated brines and multisolvent behavior
- **2020 (advanced):** adds dielectric-dependent DH + altered Born term to enable water ↔ organic transferability

Insert (optional):
- Paper summaries: `docs/papers/md/Gross, Sadowski - 2001 - PC-SAFT An equation of state based on a perturbation theory for chain molec.md`
- Paper summaries: `docs/papers/md/Cameretti, Sadowski, Mollerup - 2005 - Modeling of Aqueous Electrolyte Solutions with Perturbed-Chai.md`
- Paper summaries: `docs/papers/md/Held, Cameretti, Sadowski - 2008 - Modeling aqueous electrolyte solutions. Part 1. Fully dissociated.md`
- Paper summaries: `docs/papers/md/Held et al. - 2014 - ePC-SAFT Revised.md`
- Paper summaries: `docs/papers/md/Bülow, Ascani, Held - 2020 - ePC-SAFT advanced - Part I Physical meaning of including a concentratio.md`

---

## Slide 4 — What Actually Changed Between Versions (In One Table)
- 2001: transferable pure-component parameters for neutral mixtures (few adjustable knobs)
- 2005/2008: ions added as charged spheres; Debye–Hückel captures long-range screening effects
- 2014: ion–ion dispersion added (anion–cation) to fix high-salt deviations without breaking aqueous performance
- 2020: concentration-dependent permittivity drives both DH and Born, improving non-aqueous consistency

---

## Slide 5 — ePC-SAFT Structure (Where the Electrolyte Physics Lives)
- Neutral-fluid contributions come from PC-SAFT: hard-chain, dispersion, and association
- Electrolyte contributions are added as explicit long-range terms (DH) and optional solvation (Born)
- All properties needed for phase equilibrium come from derivatives of the residual Helmholtz energy

Equation (paste into PowerPoint equation editor):
```tex
a^{res} = a^{hc} + a^{disp} + a^{assoc} + a^{DH} + a^{Born}
```

---

## Slide 6 — Debye–Hückel vs Born (What Each Term Represents)
- **Debye–Hückel (DH):** long-range ion–ion electrostatics and screening; primary driver of MIAC trends vs salt level
- **Born:** ion–solvent solvation free energy; becomes critical when moving from water to lower-permittivity solvents
- Separating DH and Born makes debugging physical: ion–ion mismatch vs ion–solvent mismatch
- Born can be small in mostly-aqueous cases, but it is not optional for mixed-solvent transfer predictions

---

## Slide 7 — Relative Permittivity (Dielectric Constant) Is the Transferability Lever
- Permittivity sets the “strength” of electrostatics in both DH and Born terms
- Constant permittivity can work locally, but it breaks when solvent composition or salt loading shifts ε meaningfully
- The 2020 approach treats ε as composition-dependent and includes the needed derivatives for chemical potentials
- This is the key ingredient for credible water → alcohol/organic electrolyte predictions

Insert (optional):
- Dielectric fits: `data/dielc/plot_fits/test_dielc_fit_salts_in_water_fit.png`
- Dielectric derivative check: `data/dielc/plot_fits/test_dielc_diff_salts_in_water_fit.png`

---

## Slide 8 — What We Implemented in This Package (Compared to the Papers)
- Side-by-side model presets (2005 / 2008 / 2014 / 2020) to reproduce literature assumptions consistently
- A dielectric “engine” with selectable mixing rules and analytic derivatives (with safe numerical fallback)
- Altered Born contribution options to match 2020-style composition-dependent ε behavior
- Guardrails: SSM+DS-style Born solvation shell + dielectric saturation is not part of the current baseline workflow

Insert (optional):
- Model presets: `data/epcsaft_properties.py`

---

## Slide 9 — Current Restrictions and How We’re Avoiding Over-Claims
- Produced-water extraction is “near-ready”: electrolyte core is validated; full brine + LLE validation is still in progress
- Current recommendation for extraction work: **2020 preset** (DH + altered Born with composition-dependent ε)
- SSM+DS (solvation shell + dielectric saturation) is intentionally not used for current produced-water claims
- We are not reporting final extraction efficiency until two-phase LLE and partitioning are validated end-to-end

---

## Slide 10 — Fit Coverage: Water (Alkali Halides)
- Water MIAC fits cover: NaCl/NaBr/NaI, KCl/KBr/KI, LiCl/LiBr/LiI
- Osmotic-coefficient reproduction checks keep the model anchored to published aqueous benchmarks
- These fits establish the ion parameter base before moving into mixed-solvent and extraction scenarios

Insert (optional):
- Water MIAC fits: `data/MIAC_m/water/plot_fits/maic_m_water_*.png`
- Osmotic validation: `data/osmotic/water/plot_fits/validation_2014_repro_NaCl_KBr_fit.png`

---

## Slide 11 — Fit Coverage: Methanol and Ethanol (Transferability Stress Test)
- Alcohol MIAC fits demonstrate behavior in lower-permittivity solvents relevant to extraction workflows
- Methanol coverage: NaCl/NaBr/NaI, KBr/KI, LiCl/LiBr
- Ethanol coverage: LiCl/LiBr, NaBr
- The non-aqueous performance is where the 2020 dielectric-dependent Born treatment shows its value

Insert (optional):
- Methanol MIAC fits: `data/MIAC_m/methanol/plot_fits/maic_m_methanol_*.png`
- Ethanol MIAC fits: `data/MIAC_m/ethanol/plot_fits/maic_m_ethanol_*.png`

---

## Slide 12 — Why the 2020 Model Is the Best Baseline for This Program
- Best balance today: strong aqueous grounding plus credible transfer toward non-aqueous/mixed-solvent systems
- Keeps the parameter philosophy clean: mostly pure-component/ion parameters with a consistent kij strategy
- Avoids adding extra empirical knobs unless the data demands it (better for scale-up and extrapolation)
- Sets the right foundation for LLE coupling because electrostatics respond correctly to solvent composition

Insert (optional):
- DH/Born model comparisons: `data/DH_born_models/water/plot_fits/test_DH_born_models_water_*_all_fit.png`

---

## Slide 13 — LLE for Solvent Extraction (Slide Ready for When the Solver Lands)
- Water + organic phases can split; lithium distribution is controlled by activities in each liquid phase
- ePC-SAFT gives phase split + ionic activities from one thermodynamic framework (not separate correlations)
- Distribution ratio predictions can be converted into stagewise extraction efficiency in a unit model
- This capability is being finalized in parallel; this slide is positioned for the “ready to demo” milestone

Equation (optional, keep simple and manager-friendly):
```tex
D_{Li} = \frac{m_{Li}^{org}}{m_{Li}^{aq}}
```

Insert (optional):
- Multiphase equilibrium reference: `docs/papers/md/Ascani, Sadowski, Held - 2022 - Calculation of Multiphase Equilibria Containing Mixed Solvents and M.md`

---

## Slide 14 — Readiness, Risks, and Next Milestones
- Done: electrolyte model calibration + dielectric handling + fit evidence across water and alcohol solvents
- In progress: produced-water composition mapping and two-phase LLE validation for the chosen organic system
- Next deliverables: phase diagrams, distribution ratios vs salinity/temperature, and sensitivity to competing ions
- Ask: approve the 2020 baseline and the final validation path to connect thermodynamics to an extraction unit model

