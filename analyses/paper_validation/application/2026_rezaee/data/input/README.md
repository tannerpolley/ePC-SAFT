# Rezaee 2026 Input Data

This folder contains source-backed Rezaee data that can be used without rerunning a literature search.

- Source papers are stored under `papers/pdf/` and `papers/md/`, including the 2025 experimental article, 2025 supporting information, 2026 modeling article, and 2026 supporting information.
- `rezaee_2026_reaction_constants.csv`: reaction A/B standard Gibbs-energy changes and equilibrium constants from the 2026 supporting information.
- `rezaee_2026_organic_pcsaft_parameters.csv`: paper-reported DES, TOPO, RLi, and RNa PC-SAFT parameters.
- `rezaee_2026_organic_binary_interactions.csv`: paper-reported organic-phase binary interaction parameters.
- `rezaee_2026_si_density_tables.csv`: density tables transcribed from the 2026 supporting information. Table S5 is DES; Table S6 is 10 wt% TOPO in DES.
- `rezaee_2025_headline_extraction_points.csv`: source-backed headline extraction percentages from the ACS abstract for the 2025 experimental paper.
- `rezaee_2025_doe_extraction_responses.csv`: Table 5 designed-experiment Li/Na extraction response rows from the 2025 experimental paper.
- `rezaee_2025_optimum_neighborhood.csv`: Table 7 optimum-neighborhood validation rows.
- `rezaee_2025_screening_extraction.csv`: Table 3 DES-ratio screening and Table 4 coextractant screening rows.
- `rezaee_2025_real_brine_extraction.csv`: Table 8 real/synthetic brine composition and cation extraction efficiency rows.
- `rezaee_2025_extraction_equilibrium_mole_fractions.csv`: 26 designed-experiment aqueous and organic phase mole-fraction rows transcribed from the 2025 ACS supporting information Tables S1 and S2.
- `rezaee_2026_epcsaft_species_dataset.json`: analysis-local manifest of the species, source parameter tables, reaction constants, organic binary interactions, and ePC-SAFT option context used by the Rezaee validation scripts.

The 2025 experimental paper is `10.1021/acs.iecr.4c03496`; its supporting-information file is listed by ACS as `ie4c03496_si_001.pdf`. The 2026 PC-SAFT/ePC-SAFT paper is `10.1016/j.fluid.2026.114737`; its supporting-information file is `1-s2.0-S0378381226000786-mmc1`. Rezaee 2026 states that the phase-equilibrium rows come from the 2025 work, so both source years are required for reproduction. The current benchmark set is the 26 source-backed SI equilibrium rows; treat the 2026 main-text statement of 36 equilibrium points as a clerical mismatch unless a source-backed workbook is later supplied. Do not fill source input rows from model output or inferred values. Add only source-backed rows from the 2025/2026 papers, their supporting information, or from a verified digitization workflow.
