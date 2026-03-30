# 2026 Khudaida analysis provenance

- Tables 3 and 4 in the local Khudaida markdown/PDF are treated as the canonical experimental tie-line source for Figures 2-7 and for the salted points in Figures 8-9.
- `2026_Khudaida` uses the paper's Table 5 pure-component parameters, Table 7 binary interaction parameters, and an exact copy of `2025_Figiel/user_options.json` as requested.
- Figure 1 salt-free data and the no-salt points in Figures 8-9 were reconstructed from the local paper figures because the Zotero baseline source remained inaccessible in this session.
- The no-salt baseline is therefore marked as `digitized_local_paper` in the emitted CSV files.
- Tables 9 and 10 include package-generated ePC-SAFT AAD values and paper-copied eNRTL/ePC-SAFT reference values for comparison.
- The current package solver can return shallower liquid-liquid splits than the paper for this system; the analysis therefore uses an opt-in external-seed path into `pcsaft_multiphase_lle` while leaving the default solver path unchanged for other workflows.
