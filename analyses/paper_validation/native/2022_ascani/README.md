# Ascani 2022 Electrolyte LLE Validation

This lane targets the Ascani 2022 mixed-solvent mixed-electrolyte LLE case through the public `epcsaft` API and the native Ipopt electrolyte LLE route.

Current status: blocked. The local Ipopt profile imports and reports Ipopt available, but the public electrolyte LLE solve for the source-backed Case 2 feed is rejected by Ipopt before an accepted phase split is returned.

Run from the repository root:

```powershell
uv run python analyses\paper_validation\native\2022_ascani\scripts\run_all.py
```

The script writes `results/electrolyte_lle/summary.json` before exiting nonzero when the accepted-solve gate is not met. The lane must stay blocked in the benchmark registry until that command exits 0 with `accepted == true`, `solver_backend == "ipopt"`, material balance within tolerance, charge balance within tolerance, and source-derived composition checks.
