# Ascani 2022 Electrolyte LLE Validation

This lane targets an Ascani 2022-style mixed-solvent electrolyte LLE split through the public `epcsaft` API and the native Ipopt electrolyte LLE route.

Current status: executable. The lane proves an accepted native Ipopt solve with a sensible water/ion-rich phase and a butanol-rich phase. It intentionally does not force an exact match to the Ascani 2022 Case 2 mixed NaCl/KCl reported compositions.

Run from the repository root:

```powershell
uv run python analyses\paper_validation\native\2022_ascani\scripts\run_all.py
```

The script writes `results/electrolyte_lle/summary.json` and exits 0 only when `accepted == true`, `solver_backend == "ipopt"`, material balance and charge balance are within tolerance, and the phase split has the expected physical direction.
