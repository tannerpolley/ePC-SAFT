# Ascani 2022 Electrolyte LLE Validation

This lane validates Ascani 2022 Case Study 2 through the public `epcsaft` API:

- species: `H2O`, `Butanol`, `Na+`, `K+`, `Cl-`
- feed mass fractions: water `0.8094`, 1-butanol `0.1728`, NaCl `0.0054`, KCl `0.0124`
- conditions: `298.15 K`, `1 bar`

Run from the repository root:

```powershell
uv run python analyses\paper_validation\native\2022_ascani\scripts\run_all.py
```

The script exits 0 only when the public native Ipopt liquid-root electrolyte LLE route is accepted and the hard TPDF stability certificate passes. Exact paper matching is recorded as comparison data, not claimed by default.
