# Ascani 2023 Reactive Phase-Equilibrium Validation

This lane targets the Ascani 2023 esterification chemical-plus-phase equilibrium case through public generic `epcsaft` reactive phase-equilibrium APIs and native Ipopt.

Current status: `blocked_source_data`. The public native Ipopt route gates are retained for the Ascani 2023 hypothetical reactive split and the System 1 esterification route, but Figure 8/9 and reactive tie-line paper-match claims remain blocked because no machine-readable source target rows for feed basis, phase compositions, or tie-line compositions are present in `data/` yet.

The hypothetical fixture uses the Table 3 PC-SAFT parameters, Figure 3 \(K_a=2.25\), and an explicit Eq. (10) source-activity to native liquid `x*phi` standard-state conversion. It is an algorithm fixture, not an experimental paper-match validation.

Run from the repository root:

```powershell
uv run python analyses\paper_validation\native\2023_ascani\scripts\run_all.py
```

The script writes `results/reactive_phase_equilibrium/summary.json` and exits nonzero while source target rows are missing.
