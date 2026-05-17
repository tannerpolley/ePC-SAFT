# Ascani 2023 Reactive Phase-Equilibrium Validation

This lane targets the Ascani 2023 esterification chemical-plus-phase equilibrium case through public generic `epcsaft` reactive phase-equilibrium APIs and native Ipopt.

Current status: blocked. The local paper markdown contains the paper identity, pure-component parameters, binary interaction parameters, and reported equilibrium constants, but no machine-readable source target rows for feed basis, phase compositions, or tie-line compositions are present in `data/` yet. The lane must not substitute random toy fixtures or route-shape tests for literature validation.

Run from the repository root:

```powershell
uv run python analyses\paper_validation\native\2023_ascani\scripts\run_all.py
```

The script writes `results/reactive_phase_equilibrium/summary.json` and exits nonzero while source target rows are missing.
