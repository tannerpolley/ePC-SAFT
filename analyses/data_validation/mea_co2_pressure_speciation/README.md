# MEA-CO2 pressure and speciation data validation

This lane validates the package-owned MEA-CO2-H2O pressure/speciation path using copied Phase 2 MEA artifacts from `MEA-Thermodynamics`.

The lane is intentionally under `analyses/data_validation` rather than `analyses/paper_validation`: it validates fitted parameter artifacts and pressure/speciation behavior, not a single-paper reproduction workflow.

Run from the repository root:

```powershell
uv run python analyses\data_validation\mea_co2_pressure_speciation\scripts\run_all.py
```

The command writes retained inputs and diagnostics under `data/processed` and `results/pressure_speciation`. It exits `0` only after the public reactive speciation API produces an accepted native Ipopt solve. Until then, `summary.json` records the exact package blocker and the command exits nonzero.
