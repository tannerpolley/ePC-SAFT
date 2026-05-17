# Native Associating-Binary Ceres Benchmark

Date: 2026-05-16

This report preserves the useful benchmark evidence from the retired
`native-eos-association-ceres-benchmark` GoalBuddy board without keeping the
stale board state as an active package artifact.

## Benchmark Fixture

- System: ethanol + water VLE at 100 kPa.
- Dataset: `2012_Held`.
- Fixture metadata: `tests/fixtures/literature/binary_kij/ethanol_water_jced2021_100kpa.json`.
- VLE rows: `data/reference/regression/binary/ethanol_water_jced2021_vle_100kpa.csv`.
- Source notes: `data/reference/regression/binary/ethanol_water_jced2021_source_notes.md`.
- Source DOI recorded in rows: `10.1021/acs.jced.0c00686`.
- Smoke rows used by the package test: 5.
- Initial `k_ij`: `-0.049`.
- Reference paper `k_ij`: `-0.0269`.
- Bounds: `[-0.15, 0.10]`.

The repo still does not carry a normalized MEA + water binary VLE fixture, so
this ethanol + water case is the current package-owned associating-binary proof.

## Public API Call

```python
import epcsaft
from tests.helpers.binary_regression_cases import (
    ETHANOL_WATER_HELD_2012_KIJ,
    ethanol_water_jced2021_vle_records,
)

result = epcsaft.fit_binary_pair(
    ethanol_water_jced2021_vle_records(smoke_only=True),
    ("Ethanol", "H2O"),
    dataset="2012_Held",
    initial_guess={"k_ij": ETHANOL_WATER_HELD_2012_KIJ},
    bounds={"k_ij": (-0.15, 0.10)},
)
```

## Current Result

Observed with the current branch on 2026-05-16:

- `success`: `true`.
- `backend`: `ceres`.
- `optimizer_backend`: `ceres`.
- `derivative_backend`: `cppad_implicit`.
- Fitted `k_ij`: `-0.028543168162295768`.
- Absolute error from reference `k_ij`: `0.0016431681622957675`.
- Initial objective: `0.015197798818532296`.
- Final objective: `0.00023620672364734304`.
- `binary_vle_fugacity_balance` metric: `0.02173507412673548`.
- Source record count: `5`.

Acceptance bands enforced by the tracked test:

- `binary_vle_fugacity_balance < 0.04`.
- `abs(fitted k_ij - reference k_ij) < 0.01`.
- Final objective is lower than initial objective.
- Optimizer backend is native Ceres.
- Derivative backend is `cppad_implicit`.

## Verification

Current focused proof:

```powershell
uv run python run_pytest.py tests/regression/literature/test_ethanol_water_binary_vle_regression.py -q
```

The normal package quick gate also includes representative regression coverage:

```powershell
uv run python scripts/dev/validate_project.py quick
```
