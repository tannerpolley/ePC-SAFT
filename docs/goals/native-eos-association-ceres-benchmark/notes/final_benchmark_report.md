# Final Benchmark Report

Goal: `native-eos-association-ceres-benchmark`

Issue: GitHub issue #114, "Complete native EOS/AD, association implicit sensitivities, and associating-binary Ceres regression"

## Benchmark Fixture

The benchmark uses the repo-contained ethanol + water 100 kPa VLE fixture selected during Stage 0:

- Fixture metadata: `tests/fixtures/literature/binary_kij/ethanol_water_jced2021_100kpa.json`
- VLE data: `data/reference/regression/binary/ethanol_water_jced2021_vle_100kpa.csv`
- Source notes: `data/reference/regression/binary/ethanol_water_jced2021_source_notes.md`
- Source DOI recorded in rows: `10.1021/acs.jced.0c00686`
- Smoke rows used by the benchmark: 5
- Parameter set: `2012_Held`
- Initial `k_ij`: `-0.049`
- Reference `k_ij`: `-0.0269`
- Bounds: `[-0.15, 0.10]`

This is the Stage 0-approved associating alcohol-water substitute. No normalized repo-contained MEA + water binary VLE fixture was found during intake.

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
    multistart=0,
)
```

## Result

- `success`: `true`
- `backend`: `ceres`
- `optimizer_backend`: `ceres`
- `derivative_backend`: `cppad_implicit`
- `jacobian_backend`: `cppad_implicit`
- `python_objective_used`: `false`
- `hessian_backend`: `not_implemented`
- Fitted `k_ij`: `-0.02850385213142522`
- Absolute error from reference `k_ij`: `0.0016038521314252208`
- Parameter movement: `{"k_ij": 0.02049614786857478}`
- Active bounds: `[]`
- Initial objective: `0.015197798818532296`
- Final objective: `0.00023626180683943726`
- `binary_vle_fugacity_balance` metric: `0.02173760827871536`
- Row diagnostics: `[{"row_family": "binary_vle_fugacity_balance", "metric": 0.02173760827871536}]`
- Residual block norms: `{"binary_vle_fugacity_balance": 0.02173760827871536}`
- Target family summaries: `{"binary_vle_fugacity_balance": {"record_count": 5, "residual_block_norm": 0.02173760827871536}}`
- Source summaries: `{"records": {"dataset": "2012_Held", "record_count": 5, "row_families": ["binary_vle_fugacity_balance"]}}`

The fit passes the benchmark acceptance bands in the fixture:

- `binary_vle_fugacity_balance < 0.04`
- `abs(fitted k_ij - reference k_ij) < 0.01`

## Production-Path Proof

The benchmark exercises the public Python regression API, but the production optimizer loop is native Ceres. Python prepares records, validates inputs, calls native code, and formats the returned `FitResult`.

The native Ceres binary VLE residual uses the native fugacity-balance residual and consumes the native `k_ij` residual/Jacobian path. For the symmetric fitted parameter, the Ceres callback sums the forward and reverse matrix-entry derivatives. For active association, the constant-pressure derivative combines CppAD explicit EOS terms with association site-fraction implicit density sensitivities rather than differentiating through the site-solve loop.

The result metadata proves the required chain for this benchmark:

- optimizer backend is native Ceres;
- derivative and Jacobian backend are `cppad_implicit`;
- association site-fraction implicit density sensitivity participates in the constant-pressure derivative path;
- the initial objective is greater than the final objective;
- `k_ij` moves from the initial value;
- no active bound is responsible for the fitted value;
- row diagnostics, source summaries, target family summaries, and residual block norms are populated.

## Verification Commands

These commands passed during Stage 4 or Stage 5 work:

```powershell
uv run python scripts/dev/build_epcsaft.py --clean --enable-ceres --enable-cppad --parallel 4
uv run python run_pytest.py tests/native/ceres -q
uv run python run_pytest.py tests/api/regression -q
uv run python run_pytest.py tests/regression/literature/test_literature_binary_kij_regression.py tests/regression/literature/test_ethanol_water_binary_vle_regression.py -q
git diff --check
```

Observed results:

- Clean Ceres+CppAD build: native import OK.
- Native Ceres tests: `5 passed`.
- API regression tests: `35 passed`.
- Focused literature binary tests: `4 passed`.
- `git diff --check`: no whitespace errors; CRLF-normalization warnings only.

## Residual Risk

This report proves the Stage 5 ethanol + water associating-binary `k_ij` benchmark selected during intake. It does not claim a MEA + water benchmark because no normalized repo-contained MEA + water binary fixture was available.

The issue is not complete until the T011 validation ladder and T012 final audit pass. In particular, the broader regression suite, quick validation, and docs validation still need to run before any PR/merge/main-sync handling.
