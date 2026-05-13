## Issue #66 response plan

Source comment:

`https://github.com/tannerpolley/ePC-SAFT/issues/66#issuecomment-4436520857`

The response confirms the pause was appropriate and that the current pure-neutral Ceres slice is useful but insufficient to close issue #66.

## Decision Applied

The issue response gives two paths:

1. Preferred: open the current pure-neutral Ceres work as a partial foundation PR without `Closes #66`, then start a focused binary derivative-subsystem branch.
2. Acceptable: continue on the current branch while keeping the pure-neutral work logically separated from the binary `k_ij` derivative-subsystem work.

Because the user explicitly resumed `/goal` and asked to proceed with accomplishing the full goal, the board will use option 2:

```text
Continue on codex/issue-66-ceres-regression-backend-2.
Keep pure-neutral Ceres foundation logically isolated.
Implement binary k_ij derivative subsystem next.
Do not open a Closes #66 PR until pure + binary Ceres support is validated.
```

## Continuation Constraints

- Ceres owns optimizer loop.
- Python validates, serializes, and shapes results only.
- No finite differences.
- No missing Jacobian columns.
- No Python-owned production objective loop.
- No optimistic capabilities.
- Unsupported row/parameter families return `backend_unavailable`.
- Do not differentiate through density or association solver iterations; use implicit sensitivities where needed.

## Next Worker Slice

Implement a bounded native neutral binary `k_ij` derivative slice:

1. Treat `k_ij` as an independent native derivative variable for neutral dispersion mixing.
2. Produce explicit algebraic partials for each phase:
   - pressure value
   - `dP/drho`
   - `dP/dk_ij` at fixed rho
   - `drho/dk_ij = -(dP/dk_ij)/(dP/drho)`
   - `lnphi_i`
   - `(d lnphi_i/dk_ij)_rho`
   - `d lnphi_i/drho`
   - total `d lnphi_i/dk_ij` along the solved density root
3. Wire binary VLE residual/Jacobian evaluation to Ceres.
4. Keep association, ionic, liquid-electrolyte, and fully coupled reactive rows gated unless they have real derivative support.
5. Add tests proving the Ceres binary path does not use finite differences, does not use Python production objective evaluation, moves at least one `k_ij`, populates diagnostics, and gates unsupported families honestly.

## Verification Target

```powershell
uv run python scripts/build_epcsaft.py --clean --enable-ceres --enable-cppad
uv run python run_pytest.py tests/native/test_ceres_binary_regression.py tests/native/test_ceres_pure_regression.py -q
uv run python run_pytest.py tests/regression/test_ethanol_water_binary_vle_regression.py tests/api/test_regression_api.py::test_fit_binary_pair_vle_kij_default_and_rejects_temperature_models -q
uv run python scripts/validate_project.py quick
```

If derivative coverage gate tests are added:

```powershell
uv run python run_pytest.py tests/native/test_ceres_derivative_coverage_gates.py -q
```
