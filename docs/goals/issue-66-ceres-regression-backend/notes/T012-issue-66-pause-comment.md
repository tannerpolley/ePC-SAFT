## Pause handoff after Ceres pure slice and derivative-subsystem scope decision

Pausing work here at the user's request before starting the larger native derivative-subsystem implementation.

### Current branch and scope

- Active branch: `codex/issue-66-ceres-regression-backend-2`
- Issue scope remains full issue #66.
- PR #56 was not used as a base.
- Ceres must own the optimizer loop.
- Python must only validate/serialize regression inputs and outputs.
- Hard constraint for continuation: no finite differences whatsoever in the Ceres production derivative path.

### What is already implemented locally

The current branch has a validated first slice:

- `scripts/build_epcsaft.py` accepts `--enable-ceres` / `--use-system-ceres`.
- The local Ceres build path works with the MinGW/Ninja toolchain and bundled Eigen/MINIGLOG settings.
- Native pure-neutral `m`, `sigma`, and `epsilon/k` density plus pure VLE regression now has a Ceres-owned optimizer path.
- That pure-neutral Ceres path uses the existing native residual/Jacobian evaluator, not a finite-difference fallback.
- Python-side regression code is limited to validation/serialization and result shaping for that path.
- Binary and liquid-electrolyte Ceres paths are currently explicit `backend_unavailable` gates, not production support.

Validation already run before the pause:

```powershell
uv run python scripts/build_epcsaft.py --clean --enable-ceres --enable-cppad
uv run python run_pytest.py tests/native/test_ceres_pure_regression.py tests/native/test_ceres_binary_regression.py tests/native/test_ceres_liquid_electrolyte_regression.py -q
uv run python scripts/validate_project.py quick
uv run python run_pytest.py tests/api/test_regression_api.py::test_fit_binary_pair_vle_kij_default_and_rejects_temperature_models tests/regression/test_ethanol_water_binary_vle_regression.py::test_native_binary_kij_regression_matches_real_ethanol_water_vle_reference_band -q
```

### Why full issue #66 was not ready to publish

The narrow pure-neutral slice is not enough to close issue #66 because the issue asks for Ceres backends for pure and binary fitting. Production binary `k_ij` support cannot safely be wired by calling the existing generic binary optimizer.

The existing generic binary path evaluates VLE residuals like:

```text
log(x_i) + lnphi_liq_i - log(y_i) - lnphi_vap_i
```

at fixed `T`, `P`, liquid composition, and vapor composition. For Ceres to own this as a real production backend, the residual Jacobian with respect to `k_ij` must be native and real. Reusing a finite-difference Jacobian would violate the issue constraints and the user's explicit continuation constraint.

### Derivative problem found

For a binary VLE `k_ij` parameter, a correct derivative needs more than `d ares / d k_ij` at a fixed density.

The derivative needs:

- `d lnphi_i / d k_ij` for both liquid and vapor phases.
- The density-root implicit sensitivity:

```text
d rho / d k_ij = -(dP/dk_ij) / (dP/drho)
```

- A total derivative of fugacity coefficients along each solved density root:

```text
d lnphi_i/dk_total = (d lnphi_i/dk)_rho + (d lnphi_i/drho) * (d rho/dk)
```

The existing code has pieces of the required thermodynamics, but it does not currently expose the needed parameter-derivative subsystem for `k_ij` through dispersion mixing, pressure, residual chemical potentials, fugacity coefficients, and implicit density roots.

Relevant implementation owners identified:

- `src/epcsaft/native/epcsaft_core_internal.h`
  - `pair_epsilon_cpp(...)` applies `epsilon_ij *= (1.0 - k_ij[idx])`.
  - This is where the binary parameter enters neutral dispersion mixing.
- `src/epcsaft/native/epcsaft_ares.cpp`
  - Has templated scalar and CppAD-oriented residual Helmholtz machinery.
  - Current CppAD composition derivative recording is over composition variables, not over `k_ij` as an independent parameter.
- `src/epcsaft/native/epcsaft_mu.cpp`
  - Builds residual chemical potentials from `ares`, `z`, `dadx_i`, and `sum_x_dadx`.
- `src/epcsaft/native/epcsaft_fugcoef.cpp`
  - Builds `lnphi` from residual chemical potentials and the compressibility correction.
- `src/epcsaft/native/epcsaft_Z.cpp`
  - Has pressure/compressibility logic and an existing CppAD pressure-density derivative helper for `dP/drho`.
- `src/epcsaft/native/epcsaft_density.cpp`
  - Owns density root solving. Existing density validity checks include finite-difference logic, but the new Ceres derivative subsystem must not use finite differences for production derivatives.
- `src/epcsaft/native/epcsaft_regression.cpp`
  - Owns binary VLE regression residual plumbing and Ceres integration points.

### Continuation recommendation

The next implementation should be a bounded native derivative-subsystem slice, not a Python fitting change.

Recommended next task:

1. Add a neutral binary `k_ij` derivative evaluator that treats `k_ij` as an independent native derivative variable.
2. Derive `lnphi`, pressure, and density-root sensitivities without finite differences.
3. Wire binary VLE `k_ij` Ceres residuals to that evaluator.
4. Keep association, ionic, and liquid-electrolyte cases gated as `backend_unavailable` until each has real native derivative support.
5. Add focused tests that fail if the binary Ceres path falls back to finite differences or Python-owned objective evaluation.

Likely verification path after implementation:

```powershell
uv run python scripts/build_epcsaft.py --clean --enable-ceres --enable-cppad
uv run python run_pytest.py tests/native/test_ceres_binary_regression.py tests/native/test_ceres_pure_regression.py -q
uv run python run_pytest.py tests/regression/test_ethanol_water_binary_vle_regression.py tests/api/test_regression_api.py::test_fit_binary_pair_vle_kij_default_and_rejects_temperature_models -q
uv run python scripts/validate_project.py quick
```

### Current stop point

This issue should remain open. Do not open a `Closes #66` PR from the current state unless the full binary derivative path is implemented and validated, or the issue is explicitly split/narrowed.

For the next agent/thread: continue from the local GoalBuddy state under `docs/goals/issue-66-ceres-regression-backend/`, with T012 recording the pause/scope decision. The user did authorize the larger derivative-subsystem direction, but then immediately asked to pause and document the blocker first.
