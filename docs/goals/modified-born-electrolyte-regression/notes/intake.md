# Stage 0 Intake - Issue #118

Goal source: <https://github.com/tannerpolley/ePC-SAFT/issues/118>

Issue title: Complete modified Born / SSM / DS electrolyte regression and Figiel-Held benchmark proof.

Branch at intake: `issue#18`

Board-prep commit: `1dcb1fdd3879039ab9bb94b814a4bf022054016f`

Origin/main at intake refresh: `5391b4011f46f6817081a81faec0f3c5f407580d`

Worktree status at intake: clean tracked tree; only generated local GoalBuddy board assets under `.goalbuddy-board/` are untracked.

## Current Build And Test Baseline

Verified:

- `uv run python scripts/dev/doctor.py` passes. The native `_core` extension is present, CppAD is enabled and available, and the install state is current.
- Focused baseline command passed with 8 tests:
  `uv run python run_pytest.py tests/native/cppad/test_cppad_born_ssmds_liquid_derivatives.py tests/native/ceres/test_ceres_liquid_electrolyte_regression.py tests/regression/literature/test_figiel_2025_born_parameter_parity.py tests/regression/electrolyte/test_miac_liquid_electrolyte_regression.py -q`
- The current native Ceres liquid-electrolyte test passes because it checks the public unavailable-result contract, not because issue #118 is implemented.
- JetBrains MCP was probed for this worktree. The index server reported the active worktree project as unavailable; it listed the main checkout instead. Semantic IDE tooling should be used only if this worktree is opened or the user approves working against an indexed equivalent.

## Born / SSM / DS Ownership

Verified:

- Main Born, SSM, and DS contribution ownership is `src/epcsaft/native/contributions/epcsaft_contrib_born.cpp`.
- `BornSSMDSData` stores the shared shell terms in `src/epcsaft/native/epcsaft_core_internal.h`.
- Residual aggregation ownership is `src/epcsaft/native/epcsaft_ares.cpp`; `ares_born` is the term that carries the contribution into total residual Helmholtz energy and related scalar derivative paths.
- Public/native exposure of liquid Born parameter derivatives runs through:
  - `src/epcsaft/native/epcsaft_born_derivatives.cpp`
  - `src/epcsaft/native/epcsaft_state.cpp`
  - `src/epcsaft/bindings.cpp`
  - `src/epcsaft/epcsaft.py`
- Parameter setup and dielectric/radius rules are owned mainly by `src/epcsaft/native/epcsaft_parameter_setup.cpp`.
- Same-sign ionic dispersion suppression is implemented in native parameter setup/core internals and is part of the electrolyte parameter path.

Verified equation anchors:

- `docs/latex/equations.tex` and `docs/equations_registry.yaml` contain the source-backed entries for `d_born_rule`, `f_mix`, `delta_d_born`, `dterm_ssm`, `dterm_ds`, `ares_born`, `born_ares_dadrho`, `born_ares_dxi`, `born_ares_ssmds_dxi`, and `born_ares_dT`.
- `docs/latex/equations.tex` and `docs/equations_registry.yaml` also contain the relevant ion-diameter, relative-permittivity mixing, and same-sign ionic epsilon-rule anchors.

Inference:

- Stage 1 should treat `epcsaft_contrib_born.cpp`, `epcsaft_ares.cpp`, `epcsaft_born_derivatives.cpp`, `epcsaft_parameter_setup.cpp`, and focused CppAD tests as the smallest coherent derivative surface.

Unknown:

- Whether all issue-requested parameter families can be exposed through the existing public derivative payload without broadening the API shape. The existing payload is limited and should not be treated as complete evidence.

## Current Derivative Coverage

Verified:

- `d_born` and `f_solv` are represented in the public liquid Born derivative helper.
- The public liquid Born derivative helper currently reports only `d_born` and `f_solv`.
- `density` derivative for Born is implemented as zero for this contribution.
- `temperature` derivative is implemented through the Born temperature derivative path and parameter temperature-rule helpers.
- `composition` derivatives are implemented in native contribution code for the analytic liquid path.
- The CppAD composition mode for SSM/DS Born has a known unsupported path in parameter setup.
- Relative-permittivity composition derivatives exist through the parameter setup/state path, but relative-permittivity parameter derivatives are not part of the current liquid Born derivative payload.
- Ion diameter, solvated diameter, and ion-solvent dispersion energy participate in native parameter construction, but they are not exposed as explicit liquid Born parameter derivative targets in the current public payload.

Current derivative tests:

- `tests/native/cppad/test_cppad_born_ssmds_liquid_derivatives.py`
- `tests/native/cppad/test_cppad_relative_permittivity_derivatives.py`

Inference:

- Stage 1 must either expand the derivative payload and tests to cover the missing issue-requested parameter families or record a source-backed blocker before Stage 2 claims completion.

## Native Ceres Regression Ownership

Verified:

- `src/epcsaft/native/regression/` does not exist in this checkout. Native regression ownership is currently flat in `src/epcsaft/native/epcsaft_regression.cpp`.
- Python public regression API ownership is `src/epcsaft/regression.py`.
- The active native Ceres electrolyte-adjacent production path is `fit_pure_ion`, which already exercises density, osmotic coefficient, mean ionic activity coefficient, and `s`/`e`/`d_born` target families for pure-ion fits.
- `fit_liquid_electrolyte_parameters` currently returns an unavailable result instead of calling a production native Ceres liquid-electrolyte path.
- `f_solv` is user-facing as an alias, but it is not wired as a native Ceres target kind.
- `relative_permittivity` is present in public row schemas and API tests, but is not wired into the native Ceres generic residual term enum.
- Solvation or transfer Gibbs-energy rows are available as analysis/reference assets, but are not wired into native generic regression row handling.

Current Ceres/native regression tests:

- `tests/native/ceres/test_ceres_liquid_electrolyte_regression.py`
- `tests/native/ceres/test_ceres_pure_regression.py`
- `tests/api/regression/test_regression_api_native_backends.py`
- `tests/regression/literature/test_literature_pure_parameter_regression.py`

Inference:

- Stage 2 should convert the liquid-electrolyte API from an unavailable-result contract into a native Ceres call path only where target-row and target-parameter support is production-backed.
- A useful first Stage 2 slice is likely to reuse the existing native Ceres generic residual machinery, add missing target kinds such as `f_solv`, and add explicit native term support where the fixtures and property functions are already present.

## Figiel / Held Fixtures And Benchmark Assets

Verified Figiel 2025 assets:

- `tests/fixtures/literature/figiel_2025/miac_liquid_electrolyte.json`
- `data/reference/epcsaft_parameters/2025_Figiel`
- `data/reference/MIAC`
- `analyses/paper_validation/native/2025_figiel`
- Existing tests:
  - `tests/regression/literature/test_figiel_2025_born_parameter_parity.py`
  - `tests/regression/electrolyte/test_miac_liquid_electrolyte_regression.py`
  - `tests/regression/electrolyte/test_miac_liquid_electrolyte_parity.py`

Verified Held assets:

- `data/reference/osmotic/water`
- `data/reference/epcsaft_parameters/2014_Held`
- `data/reference/epcsaft_parameters/2012_Held`
- `analyses/paper_validation/native/2014_held`
- `analyses/paper_validation/native/2012_held`

Unknown:

- No Held-specific package regression test under `tests/regression/**` was found during intake. Stage 3 must add a package-owned test or record an exact fixture blocker if the Held assets are insufficient for automated tolerance checks.

Chosen final benchmark set:

- Figiel 2025 modified Born / SSM / DS liquid-electrolyte MIAC fixture and parameter parity.
- Held/Cameretti aqueous electrolyte density, osmotic coefficient, and MIAC assets where package-owned fixtures can be formed from existing reference data.
- Held alcohol/salt mixed-solvent density, osmotic coefficient, and MIAC assets where package-owned fixtures can be formed from existing reference data.

## Stage 1 Proceed Decision

Decision: proceed, but only through a source-backed derivative slice.

Required before Stage 2 claims:

- Liquid-electrolyte derivative support must be validated for every issue-required variable family or a specific family must be recorded as a source-backed blocker that keeps issue #118 open.
- Production derivative evidence must not be synthetic, mocked, or diagnostic-only.
- Documentation must state liquid-electrolyte scope exactly and must not imply vapor electrolyte Born support.
