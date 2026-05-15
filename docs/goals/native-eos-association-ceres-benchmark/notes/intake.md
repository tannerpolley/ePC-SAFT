# Intake Notes - Native EOS, Association Sensitivities, And Ceres Benchmark

This note is the Stage 0 gate required by GitHub issue #114 before source edits.

## Source Authority

- GitHub issue: https://github.com/tannerpolley/ePC-SAFT/issues/114
- Issue title: Complete native EOS/AD, association implicit sensitivities, and associating-binary Ceres regression
- Issue state at board creation: open
- Issue updated at board creation: 2026-05-15T04:15:06Z

## Board Prep Snapshot

- Prep branch: `issue#114`
- Prep HEAD: `255ca07fbce7a2842aaea353fd626ab297af6067`
- origin/main SHA at board creation: `255ca07fbce7a2842aaea353fd626ab297af6067`
- Local main SHA at board creation: `255ca07fbce7a2842aaea353fd626ab297af6067`
- Worktree status before board files: clean
- JetBrains MCP: initially unavailable, then reachable after the IDE was opened. Board-prep `ide_index_status` reported `isDumbMode=false` and `isIndexing=false`.
- JetBrains MCP board-prep file-index check:
  - Found `src/epcsaft/native/epcsaft_ares.cpp`.
  - Did not find `src/epcsaft/native/epcsaft_properties.cpp`.
  - Found `src/epcsaft/native/epcsaft_regression.cpp`.
  - Found `src/epcsaft/regression.py`.
  - Found `docs/pages/parameter_regression.rst`.
  - For `implicit_sensitivity`, found Python/API/test paths but no native file path in the first IDE-index result page.

## Stage 0 Evidence To Fill Before Source Edits

### Current branch

- `issue#114`
- Current HEAD at Stage 0: `836d27332647f2494d5d3ddfd1dbb1a09856074e`
- `origin/issue#114` at Stage 0: `836d27332647f2494d5d3ddfd1dbb1a09856074e`
- Worktree status after `git fetch origin --prune`: clean.

### origin/main SHA

- `255ca07fbce7a2842aaea353fd626ab297af6067`
- Local `main` at Stage 0: `255ca07fbce7a2842aaea353fd626ab297af6067`

### C++ EOS contribution files inspected

Inspected:

- `src/epcsaft/native/epcsaft_ares.cpp`: contribution assembly and scalar-templated EOS contribution code. Key anchors:
  - `epcsaft_ares.cpp:27` and nearby templated scalar helpers.
  - `epcsaft_ares.cpp:641` for `ares_contributions_cpp(...)`.
  - `epcsaft_ares.cpp:683`, `:757`, and `:903` for CppAD recording paths.
  - `epcsaft_ares.cpp:405` and `:683-688` currently reject active association in CppAD EOS recording because association site fractions require a solved state.
  - `epcsaft_ares.cpp:760-765` currently rejects active association in native binary pair-parameter CppAD derivatives.
- `src/epcsaft/native/contributions/epcsaft_contrib_hc.cpp`: hard-chain contribution owner; subagent anchor `:247`.
- `src/epcsaft/native/contributions/epcsaft_contrib_disp.cpp`: dispersion contribution owner; subagent anchor `:41`.
- `src/epcsaft/native/contributions/epcsaft_contrib_assoc.cpp`: association contribution owner; key anchors listed below.
- `src/epcsaft/native/contributions/epcsaft_contrib_ion.cpp`: Debye-Huckel / ionic contribution owner; subagent anchor `:23`.
- `src/epcsaft/native/contributions/epcsaft_contrib_born.cpp`: Born contribution owner; subagent anchor `:166`.
- `src/epcsaft/native/epcsaft_parameter_setup.cpp`: association and pair-parameter setup helpers; subagent anchors `:410` and `:794-798`.
- `src/epcsaft/native/contributions/epcsaft_contrib_internal.h`: native contribution state structs and prototypes; subagent anchor `:26`.
- `src/epcsaft/native/epcsaft_core_internal.h`: split property entrypoint declarations; subagent anchors `:287-294`.

`src/epcsaft/native/epcsaft_properties.cpp` does not exist. The property implementation is split across:

- `src/epcsaft/native/epcsaft_Z.cpp`
- `src/epcsaft/native/epcsaft_thermo.cpp`
- `src/epcsaft/native/epcsaft_mu.cpp`
- `src/epcsaft/native/epcsaft_fugcoef.cpp`
- `src/epcsaft/native/epcsaft_activity.cpp`
- `src/epcsaft/native/epcsaft_density.cpp`
- `src/epcsaft/native/epcsaft_state.cpp`

### Association solver files inspected

Inspected:

- `src/epcsaft/native/contributions/epcsaft_contrib_assoc.cpp`
  - `:14` mass-action site-fraction update.
  - `:94` iterative `XA` solver.
  - `:32` existing `dXA/dT` linear solve.
  - `:58` existing `dXA/dx` linear solve.
  - `:163` density derivative closure terms.
  - `:196` composition derivative closure terms.
  - `:239` assembled `AssociationIntermediateState`.
  - `:261` call to solve association site fractions.
- `src/epcsaft/native/epcsaft_ares.cpp`
  - `:646`, `:1062`, `:1086`, and `:1112` call into association intermediate state paths.
- `src/epcsaft/native/epcsaft_electrolyte.h`
  - `:493` still declares older association derivative signatures; `rg` found actual definitions in `contributions/epcsaft_contrib_assoc.cpp`.
- Theory anchors checked by the equation owner:
  - `docs/latex/equations.tex:974`, `:992`, `:1014`, `:1057`, `:1562`, `:1574`, `:1591`.

Current implementation gap:

- `AssociationIntermediateState` stores `XA`, `dXA_dt`, and `dXA_dx`, but there is no generic `dXA/dtheta` field for parameter sensitivities.
- Existing association sensitivity functions use `double` and `Eigen::MatrixXd` linear solves; there is not yet a scalar-templated or explicit implicit-sensitivity path for the issue #114 parameter set.
- The exact current code path still reports missing required association derivative support through old status wording. Do not copy that old status token into new committed text.

### Native regression files inspected

Inspected:

- `src/epcsaft/native/epcsaft_regression.cpp`
  - `:14-17` Ceres includes behind `EPCSAFT_HAS_CERES`.
  - `:963-1103` pure neutral Ceres cost function and solver ownership.
  - `:1705` binary `k_ij` Ceres validation.
  - `:1717-1737` current binary Ceres restrictions, including a rejection of association rows.
  - `:2255` `BinaryKijCeresCostFunction`.
  - `:2307` `solve_one_binary_kij_ceres_start_cpp`.
  - `:2700` top-level `fit_generic_ceres_cpp` dispatcher.
  - `:2727` reports no native analytic/CppAD/implicit derivative path for a target set that reaches the unsupported path.
- `src/epcsaft/native/epcsaft_electrolyte.h`
  - `:312` generic regression result struct.
  - `:558` native regression function declarations.
- `src/epcsaft/bindings.cpp`
  - `:1293` pybind entrypoint for generic Ceres requests.
  - `:299` generic serializer surface.
  - `:231` pure-neutral serializer comparison point.
- `src/epcsaft/epcsaft.py`
  - `:2880` Python wrapper for `_fit_generic_native_ceres(...)`.
  - `:2892` direct call into `_core._fit_generic_native_ceres(...)`.
- `src/epcsaft/regression.py`
  - `:129` maps optimization names such as `k_ij`, `l_ij`, and `k_hb_ij` to native target kinds.
  - `:1479` builds native target payloads.
  - `:1577` native generic Ceres wrapper boundary.
  - `:1926` `_fit_binary_pair_internal(...)`.
  - `:2032` routes binary `k_ij` fits to native Ceres when `optimizer_backend == "ceres"`.
  - `:2797` `fit_binary_parameters(...)` wraps `fit_binary_pair(...)`.
  - `:2848` public `fit_binary_pair(...)`.
- `src/epcsaft/__init__.py`
  - `:97` exports public regression API concepts including `fit_binary_pair`, `fit_binary_parameters`, `TargetDataset`, and `TargetRow`.
- Tests inspected:
  - `tests/native/ceres/test_ceres_binary_regression.py:8`
  - `tests/native/contracts/test_association_implicit_derivative_contract.py:15`
  - `tests/regression/literature/test_ethanol_water_binary_vle_regression.py:31`
  - `tests/regression/literature/test_literature_binary_kij_regression.py:31`
  - `tests/api/regression/test_regression_api_native_backends.py:311`

### Current failing or skipped associating-binary regression path

Focused live reproducer:

```powershell
uv run python run_pytest.py tests/native/contracts/test_association_implicit_derivative_contract.py tests/native/ceres/test_ceres_binary_regression.py tests/regression/literature/test_ethanol_water_binary_vle_regression.py -q
```

Result:

- `4 passed, 1 failed`
- Failure: `tests/regression/literature/test_ethanol_water_binary_vle_regression.py::test_native_binary_kij_regression_matches_real_ethanol_water_vle_reference_band`
- Failing public call path: `epcsaft.fit_binary_pair(...)` -> `src/epcsaft/regression.py::_fit_binary_pair_internal(...)` -> `_run_native_generic_ceres(...)` -> `src/epcsaft/epcsaft.py::_fit_generic_native_ceres(...)` -> `_core._fit_generic_native_ceres(...)`.
- Failure reason: native Ceres binary `k_ij` regression currently rejects association rows. The exact thrown message begins with the repo's old missing-backend status token, so it is paraphrased here rather than copied.
- Native gate: `src/epcsaft/native/epcsaft_regression.cpp:1734`, reached via `src/epcsaft/regression.py:2032` and `src/epcsaft/epcsaft.py:2892`.
- Existing native Ceres nonassociating binary smoke passes at `tests/native/ceres/test_ceres_binary_regression.py:8`; the failure is specific to associating rows.
- `binary_pair` V1 supports VLE `x/y` rows only; LLE/excess-property substitutes are rejected at `src/epcsaft/regression.py:1289`.

Python API owner-agent also ran:

```powershell
uv run python run_pytest.py tests/regression/literature/test_ethanol_water_binary_vle_regression.py tests/regression/literature/test_literature_binary_kij_regression.py -q
```

Result:

- `2 passed, 2 failed`
- Both failures are the associating ethanol-water binary `k_ij` Ceres path.
- No current `xfail` associating-binary regression tests were found by the API owner-agent.

### Chosen benchmark fixture

Chosen fixture: ethanol + water binary VLE, 100 kPa.

Primary files:

- `tests/fixtures/literature/binary_kij/ethanol_water_jced2021_100kpa.json`
- `data/reference/regression/binary/ethanol_water_jced2021_vle_100kpa.csv`
- `data/reference/regression/binary/ethanol_water_jced2021_source_notes.md`
- `tests/helpers/binary_regression_cases.py`
- `tests/regression/literature/test_ethanol_water_binary_vle_regression.py`
- `tests/regression/literature/test_literature_binary_kij_regression.py`

Why this fixture:

- It is repo-contained.
- It exercises the public `epcsaft.fit_binary_pair(...)` path.
- It is an associating binary, so association site-fraction sensitivities are on the required path.
- It already targets constant `k_ij` and has a literature sanity band.
- It is the issue #114 acceptable substitute because normalized MEA + water binary VLE/excess-property data were not found in the repo during Stage 0.

Limitations:

- It is ethanol + water, not MEA + water.
- It covers `k_ij`; it does not by itself prove `l_ij` or `k_hb_ij`.
- The source notes state the literature value is a broad sanity band, not an exact identity target, because repo pure parameters are not byte-identical to the paper reference parameterization.

Minimum issue #114 benchmark proof should therefore require:

- Ceres optimizer backend.
- Allowed derivative backend label.
- Association implicit sensitivity participation.
- Initial objective greater than final objective.
- `k_ij` movement from the initial value.
- Populated row diagnostics.
- Public `fit_binary_pair(...)` call succeeds on the repo-contained fixture.

## Commands Stage 0 Should Run Or Justify

```powershell
git fetch origin --prune
git status --short
git rev-parse --abbrev-ref HEAD
git rev-parse origin/main
rg --files src/epcsaft/native
rg -n "CppAD|implicit_sensitivity|association|k_ij|k_hb_ij|l_ij|Ceres" src tests docs
uv run python scripts/dev/doctor.py
```

Stage 0 must add any narrower commands it actually used, including focused tests or build checks that reveal the current failing or skipped associating-binary regression path.

## Stage 0 Commands Actually Run

```powershell
git fetch origin --prune
git status --short
git rev-parse --abbrev-ref HEAD
git rev-parse HEAD
git rev-parse origin/main
git rev-parse origin/issue#114
gh issue view 114 --repo tannerpolley/ePC-SAFT --json number,title,state,url,updatedAt,comments
rg --files src/epcsaft/native
rg --files tests data docs src scripts | rg -n "(MEA|mea|water|ethanol|alcohol|binary|vle|lle|regression|ceres|cppad|association|assoc)"
uv run python scripts/dev/doctor.py
uv run python run_pytest.py tests/native/contracts/test_association_implicit_derivative_contract.py tests/native/ceres/test_ceres_binary_regression.py tests/regression/literature/test_ethanol_water_binary_vle_regression.py -q
```

JetBrains MCP:

- `ide_index_status` reported `isDumbMode=false` and `isIndexing=false`.
- `ide_find_file` confirmed `epcsaft_ares.cpp`, `epcsaft_fugcoef.cpp`, `epcsaft_activity.cpp`, `epcsaft_mu.cpp`, `epcsaft_regression.cpp`, `test_ceres_binary_regression.py`, and `test_association_implicit_derivative_contract.py`.
- `ide_find_file` did not find `epcsaft_properties.cpp`, matching filesystem evidence.

Doctor result:

- `epcsaft_core_error: <none>`
- `cppad_status: enabled_available`
- `install_state: current`

Native solver/backend owner also ran:

```powershell
uv run python scripts/dev/build_epcsaft.py --status
uv run python run_pytest.py tests/native/ceres/test_ceres_binary_regression.py -q
uv run python run_pytest.py tests/regression/literature/test_ethanol_water_binary_vle_regression.py -q
uv run python run_pytest.py tests/regression/literature/test_literature_binary_kij_regression.py -q
```

Result:

- Build status had Ceres configured `ON`, CppAD configured `ON`, and `_core` present.
- Native Ceres binary nonassociating smoke passed.
- Both ethanol-water literature regression tests failed on the association-row blocker described above.

## Stage 0 Decision

Proceed to Judge `T002`.

Recommended first implementation scope:

- Start with the pre-Stage 1 hard requirement and the association-row blocker together only if Judge confirms the write set can stay within the issue-owned native derivative/regression paths.
- Keep the chosen benchmark as ethanol + water VLE until repo-contained MEA + water binary data exist.
- Treat the old Eigen AD production route, old backend labels, and old missing-backend status token as cleanup blockers for issue #114 rather than as acceptable transitional output.
- Do not broaden runtime capabilities until the public `fit_binary_pair(...)` benchmark passes with native Ceres and an allowed derivative backend.

## Hard Boundaries

- Source edits are blocked until Stage 0 fills every required issue #114 intake field above.
- Main thread owns clean or repair builds for `_core`.
- Do not broaden package capabilities before executable proof exists.
- Do not add downstream-application public APIs.
- If files outside the issue-owned paths are needed, Judge must explicitly approve the scope expansion before edits.
