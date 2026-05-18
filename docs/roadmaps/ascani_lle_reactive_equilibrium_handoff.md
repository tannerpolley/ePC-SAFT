# Ascani 2022 Native Ipopt Electrolyte LLE A-C Handoff

Date: 2026-05-18

Official worktree:

```text
C:\Users\Tanner\Documents\git\ePC-SAFT-ipopt-profile
```

This handoff intentionally replaces the earlier side-chat reactive/Rezaee/MEA roadmap at this path. It is an A-C-only implementation contract.

## Objective

Implement and validate only:

- Stage A: public native Ipopt liquid-root electrolyte LLE route proof.
- Stage B: hard TPDF stability certification.
- Stage C: full Ascani 2022 Case Study 2 validation.

Stop after A-C. Do not start Ascani 2023, Rezaee, Yu, MEA, or any Stage D/E/F work in this implementation thread.

## Source Anchors

Use local source assets already available in or beside this repo:

- `C:\Users\Tanner\Documents\git\ePC-SAFT\docs\papers\md\Ascani, Sadowski, Held - 2022 - Calculation of Multiphase Equilibria Containing Mixed Solvents and M.md`
- `C:\Users\Tanner\Documents\git\ePC-SAFT\docs\papers\md\Ascani, Sadowski, Held - 2022 - Supporting Information for Calculation of multiphase equilibria containing mixed solvents and mixed..md`
- `data/reference/multiphase/ascani_case2_model_comparison.csv`, if still present and source-backed.

Verified source facts to preserve:

- Ascani 2022 Case Study 2 is water + 1-butanol + NaCl + KCl at 298.15 K and 1 bar.
- Reported feed mass fractions: water 0.8094, 1-butanol 0.1728, NaCl 0.0054, KCl 0.0124.
- Distributed-ion package species for the full validation case are `H2O`, `Butanol`, `Na+`, `K+`, `Cl-`.
- Ascani uses fixed-`T,P` Gibbs minimization conditions, transformed charge-neutral electrolyte variables, per-phase electroneutrality, neutral fugacity equality, independent salt-pair/mean-ionic fugacity equality, and TPDF stability analysis.
- Exact paper matching is not the first accepted-route gate. A real accepted split with sane compositions, liquid densities, residual closure, and hard stability certification is the first package proof. Paper-match claims are separate.

## Strict Status Contract

Every top-level stage/lane summary must use exactly one of these statuses:

- `accepted_public_native_ipopt`
- `blocked_source_data`
- `blocked_solver`
- `blocked_capability`
- `failed_gate`
- `not_started`

`diagnostic_complete` may appear only as nested evidence. It must never close a top-level stage.

Benchmark registry status is separate from completion status. `EXECUTABLE` in the registry does not by itself mean the lane is complete.

## Banned Completion Evidence

Do not report A-C complete from any of these:

- fast build with `Ipopt=OFF`;
- dependency-gate-only result;
- private-native-only result;
- route-builder-only result;
- monkeypatched payload;
- synthetic/non-literature benchmark in place of Ascani evidence;
- figure-only proof;
- accepted result without hard postsolve gates;
- fallback solver or dodge flag;
- vapor-like density root accepted as liquid LLE;
- diagnostic-only TPDF result after Stage B.

Prefer loud failure over fake default behavior.

## Shared Native Ipopt Gates

Accepted A-C evidence requires all of these:

- Public API execution through `mix.equilibrium(kind="electrolyte_lle", ...)`.
- `solver_backend == "ipopt"`.
- Ipopt compiled and available in the active build.
- `derivative_backend == "cppad_implicit"`.
- Exact objective gradients and exact constraint Jacobians.
- No finite-difference derivative proof.
- `hessian_approximation == "limited-memory"` explicitly retained and reported. Exact Hessians are not required for A-C.
- `density_backend == "liquid_pressure_root"`.
- No independent phase-volume variables in the public electrolyte LLE NLP.
- Liquid density recomputation through the public liquid pressure-state path with relative error `<= 1e-8`.
- Liquid density `>= 1000 mol/m^3`.
- Material balance norm `<= 1e-8`.
- Phase charge balance norm `<= 1e-8`.
- Neutral fugacity residual `<= 1e-7`.
- Mean ionic/salt-pair residual `<= 1e-7`.
- Phase distance `>= 0.1`.
- Minimum phase fraction `>= 1e-4`.
- `ghat_delta < -1e-8`.

Trace-ion policy:

- Always retain raw trace residuals.
- Acceptance uses weighted residuals with `trace_floor = 1e-10` and tolerance `1e-7`.

All tolerances must be recorded in `analysis.yaml`, copied into retained `summary.json`, and asserted by tests. Tests must fail if `analysis.yaml`, `summary.json`, and benchmark registry tolerances disagree.

## Stage A: Liquid-Root Electrolyte LLE Foundation

Use the current source-like H2O/butanol/NaCl case only as internal route proof. It is not full Ascani Case Study 2 validation.

Stage A is complete only when public `electrolyte_lle` proves:

- native Ipopt liquid-root route is the production path;
- transformed salt-pair/formula basis is used;
- phase volumes are not NLP variables;
- both phases use liquid pressure-state roots;
- route-builder tests prove no electrolyte LLE volume variables;
- accepted-route tests prove returned densities match direct `phase="liq"` state evaluations;
- postsolve rejects nonfinite, vapor-like, or low-density roots;
- postsolve rejects material-balance, charge-balance, fugacity-residual, `ghat_delta`, or phase-distance failures;
- dependency-only payloads cannot be accepted.

Required retained diagnostics:

- solver backend, derivative backend, density backend, Hessian approximation;
- variable model;
- species labels/order and charge vector;
- selected salt/counter-ion pairs and basis rank;
- formula feed and transformed variables;
- phase fractions and phase compositions;
- liquid-root densities and derived phase volumes;
- `ln_phi`;
- neutral log-fugacity residuals;
- mean ionic/salt-pair log-fugacity residuals;
- material residual vector and norm;
- phase charge residuals and norm;
- phase distance and minimum phase fraction;
- `ghat_feed`, `ghat_split`, `ghat_delta`;
- Ipopt status, application status, and callback exception string.

Suggested Stage A validation:

```powershell
uv run python run_pytest.py tests/native/equilibrium/test_route_builders.py tests/native/equilibrium/test_electrolyte_lle_residual_surface.py tests/native/equilibrium/test_electrolyte_lle_residual_jacobian.py tests/equilibrium/electrolyte -q
```

## Stage B: Hard TPDF Stability Certification

Stage B is a hard acceptance gate, not a diagnostic layer.

Run electrolyte TPDF stability certification on:

- the feed;
- each final returned phase.

Fixed multistart set:

- feed-like;
- returned opposite phase;
- water-rich;
- organic-rich;
- salt-rich;
- each Stage C deterministic seed endpoint.

Gate:

- Feed is unstable if any fixed-start minimized TPDF value is `< -1e-8`.
- Final phases are accepted only if every fixed-start minimized TPDF value is `>= -1e-8`.

If blocked or failed, retain:

- parent phase label;
- trial composition;
- minimized TPDF value;
- tolerance;
- failure reason;
- whether the blocker is source data, solver convergence, missing capability, or failed thermodynamic gate;
- a negative-trial composition suitable for later phase-count-discovery work.

Stage B is complete only when public electrolyte LLE accepted results include passing hard TPDF certification for the Stage A proof case and the Stage C full case.

## Stage C: Ascani 2022 Full Case Study 2

Run full Ascani 2022 Case Study 2 through the public package API.

Species:

```text
H2O
Butanol
Na+
K+
Cl-
```

Conditions:

```text
T = 298.15 K
P = 1 bar
feed mass fractions:
  water      0.8094
  1-butanol 0.1728
  NaCl      0.0054
  KCl       0.0124
```

The validation lane must normalize the source feed into both mass and ion mole-fraction bases and retain the exact conversion table.

Required folder:

```text
analyses/paper_validation/native/2022_ascani/
  README.md
  analysis.yaml
  data/input/
  data/processed/
  results/electrolyte_lle/
  figures/figure_4/
  figures/table_5_fugacity/
  figures/stability_summary/
  figures/density_summary/
  figures/residual_summary/
  scripts/run_all.py
```

Attempt matrix:

1. Primary paper-era no SSM+DS Born model basis.
2. If that fails or mismatches, run SSM+DS Born comparison basis.

Each basis uses the same named deterministic seed set, materialized before solving and retained in `analysis.yaml` and outputs:

- `source_expected_or_table_seed`: paper/source-backed phase split if available.
- `water_rich_salt_rich_seed`: water/salt-rich aqueous phase, butanol-rich organic phase.
- `butanol_rich_trace_salt_seed`: strong organic split with trace ions in organic phase.
- `balanced_feed_perturbation_seed`: charge-neutral perturbation around the feed.
- `stage_a_accepted_mapped_seed`: Stage A accepted split mapped to Na/K by feed cation ratio.

The implementation must record the exact numeric seed payloads before the solve attempt. Do not choose seeds after seeing which result converges.

Either model basis may satisfy accepted-route completion if all route gates pass. The summary must state which basis passed. It must not imply paper-era basis matched if only SSM+DS passed.

Paper-match claim is separate from accepted-route completion. Paper-match requires source-backed values and:

- phase composition absolute error `<= 0.02`;
- Table 5 `ln(f/bar)` absolute error `<= 0.1`.

Figures are required only when source tables or QA digitization exist. Figure evidence must record one of:

- `source_table`;
- `digitized_with_qa`;
- `blocked_source_data`.

Blocked figures must not be shown as validation evidence.

Stage C is complete only when:

- `scripts/run_all.py` exits 0 from the repo root;
- public API result has `status == "accepted_public_native_ipopt"`;
- all shared native Ipopt gates pass;
- Stage B hard TPDF certification passes;
- retained outputs include phase split, phase compositions, densities, `ln_phi`, residuals, `ghat`, and source comparison tables;
- no banned completion evidence is used.

Suggested Stage C validation:

```powershell
uv run python analyses\paper_validation\native\2022_ascani\scripts\run_all.py
uv run python run_pytest.py tests/workflows/paper_validation/test_ascani_2022_lle_validation.py tests/workflows/benchmarks/test_benchmark_literature_suite.py -q
uv run python scripts/benchmarks/benchmark_literature_suite.py --case ascani_2022_distributed_ion_lle --json build/validation/ascani_2022_lle.json
```

## Required Summary Schema

Each retained top-level summary must include:

- `schema_version`;
- `stage`;
- `lane_id`;
- strict top-level `status`;
- `status_reason`;
- source assets;
- command list;
- resolved tolerances copied from `analysis.yaml`;
- Ipopt/runtime diagnostics;
- derivative diagnostics;
- Hessian approximation diagnostics;
- density diagnostics;
- material/charge/fugacity residuals;
- TPDF stability results after Stage B;
- `ghat_feed`, `ghat_split`, `ghat_delta`;
- retained outputs;
- blockers;
- claim boundary.

## Ipopt Environment

Use the documented Windows Ipopt profile for any accepted Ipopt claim:

```powershell
$env:PATH = "C:\ProgramData\miniconda3\envs\ePC-SAFT\Library\bin;$env:PATH"
$env:EPCSAFT_RUNTIME_DLL_DIRS = "C:\ProgramData\miniconda3\envs\ePC-SAFT\Library\bin"
uv run python scripts/dev/build_epcsaft.py --clean --profile ipopt --ipopt-root C:\ProgramData\miniconda3\envs\ePC-SAFT\Library --parallel 4
uv run python scripts/dev/doctor.py --require-ipopt
```

Normal fast-profile checks may supplement development, but they must not be cited as accepted Ipopt production evidence.

## Final A-C Acceptance Commands

```powershell
uv run python scripts/dev/doctor.py --require-ipopt
uv run python run_pytest.py tests/native/equilibrium/test_route_builders.py tests/native/equilibrium/test_electrolyte_lle_residual_surface.py tests/native/equilibrium/test_electrolyte_lle_residual_jacobian.py tests/equilibrium/electrolyte -q
uv run python analyses\paper_validation\native\2022_ascani\scripts\run_all.py
uv run python run_pytest.py tests/workflows/paper_validation/test_ascani_2022_lle_validation.py tests/workflows/benchmarks/test_benchmark_literature_suite.py -q
uv run python scripts/benchmarks/benchmark_literature_suite.py --case ascani_2022_distributed_ion_lle --json build/validation/ascani_2022_lle.json
uv run python scripts/dev/check_text_gates.py
uv run python scripts/dev/validate_project.py quick
pwsh.exe -NoProfile -ExecutionPolicy Bypass -File "$env:USERPROFILE\.codex\hooks\codex-cleanup.ps1" -RepoRoot .
git status --short --branch
```

## Hard Stop

After Stages A-C are accepted or explicitly blocked, stop. Do not implement Stage D/E/F in the same thread.

The final A-C implementation response must include:

- current branch/status;
- Ipopt setup/build result;
- exact commands run;
- tests and pass/fail counts;
- Stage A status;
- Stage B status;
- Stage C status;
- whether an accepted public native Ipopt liquid-root Ascani 2022 solve was proven;
- remaining blocker, if any;
- explicit note that Stage D/E/F must start in a new thread/agent.
