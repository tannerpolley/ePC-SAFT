# ePC-SAFT Equilibrium V4 Handoff: Electrolyte LLE

Branch: `codex/equilibrium-v4-electrolyte-lle`

## Scope

V4 adds a Python-first electrolyte liquid-liquid equilibrium path beside the existing neutral VLE, neutral LLE, and neutral TPD workflows. The target use case is a mixed-solvent electrolyte split such as the Ascani/Sadowski/Held 2022 water + 1-butanol + salt examples, where an aqueous-rich liquid and an organic-rich liquid both contain distributed ions.

This pass intentionally does not touch generated plot/gallery assets and does not migrate phase-equilibrium kernels to C++ yet.

## Reference Algorithm Context

Ascani, Sadowski, and Held, "Calculation of Multiphase Equilibria Containing Mixed Solvents and Mixed Electrolytes", J. Chem. Eng. Data 2022, DOI `10.1021/acs.jced.1c00866`, frames electrolyte phase equilibrium differently from the neutral LLE/VLE helpers in this repo:

- ions are explicit species, not neutral salt pseudo-components at the phase-result level;
- ions may distribute across multiple liquid phases;
- each phase must remain electroneutral;
- the feed may be specified on molality or mole-fraction basis;
- ePC-SAFT advanced supplies species fugacity/activity surfaces for each candidate phase.

The current v4 implementation uses that structure as a Python charge-constrained prototype. Internally it maps explicit ions into neutral 1:1 salt-pair formula variables, expands solved phases back to explicit ions for native ePC-SAFT state evaluation, and keeps the public phase payloads on the explicit-ion basis.

## Implemented In This Pass

- New public equilibrium kind: `kind="electrolyte_lle"` / `kind="electrolyte_lle_flash"`.
- New public stability kind: `kind="electrolyte_stability"`.
- New coordinator route: `kind="auto"` routes ion-containing mixtures to electrolyte LLE.
- New backend: `backend="electrolyte_lle"`.
- New feed helper: `epcsaft.electrolyte_feed_from_molality(mixture, solvent_feed=..., salt_molality=...)`.
- Direct ionic mole-fraction feed support through `z=[...]`.
- Mixed-solvent molality feed support, for example:

```python
result = mix.equilibrium(
    kind="electrolyte_lle",
    T=298.15,
    P=1.013e5,
    solvent_feed={"H2O": 0.58, "Butanol": 0.42},
    salt_molality={"NaCl": 1.0},
)
```

- Aqueous/organic phase labeling: result phases are `aq` and `org`.
- Internal `src/epcsaft/equilibrium_core/` helpers:
  - `classify.py` for equilibrium route selection;
  - `electrolyte_basis.py` for the Ascani-style independent counterion-pair matrix.
- SciPy-backed predictive solve path:
  - project dependency `scipy>=1.13.1`;
  - electrolyte TPD uses `scipy.optimize.differential_evolution` with local `minimize` polishing;
  - Gibbs phase-amount seeding uses `scipy.optimize.minimize_scalar`;
  - final nonlinear residual solve uses `scipy.optimize.least_squares`.
- Charge-constrained formula-basis electrolyte split path:
  - `phase_equilibrium_model: electrolyte_lle_v4_charge_constrained_solve`
  - explicit-ion phase payloads with internal neutral-salt pair bookkeeping
  - `variable_model: ascani_transformed_salt_pairs`
  - `basis_rank`, `e_matrix`, and `salt_pairs` diagnostics
  - no fixture-accepted or partition-accepted result path remains.
- Predictive nonlinear electrolyte solve:
  - phase material balance is enforced by reconstructing the dependent phase from feed, phase fraction, and one independent formula composition;
  - residuals include neutral fugacity equality and mean-ionic/salt-pair fugacity equality;
  - accepted splits require residual, material balance, charge balance, Gibbs decrease, phase distance, and non-boundary phase-fraction gates;
  - rejected attempts raise `SolutionError` with JSON-serializable diagnostics instead of returning a fake split.
- Electrolyte transformed-variable TPD diagnostics:
  - `stability_analysis: electrolyte_tpd`
  - `stability_checked`
  - `stability_min_tpd`
  - `stability_trial_composition`
  - `repeated_stability_iterations`
- Electroneutral phase construction and diagnostics:
  - `charge_balance_error`
  - `phase_charge_balance`
  - `electrolyte_balance_error`
  - `material_balance_error`
  - `solver_residual_norm`
  - `neutral_fugacity_residuals`
  - `mean_ionic_fugacity_residuals`
  - `salt_pair_residuals`
  - `gibbs_feed`, `gibbs_split`, and `gibbs_delta`
  - `solver_seed_name`
  - `solver_method`, `tpd_method`, and `gibbs_seed_method`
  - `acceptance_gate`
  - Ascani 2022 DOI in `algorithm_reference`

## Current Status And Limitations

- This is still Python-first v4, not the native v5 multiphase solver.
- The old non-predictive acceptance paths were removed:
  - no `ascani_case2_fixture_regression` acceptance;
  - no `v4_partition_seed_api_compatibility` acceptance.
- The Ascani case-2 CSV remains reference data for tests and comparison only. It is not used to accept or construct a runtime split.
- With the corrected `2B` association-scheme datasets, the predictive solver now returns real two-liquid splits for the Ascani 2022 direct ionic and salt-molality smoke cases. These are accepted through `acceptance_gate: predictive_nonlinear_solve`; no fixture or partition-seed acceptance path is used.
- The practical next scientific blocker is broad benchmark validation and solver robustness across the full Khudaida/Ascani matrices, not restoring fake acceptance plumbing.
- The local Ascani model fixture has an organic branch that is butanol-enriched relative to the aqueous branch, but not butanol-majority by mole fraction. Tests intentionally pin the local model values instead of forcing the paper-majority expectation.
- The current molality parser supports salts whose cation/anion labels can be uniquely mapped from the mixture species, such as `NaCl` with `Na+` and `Cl-`.
- V4 requires 1:1 salts with one shared anion when multiple salts are present.
- The existing neutral LLE solver still rejects ion-containing mixtures. That is intentional; electrolyte LLE must go through the explicit electrolyte path.
- Association-scheme dataset parity with the older PC-SAFT repository is now guarded by `tests/api/test_parameter_dataset_contracts.py`: associating species in parameter CSVs must use `2B`; non-associating hydrocarbons and ions must leave `assoc_scheme` blank. This fixed both `2022_Ascani` and `2026_Khudaida`.
- Khudaida 2026 parameter parity with the older PC-SAFT repository is now restored for the cached package tie-lines: the ePC-SAFT dataset uses `2B` association for water, ethanol, and butanol, no association sites for the ions, and the diagnostic pressure is `1.0e5 Pa` to match the legacy generation script.

## Khudaida 2026 Thermodynamic Diagnostics

This pass adds a diagnosis-first benchmark harness around the existing Khudaida 2026 artifacts for water + ethanol + isobutanol + NaCl at 293.15 K, 303.15 K, and 313.15 K with 5 wt% and 10 wt% NaCl.

New internal helper:

- `src/epcsaft/equilibrium_core/thermo_diagnostics.py`

New tests:

- `tests/equilibrium/test_electrolyte_thermo_diagnostics.py`

The harness is non-mutating: it reads the cached data under `scripts/paper_validation/2026_Khudaida_analysis` and the parameter dataset under `data/epcsaft_parameters/2026_Khudaida`, expands NaCl formula-basis rows to explicit `Na+`/`Cl-` phases, and recomputes fixed-phase diagnostics through the current local package.

Representative 5 wt% NaCl, 293.15 K, tie-line 1 results in this worktree:

- cached legacy package residual from `model_tielines.csv`: `3.458188586228952e-09`;
- recomputed current fixed-phase residual norm after PC/ePC parameter parity fix: `8.495497638705274e-10`;
- recomputed current Gibbs delta for the cached package split: `-0.1445783714041866`;
- experimental tie-line recomputed residual norm: `3.0478838400364907`;
- decision: `package_fixed_tieline_internally_consistent` for the cached package tie-line and `thermodynamic_surface_differs_from_reference_tieline` for the experimental tie-line;
- solver-gate decision from the package feed: `fixed_tieline_consistent_solver_suspect`.

The full fixture matrix currently loads 39 cached tie-lines with zero charge-balance error and no missing data. Tables 9 and 10 show that the current package AAD is much worse than the paper ePC-SAFT reference: max package grand AAD `0.1851` versus max paper ePC-SAFT grand AAD `0.0357`.

Interpretation: for this benchmark, the fixed-phase thermodynamic surface now agrees with the legacy PC-SAFT cached package split. The remaining v4 failure is in the predictive solver path from the feed: the current solver rejects with `predictive_solve_failed`, `best_failure_reason: nonlinear residual did not converge`, and `solver_residual_norm: 0.11025562180070825` even though the cached fixed split is fugacity-stationary. The next v4 work should seed or parameterize the nonlinear solve so it can recover this known stationary split, then broaden to the full Khudaida matrix.

## Version 5 Direction

V5 should move the expensive and scientifically central pieces into C++:

- charge-constrained phase-composition parameterization;
- electrolyte TPD/stability analysis;
- multiphase Gibbs/fugacity residual solve with electroneutrality constraints;
- C++ IPOPT or another constrained-NLP backend as an optional Gibbs minimization path after v4 thermodynamic-surface parity is validated;
- analytic or native finite-difference Jacobian assembly;
- bulk phase candidate enumeration for aqueous, organic, vapor, and additional liquid phases.

The Python API added here should remain the public control surface while the residual evaluation and nonlinear solve migrate below it.

## Focused Validation

```powershell
uv run python run_pytest.py tests\equilibrium\test_electrolyte_lle.py -q
uv run python run_pytest.py tests\equilibrium\test_electrolyte_thermo_diagnostics.py -q
uv run python run_pytest.py tests\equilibrium -q
```

Expected: focused electrolyte LLE API tests and the full equilibrium suite pass.

Latest local validation in this worktree:

```powershell
uv run python scripts\codex_doctor.py
uv run python run_pytest.py tests\equilibrium\test_electrolyte_lle.py -q
uv run python run_pytest.py tests\equilibrium\test_electrolyte_thermo_diagnostics.py -q
uv run python run_pytest.py tests\equilibrium -q
uv run python run_pytest.py --confidence -q
uv run python scripts\sync_equation_registry.py --check
```

Observed after the Khudaida diagnostic addition:

- doctor passed;
- `tests\equilibrium\test_electrolyte_lle.py -q`: 12 passed;
- `tests\equilibrium\test_electrolyte_thermo_diagnostics.py -q`: 6 passed;
- `tests\equilibrium -q`: 71 passed;
- `--confidence -q`: 117 passed, 2 skipped because the docs submodule is absent;
- direct registry check exited nonzero with the intended missing-submodule message and no traceback.
