# ePC-SAFT Equilibrium V4 Handoff

## Active Worktree

Use this checkout for the next equilibrium thread:

```text
C:\Users\Tanner\Documents\git\ePC-SAFT-equilibrium-v3
```

Current branch:

```text
codex/equilibrium-v4-electrolyte-lle
```

This worktree was moved out of the repo-local `.worktrees` folder and recreated as a sibling of `C:\Users\Tanner\Documents\git\ePC-SAFT`. The dirty v4 files were copied from the old checkout and verified byte-for-byte with SHA256 before the old Git worktree registration was removed. Git `safe.directory` was also updated for the sibling path.

One leftover cleanup note: the old path may still contain an empty directory at:

```text
C:\Users\Tanner\Documents\git\ePC-SAFT\.worktrees\ePC-SAFT-equilibrium-v3
```

Windows would not delete it during the move because the running Codex session still had a handle there. It should be safe to remove after no process is using it.

## What Changed In The Equilibrium Workflow

### V3 API Contract Stabilization

The equilibrium API was tightened before the v4 electrolyte work:

- Skipped stability now reports that stability was not checked instead of implying the feed is physically stable.
- `EquilibriumPhase.ln_fugacity_coefficient` was added as the explicit natural-log fugacity coefficient field.
- `EquilibriumPhase.fugacity_coefficient` remains as a backwards-compatible alias for the same natural-log values.
- `to_dict()` emits both keys for compatibility.
- Equilibrium diagnostics and JSON serialization tests were expanded.
- `scripts/sync_equation_registry.py --check` now gives a clear missing-docs-submodule message instead of a traceback when `docs/latex/equations.tex` is absent.

### V4 Modular Electrolyte Equilibrium Direction

The v4 work added a Python-first modular structure for electrolyte LLE work while keeping neutral VLE/LLE behavior backward-compatible. New internal helpers were introduced under:

```text
src/epcsaft/equilibrium_core/
```

Important new files include:

- `classify.py` - route classification helpers.
- `electrolyte_basis.py` - transformed electrolyte basis helpers.
- `thermo_diagnostics.py` - fixed-phase thermodynamic diagnostics for electrolyte LLE benchmarks.

The public API remains centered on `mix.equilibrium(...)`, including electrolyte paths such as `kind="electrolyte_lle"` and salt-molality feed support.

## Main Thermodynamic Diagnosis Result

A key blocker was traced to parameter dataset parity, not only the nonlinear solver.

The old PC-SAFT Khudaida package cache used PC-style `assoc_scheme=2B` for associating solvents. The current ePC-SAFT parameter CSVs had some associating species using `assoc_scheme=1`, which caused the fixed-phase fugacity residuals to disagree badly with the cached legacy/package tie-lines.

The following datasets were corrected to use the expected association workflow:

```text
data/epcsaft_parameters/2022_Ascani/pure/any_solvent.csv
data/epcsaft_parameters/2026_Khudaida/pure/any_solvent.csv
```

Current rule:

- Associating solvents with `e_assoc` / `vol_a` use `assoc_scheme=2B`.
- Ions and non-associating species keep `assoc_scheme` blank.
- Do not add unsupported PC-only fields like `dipm` / `dip_num` to ePC parameter files.

A new dataset contract test guards this:

```text
tests/api/test_parameter_dataset_contracts.py
```

## Khudaida Diagnostic Harness

A new diagnosis-first benchmark layer was added around the existing Khudaida 2026 artifacts:

```text
tests/equilibrium/test_electrolyte_thermo_diagnostics.py
src/epcsaft/equilibrium_core/thermo_diagnostics.py
```

It evaluates fixed aqueous/organic tie-lines and reports:

- neutral component log-fugacity residuals;
- mean ionic NaCl residual;
- phase charge balance;
- reconstructed feed material balance;
- Gibbs feed/split/delta;
- density branch and fugacity contribution snapshots, including Born, Debye-Huckel, and dielectric-related terms where available.

After the `2B` association correction, the cached package-generated Khudaida tie-line and recomputed current ePC-SAFT residuals are aligned at very low residuals. This strongly indicates the thermodynamic surface is now consistent with the cached package result for that diagnostic case.

The v4 predictive solver still does not reliably recover the known Khudaida split from the feed. The fixed-phase thermodynamics are now good enough that this points back toward solver seeding/globalization/recovery for Khudaida, not an obvious fugacity-value mismatch.

## Ascani Electrolyte LLE State

The Ascani 2022 water/butanol/NaCl/KCl path improved after the same association-scheme fix.

The electrolyte LLE tests now assert predictive acceptance gates instead of accepting fixture or partition-seed fallback strings. The important accepted gate is:

```text
predictive_nonlinear_solve
```

The electrolyte tests cover:

- Ascani mixed-salt case;
- direct explicit-ion feeds;
- salt-molality feeds;
- one-salt NaCl smoke behavior;
- non-neutral feed rejection;
- JSON-serializable solver diagnostics.

## Validation Already Run

The last known validation set passed from the sibling-equivalent code state:

```powershell
uv run python scripts\codex_doctor.py
uv run python run_pytest.py tests\api\test_parameter_dataset_contracts.py -q
uv run python run_pytest.py tests\equilibrium\test_electrolyte_lle.py -q
uv run python run_pytest.py tests\equilibrium\test_electrolyte_thermo_diagnostics.py -q
uv run python run_pytest.py tests\equilibrium -q
uv run python run_pytest.py --confidence -q
```

Observed results at that point:

- `tests/api/test_parameter_dataset_contracts.py`: 1 passed.
- `tests/equilibrium/test_electrolyte_lle.py`: 12 passed.
- `tests/equilibrium/test_electrolyte_thermo_diagnostics.py`: 6 passed.
- `tests/equilibrium`: 71 passed.
- `--confidence`: 117 passed, 2 skipped because the private docs submodule was absent.

The direct registry check was also verified to fail cleanly when the docs submodule is absent:

```powershell
uv run python scripts\sync_equation_registry.py --check
```

Expected behavior without `docs/latex/equations.tex`: nonzero exit with a concise missing-submodule message, no Python traceback.

## Known Dirty / Added Files

Important modified files include:

```text
data/epcsaft_parameters/2022_Ascani/pure/any_solvent.csv
data/epcsaft_parameters/2026_Khudaida/pure/any_solvent.csv
pyproject.toml
uv.lock
scripts/sync_equation_registry.py
src/epcsaft/equilibrium.py
src/epcsaft/epcsaft.py
src/epcsaft/_types.py
src/epcsaft/__init__.py
tests/equilibrium/test_api.py
tests/equilibrium/test_lle.py
tests/equilibrium/test_stability.py
tests/equilibrium/test_vle.py
tests/native/test_equation_registry.py
```

Important added/untracked files include:

```text
codex-handoff-equilibrium-v3.md
codex-handoff-equilibrium-v4-electrolyte-lle.md
src/epcsaft/equilibrium_core/__init__.py
src/epcsaft/equilibrium_core/classify.py
src/epcsaft/equilibrium_core/electrolyte_basis.py
src/epcsaft/equilibrium_core/thermo_diagnostics.py
tests/api/test_parameter_dataset_contracts.py
tests/equilibrium/test_electrolyte_lle.py
tests/equilibrium/test_electrolyte_thermo_diagnostics.py
```

Do not revert unrelated dirty files without checking their purpose. This worktree has accumulated changes from multiple equilibrium passes.

## Next Best Steps

1. Start all future work from:

   ```text
   C:\Users\Tanner\Documents\git\ePC-SAFT-equilibrium-v3
   ```

2. Re-run a quick validation from the new location:

   ```powershell
   uv run python scripts\codex_doctor.py
   uv run python run_pytest.py tests\equilibrium -q
   uv run python run_pytest.py --confidence -q
   ```

3. Focus the next v4 solver pass on Khudaida feed-to-split recovery:

   - seed the electrolyte solver from the package-generated fixed tie-line;
   - inspect why the predictive solver rejects that known-good split;
   - improve variable scaling, phase-label handling, TPD/Gibbs seed selection, or least-squares globalization as needed;
   - keep Ipopt/C++ as v5, not v4.

4. Watch one likely small bug: phase labels may be reversed for at least one one-salt molality smoke case. The physics can be right while `aq` / `org` assignment is wrong.

5. Keep generated plot/gallery assets out of this equilibrium pass unless the user explicitly asks for them.
