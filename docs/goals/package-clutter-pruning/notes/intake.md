# Intake Notes - Package Clutter Pruning

This note records the Phase 0 safety baseline required by issue #120 before cleanup edits.

## Source Authority

- GitHub issue: https://github.com/tannerpolley/ePC-SAFT/issues/120
- Checklist comment: https://github.com/tannerpolley/ePC-SAFT/issues/120#issuecomment-4456462579
- Issue title: Prune generated artifacts, runtime clutter, and staged/optional API noise before vertical implementation work
- Issue state at board creation: open

## Board Prep Snapshot

- Prep branch: `codex/package-cleanup`
- Prep HEAD: `8f0ea57`
- Prep observation: local `main`, `origin/main`, and `origin/HEAD` pointed at `8f0ea57` when the board was created.
- Prep worktree status before board files: clean.

## Phase 0 Evidence To Fill

Commands run:

```powershell
git status --short
git fetch origin --prune
git merge --ff-only origin/main
git rev-parse HEAD
git rev-parse origin/main
git branch --show-current
git ls-files | rg "(__pycache__|\.pyc$|\.pyo$|\.pyd$|\.so$|\.dll$|\.dylib$|\.egg-info/|(^|/)build/|(^|/)dist/)"
rg --files src/epcsaft/benchmarks
rg -n "except Exception|except BaseException|except:" src/epcsaft
```

- origin/main SHA: `8f0ea57c8e5fed3b82ced58f8da9b86b31426932`
- current branch: `codex/package-cleanup`
- current HEAD: `e524248909f41f32de904b2237cf06c78ffb1afe`
- branch/base state: `git fetch origin --prune` succeeded and `git merge --ff-only origin/main` reported already up to date.
- worktree state before cleanup edits: clean after the board-prep commit.

### Current tracked artifact scan output

The exact issue #120 scan currently reports these tracked paths:

```text
tests/workflows/build/__init__.py
tests/workflows/build/test_build_backend.py
tests/workflows/build/test_build_epcsaft.py
tests/workflows/build/test_build_epcsaft_script.py
```

Interpretation: these are build workflow tests, not generated build artifacts. Phase 1 or the final audit must not delete them as generated files; it should either refine the artifact scan for package/source artifacts or explicitly justify this false-positive path segment.

### Current `epcsaft.__init__` import/export surface

`src/epcsaft/__init__.py` currently imports and re-exports package runtime objects from:

```text
dataset_validation
electrolyte_bubble
epcsaft
equilibrium
implicit_sensitivity
parameter_schema
parameter_templates
parameters
properties
reactive_electrolyte
reactive_regression
reactive_speciation
reactive_staged
regression
runtime
```

Notable public-surface risks visible in `__all__`:

```text
ElectrolyteBubbleOptions
ElectrolyteBubbleResult
ReactiveElectrolyteBubbleOptions
ReactiveElectrolyteBubbleResult
ReactiveStagedEquilibriumResult
<implicit-sensitivity export embedding the old backend-status phrase>
solve_reactive_electrolyte_bubble
solve_reactive_electrolyte_bubble_sweep
solve_reactive_staged_equilibrium
```

These names overlap issue #120 concerns around top-level staged, optional, or old-status surfaces, including one implicit-sensitivity export whose name embeds the old backend-status phrase. Later phases must decide whether they remain stable public API, compatibility shims, internal paths, or removable exports.

### Current benchmark modules under `src/epcsaft`

```text
src/epcsaft/benchmarks/__init__.py
src/epcsaft/benchmarks/literature.py
src/epcsaft/benchmarks/neutral_equilibrium.py
src/epcsaft/benchmarks/reactive_regression.py
```

### Current top-level staged/optional modules

```text
src/epcsaft/electrolyte_bubble.py
src/epcsaft/ipopt_backend.py
src/epcsaft/reactive_electrolyte.py
src/epcsaft/reactive_staged.py
```

### Current broad exception sites in `src/epcsaft`

```text
src/epcsaft/benchmarks/neutral_equilibrium.py:250
src/epcsaft/benchmarks/reactive_regression.py:100
src/epcsaft/benchmarks/reactive_regression.py:214
src/epcsaft/benchmarks/reactive_regression.py:1045
src/epcsaft/epcsaft.py:1111
src/epcsaft/epcsaft.py:1199
src/epcsaft/epcsaft.py:1226
src/epcsaft/epcsaft.py:1244
src/epcsaft/epcsaft.py:1499
src/epcsaft/epcsaft.py:1518
src/epcsaft/epcsaft.py:1557
src/epcsaft/epcsaft.py:1581
src/epcsaft/ipopt_backend.py:55
src/epcsaft/ipopt_backend.py:218
src/epcsaft/ipopt_backend.py:290
src/epcsaft/equilibrium.py:986
src/epcsaft/equilibrium.py:1525
src/epcsaft/equilibrium.py:1562
src/epcsaft/reactive_regression.py:781
src/epcsaft/reactive_speciation.py:715
src/epcsaft/reactive_speciation.py:897
src/epcsaft/reactive_speciation.py:1035
src/epcsaft/reactive_speciation.py:1444
src/epcsaft/equilibrium_core/thermo_diagnostics.py:380
src/epcsaft/equilibrium_core/thermo_diagnostics.py:444
src/epcsaft/equilibrium_core/thermo_diagnostics.py:456
src/epcsaft/runtime.py:103
src/epcsaft/runtime.py:111
src/epcsaft/runtime.py:140
src/epcsaft/__main__.py:30
```

### Baseline validation

Validation was delegated to a validation-only `command_runner` during Phase 0.

- baseline `uv run python scripts/dev/validate_project.py quick` result: PASS in 5.58 s. Evidence: doctor reported `epcsaft_core_error: <none>` and `install_state: current`; bundled quick pytest slice reported `31 passed in 3.27s`.
- baseline `uv run python scripts/dev/validate_project.py docs` result: PASS in 4.34 s. Evidence: Sphinx reported `build succeeded` and HTML output under `build\docs-html`.
- baseline `uv run python run_pytest.py tests/api/package -q` result: PASS in 1.70 s. Evidence: `6 passed in 1.03s`.
- baseline `uv run python run_pytest.py tests/api/runtime -q` result: PASS in 1.37 s. Evidence: `43 passed in 0.82s`.
- baseline `uv run python run_pytest.py tests/workflows/repo -q` result: FAIL in 0.76 s during collection. Evidence: `tests/workflows/repo/test_dependency_issue_triage.py` imports `yaml`; collection fails with `ModuleNotFoundError: No module named 'yaml'`.

Preliminary classification: the repo-workflow failure is pre-existing because it occurs before cleanup product edits. A targeted search found `yaml` usage only in `tests/workflows/repo/test_dependency_issue_triage.py`; no `yaml`, `PyYAML`, or `ruamel` dependency declaration was found in `pyproject.toml` or `uv.lock` by the Phase 0 search. Judge should decide whether to add a board task for the missing test dependency before or during cleanup, because final #120 validation requires `tests/workflows/repo -q` to pass.

## Hard Boundaries

- No scientific implementation changes.
- No #114-#119 implementation work.
- No inventory-only closure.
- No staged/debug-route proof counted as production completion.
- No exact banned backend or derivative-status tokens in committed text.
