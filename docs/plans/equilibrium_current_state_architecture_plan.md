# Equilibrium Current-State Architecture Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Revise the older equilibrium architecture handoff against the live `epcsaft` package and define the best next architecture for the package as it exists now.

**Architecture:** Keep the current native-backed thermodynamic surface working while replacing the user-facing string-dispatched equilibrium API with explicit methods and splitting the large native equilibrium implementation by responsibility. Treat chemical speciation, fixed-liquid electrolyte bubble pressure, electrolyte LLE, and reactive stability as existing capabilities, not future placeholders.

**Tech Stack:** Python 3.13, `uv`, pytest, scikit-build-core, CMake, pybind11, C++17 native backend, NumPy, pandas.

---

## Evidence Used

This plan is based on the live checkout at commit `0f280dcab93f` in `C:\Users\Tanner\.codex\worktrees\b2b3\ePC-SAFT`, plus the handoff document at `C:\Users\Tanner\Downloads\epcsaft_equilibrium_architecture_handoff.md`.

Verified commands:

```powershell
uv run python -c "import epcsaft, json; print(epcsaft.__version__); print(json.dumps(epcsaft.runtime_build_info(), indent=2)); print(json.dumps(epcsaft.capabilities(), indent=2))"
uv run python run_pytest.py tests/native/test_equilibrium_native_contracts.py tests/native/test_chemical_equilibrium_native.py -q
uv run python scripts/validate_project.py quick
```

Observed results:

- Package version: `1.5.0`.
- Native extension present: `src\epcsaft\_core.cp313-win_amd64.pyd`.
- Native equilibrium and chemical-equilibrium contract tests passed: `22 passed`.
- Quick project validation passed: `22 passed`.
- One focused API/runtime slice failed because `tests/api/test_runtime.py` still expects `electrolyte_bubble_pressure` and `reactive_electrolyte_bubble` to be unavailable, while `epcsaft.capabilities()` now reports both as available. That is a stale contract/docs/test issue, not evidence of a missing native backend.

## What The Handoff Got Right

- The central principle is still correct: expose problem-specific public methods while sharing native thermodynamic primitives internally.
- The native core should remain responsible for density, residual Helmholtz terms, fugacity/activity, TPD, LLE/VLE solves, chemical equilibrium, and performance-critical numerical work.
- The current `equilibrium(kind=...)` dispatcher should remain for compatibility, but it should stop being the primary API.
- Bubble, dew, flash, stability, chemical equilibrium, and reactive equilibrium should reuse common state, residual, composition, and diagnostics machinery.
- Reactive flash must eventually solve phase equilibrium and reaction equilibrium together; the current sequential reactive workflows should not be renamed as rigorous reactive flash.

## What Is Now Stale

The older handoff treats several items as future or placeholder work that are now implemented at least in scoped form:

- `src/epcsaft/equilibrium.py` already has `EquilibriumOptions`, `EquilibriumPhase`, `EquilibriumResult`, `StabilityTrial`, and `StabilityResult`.
- `ePCSAFTMixture.equilibrium(...)` supports `tp_flash`, `auto`, `lle_flash`, `electrolyte_lle`, `electrolyte_bubble_pressure`, `electrolyte_stability`, `stability`, `chemical_equilibrium`, `reactive_stability`, and `reactive_electrolyte_bubble_pressure`.
- `src/epcsaft/electrolyte_bubble.py` implements fixed-liquid electrolyte bubble pressure through `_core._solve_electrolyte_bubble_native`.
- `src/epcsaft/reactive_speciation.py` implements homogeneous activity-coupled reactive speciation through `_core._solve_chemical_equilibrium_native`.
- `src/epcsaft/reactive_electrolyte.py` implements a composed native workflow: chemical speciation followed by fixed-liquid electrolyte bubble pressure.
- `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp` is already split out from the main equilibrium file.
- `src/epcsaft/equilibrium_core/electrolyte_basis.py`, `electrolyte_seeds.py`, `confidence.py`, and `thermo_diagnostics.py` already hold electrolyte workflow helpers and validation diagnostics.

The outdated docs/tests must be updated before broad architecture work, otherwise future agents will keep planning around false negatives.

## Current Architectural Problems

1. `ePCSAFTMixture.equilibrium(...)` is too broad.

It imports and dispatches every equilibrium mode, validates unrelated arguments in repeated branches, and requires callers to remember string mode names. It is useful as compatibility glue, but it is the wrong long-term primary API.

2. `src/epcsaft/native/epcsaft_equilibrium.cpp` is too large.

It is about 128 KB and owns phase-state construction, feed normalization, electrolyte basis logic, TPD, neutral VLE, neutral LLE, fixed-liquid electrolyte bubble pressure, electrolyte LLE, seed generation, acceptance gates, and diagnostics. This is now the highest-risk maintenance bottleneck.

3. Native file-splitting cannot blindly use the old nested layout.

`CMakeLists.txt` currently collects native sources with:

```cmake
file(GLOB EPCSAFT_NATIVE_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/contributions/*.cpp"
)
```

A new `src/epcsaft/native/equilibrium/*.cpp` folder would not be built unless CMake changes too. The safer first split is same-directory shards such as `epcsaft_equilibrium_tpd.cpp`, then move to a subfolder only after the build contract is updated and tested.

4. Neutral bubble/dew is still missing as a public capability.

The package has electrolyte fixed-liquid bubble pressure, but no public neutral `bubble_p`, `bubble_t`, `dew_p`, or `dew_t`. This reverses the original sequencing. The revised priority is to add neutral bubble/dew now, using the mature native phase-state and fugacity machinery.

5. Chemical standard states are not explicit enough.

`ReactionDefinition` has `stoichiometry`, `log_equilibrium_constant`, and `name`, but no explicit `standard_state`. Native diagnostics expose `activity_model`, and tests cover ePC-SAFT activity coupling, but the public reaction object still lets callers omit the standard-state convention.

6. Sequential reactive workflows need clearer names and boundaries.

`reactive_stability` and `reactive_electrolyte_bubble_pressure` are useful workflows, but they are sequential composition paths. They must stay distinct from a future rigorous `reactive_flash_tp`.

7. Runtime metadata, docs, and tests are inconsistent.

`src/epcsaft/runtime.py` says electrolyte bubble pressure is available. `docs/pages/downstream_local_installs.rst` and `tests/api/test_runtime.py` still describe or assert it as unavailable.

## Revised Public API Target

Add explicit methods on `ePCSAFTMixture` as thin facades over existing functions. The implementation should not change numerical behavior in this phase.

Current dispatcher names should remain:

```python
mix.equilibrium(kind="tp_flash", ...)
mix.equilibrium(kind="lle_flash", ...)
mix.equilibrium(kind="electrolyte_lle", ...)
mix.equilibrium(kind="electrolyte_bubble_pressure", ...)
mix.equilibrium(kind="chemical_equilibrium", ...)
mix.equilibrium(kind="reactive_stability", ...)
mix.equilibrium(kind="reactive_electrolyte_bubble_pressure", ...)
```

Add preferred explicit methods:

```python
mix.flash_tp(T, P, z, *, options=None)
mix.lle_tp(T, P, z, *, options=None, initial_phases=None)
mix.stability_tp(T, P, z, *, options=None, parent_phase=None, trial_phases=None)
mix.electrolyte_stability_tp(T, P, *, z=None, solvent_feed=None, salt_molality=None, options=None)
mix.electrolyte_lle_tp(T, P, *, z=None, solvent_feed=None, salt_molality=None, initial_phases=None, options=None)
mix.electrolyte_bubble_p(T, *, x_liq=None, z=None, vapor_species=None, volatile_species=None, nonvolatile_species=None, options=None)
mix.chemical_equilibrium(T, P, *, balances, totals, reactions, initial_x=None, z=None, options=None)
mix.reactive_stability_tp(T, P, *, balances, totals, reactions, initial_x=None, z=None, options=None, phase_options=None)
mix.reactive_electrolyte_bubble_p(T, *, P_seed, balances, totals, reactions, initial_x=None, z=None, vapor_species=None, options=None)
```

Add later, after native support:

```python
mix.bubble_p(T, x, *, vapor_species=None, options=None, initial_guess=None)
mix.bubble_t(P, x, *, vapor_species=None, options=None, initial_guess=None)
mix.dew_p(T, y, *, liquid_species=None, options=None, initial_guess=None)
mix.dew_t(P, y, *, liquid_species=None, options=None, initial_guess=None)
mix.reactive_flash_tp(T, P, z, reactions, *, phases=("liquid", "vapor"), options=None, initial_guess=None)
```

Do not add a public `reactive_flash_tp` that only performs chemical equilibrium followed by ordinary flash. If a sequential helper is needed, name it as a staged workflow.

## Method-by-Method Status

| Method or workflow | Current state | Best next action |
| --- | --- | --- |
| `state(...)` | Native-backed and tested. | Keep stable; reuse for diagnostics and smoke tests. |
| `equilibrium(kind="tp_flash")` | Native-backed neutral TP flash. | Add `flash_tp(...)` wrapper and route legacy dispatcher to it. |
| `equilibrium(kind="lle_flash")` | Native-backed neutral LLE. | Add `lle_tp(...)` wrapper; keep `initial_phases` continuation. |
| `equilibrium(kind="stability")` | Native neutral TPD. | Add `stability_tp(...)` wrapper. |
| `equilibrium(kind="electrolyte_stability")` | Native transformed-basis electrolyte TPD. | Add `electrolyte_stability_tp(...)` wrapper. |
| `equilibrium(kind="electrolyte_lle")` | Native electrolyte LLE with charge-constrained transformed basis and continuation helpers. | Add `electrolyte_lle_tp(...)`; continue improving diagnostics, not API churn. |
| `equilibrium(kind="electrolyte_bubble_pressure")` | Native fixed-liquid electrolyte bubble pressure, neutral vapor species only. | Update stale runtime tests/docs; add `electrolyte_bubble_p(...)` wrapper. |
| `solve_reactive_speciation(...)` and `kind="chemical_equilibrium"` | Native homogeneous activity-coupled speciation. | Add `standard_state` to reaction API with backward-compatible default; document current convention. |
| `kind="reactive_stability"` | Sequential chemical equilibrium then native TPD. | Keep, but document as staged reactive stability, not rigorous reactive flash. |
| `solve_reactive_electrolyte_bubble(...)` and `kind="reactive_electrolyte_bubble_pressure"` | Sequential chemical speciation then fixed-liquid electrolyte bubble pressure. | Keep, harden docs/tests, and do not label as rigorous reactive VLE. |
| Neutral `bubble_p`, `bubble_t`, `dew_p`, `dew_t` | Missing. | Implement before true reactive flash. |
| Electrolyte dew or full electrolyte VLE flash | Missing. | Defer until neutral bubble/dew and fixed-liquid electrolyte bubble contracts are stable. |
| Rigorous `reactive_flash_tp` | Missing. | Defer until standard states and reusable residual blocks are explicit. |

## Implementation Plan

### Task 1: Repair Current Capability Contracts

**Files:**

- Modify: `tests/api/test_runtime.py`
- Modify: `docs/pages/downstream_local_installs.rst`
- Review: `docs/pages/api_reference.rst`
- Review: `docs/pages/electrolyte_vle_reactive_workflow.rst`
- Test: `tests/api/test_runtime.py`

- [ ] Update `tests/api/test_runtime.py` so `electrolyte_bubble_pressure` and `reactive_electrolyte_bubble` are expected available with native backend metadata.
- [ ] Update docs that still call those paths placeholders.
- [ ] Preserve the important caveat: reactive electrolyte bubble is staged speciation plus fixed-liquid bubble pressure, not rigorous reactive flash.
- [ ] Run:

```powershell
uv run python run_pytest.py tests/api/test_runtime.py -q
```

Expected: pass.

- [ ] Run:

```powershell
uv run python scripts/validate_project.py quick
```

Expected: pass.

### Task 2: Add Explicit Public Wrappers Without Numerical Changes

**Files:**

- Modify: `src/epcsaft/epcsaft.py`
- Modify: `docs/pages/api_reference.rst`
- Test: `tests/equilibrium/test_api.py`
- Test: `tests/api/test_reactive_speciation.py`
- Test: `tests/api/test_reactive_electrolyte_bubble.py`

- [ ] Add `flash_tp`, `lle_tp`, `stability_tp`, `electrolyte_stability_tp`, `electrolyte_lle_tp`, `electrolyte_bubble_p`, `chemical_equilibrium`, `reactive_stability_tp`, and `reactive_electrolyte_bubble_p` methods to `ePCSAFTMixture`.
- [ ] Make each wrapper call the existing module-level helper or staged workflow.
- [ ] Refactor `equilibrium(kind=...)` branches to call the wrappers after preserving current argument validation and error messages.
- [ ] Add parity tests comparing each wrapper to the legacy dispatcher for the same request.
- [ ] Do not emit deprecation warnings yet. The compatibility path is still actively used in tests and downstream scripts.
- [ ] Run:

```powershell
uv run python run_pytest.py tests/equilibrium/test_api.py tests/api/test_reactive_speciation.py tests/api/test_reactive_electrolyte_bubble.py -q
```

Expected: pass.

### Task 3: Split Python Dispatcher Validation Into Focused Helpers

**Files:**

- Modify: `src/epcsaft/epcsaft.py`
- Optionally create: `src/epcsaft/equilibrium_api.py`
- Test: `tests/equilibrium/test_api.py`

- [ ] Extract repeated "argument only supported for kind X" checks into small helper functions.
- [ ] Keep public error strings stable where existing tests assert them.
- [ ] Keep the wrapper methods readable enough that each method's required and rejected arguments are obvious.
- [ ] Run:

```powershell
uv run python run_pytest.py tests/equilibrium/test_api.py -q
```

Expected: pass.

### Task 4: Split Native Equilibrium In Same Directory First

**Files:**

- Modify: `src/epcsaft/native/epcsaft_equilibrium.cpp`
- Modify: `src/epcsaft/native/epcsaft_equilibrium.h`
- Create: `src/epcsaft/native/epcsaft_equilibrium_types.h`
- Create: `src/epcsaft/native/epcsaft_equilibrium_phase_state.cpp`
- Create: `src/epcsaft/native/epcsaft_equilibrium_tpd.cpp`
- Create: `src/epcsaft/native/epcsaft_equilibrium_neutral_lle.cpp`
- Create: `src/epcsaft/native/epcsaft_equilibrium_electrolyte_lle.cpp`
- Create: `src/epcsaft/native/epcsaft_equilibrium_electrolyte_bubble.cpp`
- Test: `tests/native/test_equilibrium_native_contracts.py`

- [ ] Move shared structs to `epcsaft_equilibrium_types.h`.
- [ ] Move phase-state construction and density-scoped helper calls to `epcsaft_equilibrium_phase_state.cpp`.
- [ ] Move neutral and electrolyte TPD routines to `epcsaft_equilibrium_tpd.cpp`.
- [ ] Move neutral LLE routines to `epcsaft_equilibrium_neutral_lle.cpp`.
- [ ] Move electrolyte LLE routines to `epcsaft_equilibrium_electrolyte_lle.cpp`.
- [ ] Move fixed-liquid electrolyte bubble pressure routines to `epcsaft_equilibrium_electrolyte_bubble.cpp`.
- [ ] Leave public native entrypoints and thin orchestration in `epcsaft_equilibrium.cpp`.
- [ ] Do not create `src/epcsaft/native/equilibrium/*.cpp` in this task because current CMake would not compile nested sources without a build-system change.
- [ ] Run:

```powershell
uv run python scripts/build_epcsaft.py
uv run python run_pytest.py tests/native/test_equilibrium_native_contracts.py tests/equilibrium/test_lle.py tests/equilibrium/test_electrolyte_lle.py tests/equilibrium/test_electrolyte_bubble.py -q
```

Expected: pass.

### Task 5: Add Neutral Bubble/Dew After The Split

**Files:**

- Modify: `src/epcsaft/equilibrium.py`
- Modify: `src/epcsaft/epcsaft.py`
- Modify: `src/epcsaft/bindings.cpp`
- Modify: `src/epcsaft/native/epcsaft_equilibrium.h`
- Modify or create native bubble/dew source from Task 4 split
- Create: `tests/equilibrium/test_bubble_dew.py`
- Modify: `docs/pages/api_reference.rst`

- [ ] Implement `bubble_p(T, x)` first as the smallest useful neutral bubble calculation.
- [ ] Use a pressure variable transform or bracketed pressure search so pressure stays positive.
- [ ] Compute incipient vapor composition only over molecular species in the neutral mixture.
- [ ] Return residual diagnostics containing pressure, vapor composition, fugacity residuals, residual norm, and iteration count.
- [ ] Add `bubble_t(P, x)` only after `bubble_p` is stable.
- [ ] Add `dew_p(T, y)` and `dew_t(P, y)` after bubble tests pass.
- [ ] Validate against a simple binary and against nearby TP flash consistency.
- [ ] Run:

```powershell
uv run python scripts/build_epcsaft.py
uv run python run_pytest.py tests/equilibrium/test_bubble_dew.py tests/equilibrium/test_vle.py tests/native/test_equilibrium_native_contracts.py -q
```

Expected: pass.

### Task 6: Make Chemical Standard States Explicit

**Files:**

- Modify: `src/epcsaft/reactive_speciation.py`
- Modify: `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp`
- Modify: `src/epcsaft/bindings.cpp`
- Modify: `docs/pages/electrolyte_vle_reactive_workflow.rst`
- Modify: `docs/pages/api_reference.rst`
- Test: `tests/native/test_chemical_equilibrium_native.py`
- Test: `tests/api/test_reactive_speciation.py`

- [ ] Add `standard_state: str = "mole_fraction_activity"` to `ReactionDefinition`.
- [ ] Preserve existing behavior when callers omit `standard_state`.
- [ ] Pass the standard-state label through Python diagnostics before changing native math.
- [ ] Reject unknown standard-state labels with `InputError`.
- [ ] Document that current `log_equilibrium_constant` is interpreted against the selected activity convention.
- [ ] Add tests proving existing reactions still solve and diagnostics expose the convention.
- [ ] Run:

```powershell
uv run python run_pytest.py tests/native/test_chemical_equilibrium_native.py tests/api/test_reactive_speciation.py -q
```

Expected: pass.

### Task 7: Extract Reusable Chemical Residual Blocks

**Files:**

- Modify: `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp`
- Modify: `src/epcsaft/native/epcsaft_chemical_equilibrium.h`
- Optionally create: `src/epcsaft/native/epcsaft_reaction_residuals.cpp`
- Optionally create: `src/epcsaft/native/epcsaft_reaction_residuals.h`
- Test: `tests/native/test_chemical_equilibrium_native.py`

- [ ] Extract reaction-affinity residual evaluation from the standalone solver loop.
- [ ] Extract balance residual evaluation from the standalone solver loop.
- [ ] Keep state/activity evaluation centralized so future reactive flash does not duplicate fugacity/activity logic.
- [ ] Preserve all current diagnostics: activity model, state failure counts, residual families, Jacobian fallback reason, and phase equilibrium handoff.
- [ ] Run:

```powershell
uv run python scripts/build_epcsaft.py
uv run python run_pytest.py tests/native/test_chemical_equilibrium_native.py -q
```

Expected: pass.

### Task 8: Formalize Staged Reactive Workflows

**Files:**

- Modify: `src/epcsaft/reactive_electrolyte.py`
- Modify: `src/epcsaft/epcsaft.py`
- Modify: `docs/pages/electrolyte_vle_reactive_workflow.rst`
- Test: `tests/api/test_reactive_electrolyte_bubble.py`
- Test: `tests/native/test_chemical_equilibrium_native.py`

- [ ] Keep `reactive_stability` described as chemical equilibrium then TPD.
- [ ] Keep `reactive_electrolyte_bubble_pressure` described as chemical equilibrium then fixed-liquid electrolyte bubble pressure.
- [ ] Add diagnostics field `reactive_workflow_class = "staged"` to staged workflows.
- [ ] Add docs warning that these are useful workflow approximations or handoffs, not coupled reactive flash.
- [ ] Run:

```powershell
uv run python run_pytest.py tests/api/test_reactive_electrolyte_bubble.py tests/native/test_chemical_equilibrium_native.py -q
```

Expected: pass.

### Task 9: Design, Then Implement Rigorous Reactive Flash

**Files:**

- Create: `docs/plans/reactive_flash_tp_design.md`
- Later modify: `src/epcsaft/native/epcsaft_reactive_flash.cpp`
- Later modify: `src/epcsaft/native/epcsaft_reactive_flash.h`
- Later modify: `src/epcsaft/bindings.cpp`
- Later modify: `src/epcsaft/epcsaft.py`
- Later create: `tests/equilibrium/test_reactive_flash.py`

- [ ] Do not start this until Tasks 5, 6, and 7 are done.
- [ ] Define unknowns: phase amounts, phase compositions, reaction extents or per-phase species amounts.
- [ ] Define residual families: material balances, phase fugacity/chemical-potential equality, reaction affinity, charge balance, phase normalization, nonnegativity.
- [ ] Add a no-reaction reduction test that matches `flash_tp`.
- [ ] Add a single-phase reduction test that matches `chemical_equilibrium`.
- [ ] Add electrolyte charge-balance tests before calling the electrolyte variant complete.

Expected first deliverable: a design doc, not code.

## Documentation Plan

Update docs in this order:

1. Runtime capabilities and downstream install docs.
2. API reference with explicit wrappers.
3. Equilibrium workflow page that separates neutral flash, LLE, electrolyte LLE, fixed-liquid electrolyte bubble pressure, homogeneous reactive speciation, staged reactive stability, staged reactive electrolyte bubble pressure, and future rigorous reactive flash.
4. Developer native architecture page explaining the split native ownership and validation commands.

Docs must avoid saying "placeholder" for implemented native paths. They should instead state exact scope and limitations.

## Validation Ladder

Use the smallest relevant command first:

```powershell
uv run python run_pytest.py tests/api/test_runtime.py -q
uv run python run_pytest.py tests/equilibrium/test_api.py -q
uv run python run_pytest.py tests/native/test_equilibrium_native_contracts.py -q
uv run python run_pytest.py tests/native/test_chemical_equilibrium_native.py -q
```

After native edits:

```powershell
uv run python scripts/build_epcsaft.py
uv run python run_pytest.py tests/equilibrium/test_vle.py tests/equilibrium/test_lle.py tests/equilibrium/test_electrolyte_lle.py tests/equilibrium/test_electrolyte_bubble.py -q
```

Before handing off:

```powershell
uv run python scripts/validate_project.py quick
```

Use the full equilibrium or confidence suites only after solver behavior changes:

```powershell
uv run python run_pytest.py tests/equilibrium -q
uv run python run_pytest.py --confidence -q
```

## Decisions To Keep Explicit

1. Should explicit wrappers be added before native refactoring?

Recommended answer: yes. This gives downstream users a stable API immediately and reduces the need to touch user-facing behavior during native file splitting.

2. Should the native split use a nested `native/equilibrium/` folder immediately?

Recommended answer: no. Split into same-directory `epcsaft_equilibrium_*.cpp` files first because current CMake already picks those up. Move to a nested folder only as a separate build-system task.

3. Should electrolyte bubble/dew wait for neutral bubble/dew?

Recommended answer: electrolyte fixed-liquid bubble pressure already exists, so do not undo it. But full electrolyte dew/VLE should wait until neutral bubble/dew is stable and the volatility/nonvolatile masks are thoroughly tested.

4. Should `reactive_electrolyte_bubble_pressure` become `reactive_flash_tp`?

Recommended answer: no. It is a staged workflow. Keep it, document it, and reserve `reactive_flash_tp` for a coupled solve.

5. Should old `equilibrium(kind=...)` calls be deprecated immediately?

Recommended answer: no. Add explicit wrappers first, keep compatibility, then consider warnings only after docs and downstream scripts have moved.

## Stop Conditions

Stop and re-plan if any of these occur:

- A wrapper parity test changes numerical results.
- Native file splitting changes solver diagnostics or acceptance gates.
- A bubble/dew implementation passes residual tests but fails TP flash consistency in a simple binary case.
- Chemical standard-state changes alter existing reaction solutions without an intentional migration path.
- Reactive workflow naming starts hiding the difference between staged composition and rigorous coupled equilibrium.
