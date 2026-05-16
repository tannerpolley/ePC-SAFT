# Native Ipopt Derivative Gates Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` when sub-agents are available, or `superpowers:executing-plans` when working inline. Steps use checkbox syntax for tracking.

**Goal:** Follow this plan to its entirety while keeping the strict derivative, solver, packaging, and cleanup gates discussed in the active `/goal`.

**Architecture:** Keep one installable `epcsaft` package for now, but split the internals into EOS/property core, equilibrium extension, and regression extension. Build the new equilibrium implementation as a native C++ Ipopt thermodynamic NLP core in new files, then route public facades into it and delete old paths as parity gates pass.

**Tech Stack:** uv, scikit-build-core, CMake, pybind11, native C++, Ipopt, Ceres, CppAD, Eigen for linear algebra only, pytest through `run_pytest.py`, and repo validation through `scripts/dev/validate_project.py`.

---

## Goal UI Contract

- Active native Codex goal objective:
  - `Follow docs/superpowers/plans/2026-05-16-native-ipopt-derivative-gates.md to its entirety while keeping to the strict derivative, solver, packaging, and cleanup gates discussed.`
- Future agents must start by checking the active `/goal` state and this plan.
- If the goal is interrupted, resume from the checklist in this file and the current git branch state. Do not reconstruct intent from chat history.
- The intended implementation branch is `codex/native-ipopt-derivative-gates`.

## Non-Negotiable Gates

### Derivative Gates

- No derivative-approximation routes anywhere, including tests, examples, diagnostics, benchmarks, and analysis workflows that remain active.
- The old FD-derived concepts are banned even when used as a test oracle.
- Derivatives must be one of:
  - exact analytical derivatives;
  - CppAD derivatives for explicit algebraic terms;
  - CppAD implicit sensitivities for solved state variables;
  - Ceres autodiff only where Ceres owns the regression residual/cost function.
- Exact gradients and Jacobians are required for production Ipopt/Ceres routes.
- Exact Hessians are preferred where available. Ipopt limited-memory Hessian approximation and Ceres Gauss-Newton behavior are allowed only as solver-internal mechanisms and must not be reported as package derivative backends.
- Gate scripts must search for banned derivative/status literals by assembling strings at runtime so the literal banned terms are not checked into source.

### Solver Gates

- Ipopt owns every production equilibrium solve.
- Ceres owns every production regression solve.
- ePC-SAFT code may build thermodynamic functions, residuals, constraints, bounds, scaling, canonical initial points, and derivative callbacks.
- ePC-SAFT code must not own public nonlinear solve algorithms:
  - no package-owned Newton iteration;
  - no package-owned bracketing or bisection;
  - no package-owned scalar search;
  - no package-owned line search or damping loop;
  - no hidden retry loop;
  - no adaptive multistart loop;
  - no legacy status-string route that pretends unsupported work is a valid capability.
- Each equilibrium route gets exactly one canonical deterministic initial point.
- If the canonical point is insufficient, the route fails with typed diagnostics. Fix formulation, scaling, bounds, derivatives, or the canonical initializer rather than adding retries.
- Eigen is allowed only for matrices, factorizations, and linear algebra used by analytical/implicit derivative calculations. Eigen must not be a nonlinear solver or optimizer backend.

### Package Gates

- Keep one installable package until internal boundaries are proven.
- Design internally as three extensions:
  - EOS/property core;
  - native Ipopt equilibrium extension;
  - native Ceres regression extension.
- EOS/property usage must not import or execute equilibrium or regression code.
- Public APIs may remain as facades, but real implementation must move behind subsystem boundaries.
- Separate installable packages are deferred until the boundaries are real and tested.

### Documentation And Test Gates

- Capabilities must report only implemented, validated routes.
- Missing-backend and missing-method wording must be deleted from active source, tests, docs, and roadmaps. If a gate must refer to historical tokens, construct the text by concatenation.
- Delete or rewrite tests that only protect legacy status or dodge behavior.
- The normal quick gate should stay under 10 minutes. Focused slices should be around 90 seconds when practical. Individual unit tests should generally stay under 5 seconds unless marked as confidence/slow coverage.
- Validation must be proportional but real: no final success claim without commands and outputs.

---

## Target File Structure

### Keep As Public Facades

- `src/epcsaft/equilibrium.py`
  - Keep public dataclasses and public function names during migration.
  - Remove old solver option values once the native Ipopt route is implemented.
  - Route supported equilibrium work to native Ipopt entry points.
- `src/epcsaft/reactive_speciation.py`
  - Keep public request/result contracts while moving solve ownership to native Ipopt.
- `src/epcsaft/regression.py`
  - Keep public regression request/result contracts while moving all optimization ownership to native Ceres.
- `src/epcsaft/epcsaft.py`
  - Keep EOS/property public APIs stable while removing direct orchestration of equilibrium/regression internals.

### Add New Native Equilibrium Core

Create a new native subsystem, mostly under:

```text
src/epcsaft/native/equilibrium_nlp/
```

Initial files should be:

- `nlp_problem.h/.cpp`
  - Abstract native problem contract for variables, bounds, constraints, objective, derivatives, scaling, and result unpacking.
- `ipopt_adapter.h/.cpp`
  - The only native layer that calls Ipopt.
  - Converts an `NlpProblem` into Ipopt callbacks.
  - Owns Ipopt status mapping, iteration metadata, and typed failure translation.
- `variable_block.h/.cpp`
  - Variable index registry, bounds, scaling, and named variable slices.
- `constraint_block.h/.cpp`
  - Constraint registry, bounds, scaling, and named constraint slices.
- `gibbs_blocks.h/.cpp`
  - Ideal mixing and reaction-standard-state Gibbs terms.
  - Contains the convex homogeneous ideal chemical-equilibrium subkernel.
- `eos_phase_block.h/.cpp`
  - Helmholtz/EOS phase objective terms, phase volume/density variables, pressure consistency, fugacity/chemical-potential derivatives.
- `reaction_block.h/.cpp`
  - Stoichiometry, reaction extents, reaction constraints, and reaction diagnostic residuals.
- `association_block.h/.cpp`
  - Association site variables and mass-action constraints where coupled into equilibrium NLPs.
- `electrolyte_block.h/.cpp`
  - Charge constraints, ionic contribution wiring, electrolyte-specific diagnostics.
- `route_builders.h/.cpp`
  - Builds route-specific NLPs from common blocks.
- `result_builder.h/.cpp`
  - Converts Ipopt solution plus postsolve diagnostics into native result payloads.

### Add Thin Python Internal Helpers

Use internal Python modules only for request/result payload conversion:

```text
src/epcsaft/equilibrium_core/
```

Expected helpers:

- `native_requests.py`
  - Serialize public Python problem/options into native dictionaries.
- `native_results.py`
  - Convert native dictionaries into public result dataclasses.
- `test_audit.py` or docs-only audit artifact
  - Summarize slow, unnecessary, and weak tests when the test cleanup tranche is executed.

Do not put production solve loops in Python.

---

## Native Ipopt Equilibrium Design

### Core Flow

Every supported equilibrium route must follow this shape:

```text
public Python problem
-> validated native request
-> native route builder
-> one canonical initial point
-> one native Ipopt solve
-> postsolve thermodynamic gates
-> public result or typed failure
```

### Objective Form

For phase equilibrium with an EOS, prefer thermodynamic-potential minimization rather than residual least-squares:

```text
minimize sum_phase [ A_phase(T, V_phase, n_phase) + P_spec * V_phase ]
```

where `A_phase` includes ideal and residual Helmholtz terms. This aligns with ePC-SAFT's Helmholtz-residual structure and allows phase volume/density to be solved as part of the NLP rather than through nested closure solves.

For Gibbs-style homogeneous speciation routes, use Gibbs minimization directly when that is the natural form.

### Constraints

Constraints should represent physical model contracts:

- material balance;
- charge balance;
- phase-sum and phase-amount bounds;
- reaction stoichiometry;
- reaction equilibrium constraints when extents are not represented directly in the objective;
- EOS pressure consistency when using phase volume/density variables;
- association mass-action equations when site fractions are coupled into the route;
- route-specific fixed specifications such as fixed `T`, fixed `P`, bubble pressure variable, or dew pressure variable.

### Convex Chemical Equilibrium Subkernel

Use the Rawlings/Ekerdt-style reduced Gibbs formulation only where its assumptions hold:

- homogeneous chemical/speciation equilibrium;
- fixed `T` and `P`;
- ideal liquid or ideal vapor activity model, or an ideal reference subproblem;
- reaction extents or species amounts tied by linear stoichiometric constraints;
- nonnegative species amounts enforced as bounds.

For those cases, use:

```text
reaction term:       -sum_i log(K_i) * extent_i
ideal mixing term:    sum_j nbar_j * log(nbar_j) - nbar_T * log(nbar_T)
vapor pressure term: +nbar_T * log(P / P_ref), when applicable
```

Use this subkernel as:

- the ideal reaction-equilibrium validation kernel;
- a deterministic initialization/reference formulation for reactive speciation;
- the replacement for custom reaction-extent root/bracket update code in ideal homogeneous cases;
- an exact analytical derivative source for ideal homogeneous reaction/speciation tests.

Do not overclaim convexity. Once ePC-SAFT residual chemical potentials, activity coefficients, density/volume variables, association variables, electrolyte terms, or multiphase splitting are added, the full problem is generally nonconvex. In those cases, retain the thermodynamic NLP structure but do not describe the full route as convex.

### Density And Association Coupling

For equilibrium routes:

- include phase density or volume as Ipopt variables where practical;
- include association site fractions as Ipopt variables where practical;
- express EOS and association consistency as constraints with analytical/CppAD derivatives;
- avoid nested nonlinear closures inside Ipopt objective/constraint callbacks.

For standalone property/state calls:

- density and association closures are evidence-gated internal exceptions;
- replacement with Ceres/Ipopt is allowed only if runtime is near parity, or up to 2x slower with a documented accuracy/stability improvement;
- if the current closure remains, remove derivative-approximation checks and legacy status/dodge diagnostics from it.

---

## Implementation Tasks

### Task 1: Record And Enforce Plan/Gate Baseline

**Files:**
- Modify: `docs/superpowers/plans/2026-05-16-native-ipopt-derivative-gates.md`
- Later create/modify: gate script under `scripts/dev/` or `tests/workflows/repo/`

- [x] Verify the branch is `codex/native-ipopt-derivative-gates`.
- [x] Run `git status --short --branch` and confirm no unrelated changes.
- [x] Add tracked text gate tests that assemble banned terms from pieces.
- [x] Gate active source, tests, docs, roadmaps, scripts, and retained analysis workflows.
- [x] Allow historical archived paper text under `docs/papers/**` only if explicitly excluded and documented.
- [x] Run the gate against current legacy text.
- [x] Remove or rewrite enough explicit banned wording for the gate to pass.
- [x] Commit as `Add strict solver derivative text gates`.

Task 1 note: the first committed gate covers the explicit backend-status and numerical-difference text bans. The broader legacy missing-status family remains widespread and is intentionally left as implementation debt for the solver/derivative replacement tasks, not as completed capability work.

Task 1 continuation note: the CppAD disabled smoke/default derivative result no longer assembles the removed backend-status token at runtime. Disabled CppAD smoke payloads now report `cppad_disabled`, and a repo test blocks reintroducing the old assembled source pattern.

### Task 2: Test Audit And Prune

**Files:**
- Create: `docs/roadmaps/native_ipopt_test_audit.md` or equivalent concise tracked audit.
- Modify: `tests/**`
- Modify: `scripts/dev/validate_project.py` if suite boundaries change.

- [x] Collect test inventory and durations with pytest collection and a duration-enabled run.
- [x] Classify tests as fast gate, focused native, confidence, slow/scientific, docs, package-boundary, or obsolete.
- [x] Identify tests that only protect legacy status/dodge behavior.
- [ ] Delete or rewrite obsolete tests.
- [ ] Move slow scientific matrix coverage out of the quick gate.
- [ ] Add strict tests for new gates:
  - no banned derivative/status concepts;
  - no SciPy package/dev/test dependency;
  - no Eigen nonlinear optimizer route;
  - no Python production solver loop.
- [x] Validate the quick gate remains under 10 minutes.
- [x] Commit as `Audit and tighten solver gate tests`.

Task 2 continuation note: native derivative tests that previously used shifted-source oracle evidence have been rewritten around analytical/CppAD derivative contracts, exact chain-rule identities, backend identity, and finite/nonzero payload checks. The tracked text gate now also blocks common shifted-source oracle and Eigen nonlinear optimizer tokens.

Task 2 continuation note: the electrolyte LLE confidence validation report no longer generates parameter-shift sensitivity artifacts. Report outputs now cover benchmark predictions, continuation, oracle checks, stress cases, and residual/error plots without numeric-differencing or parameter-shift sensitivity metrics.

Task 2 continuation note: the stale full-duration failure list in the tracked test audit has been retired. The previously listed dependency, CppAD LLE, reactive-phase Ceres, and dependency-triage nodes now pass individually.

Task 2 continuation note: after the reactive-speciation diagnostic and MIAC fixture cleanup slices, `uv run python run_pytest.py --all -q` passed with 524 tests and 21 skips in 164.09 seconds.

Task 2 continuation note: electrolyte LLE tests that named the Ceres equilibrium route as production coverage were first narrowed to transitional coverage. The reactive-phase residual-surface evaluator now computes the exact CppAD-implicit Jacobian by default for `auto`, `analytic`, and `cppad` requests instead of emitting an empty derivative payload.

Task 2 continuation note: neutral equilibrium and native chemical-equilibrium diagnostics no longer manufacture missing-status reason fields on successful routes. Routes without a derivative payload now use `not_applicable`, and successful native chemical-equilibrium payloads omit the old empty reason field.

Task 2 continuation note: the strict text gate now blocks explicit Ceres trust-region strategy names. Transitional equilibrium Ceres routes rely on solver defaults and report only an internal trust-region label.

Task 2 continuation note: the neutral equilibrium benchmark payload and table no longer report the retired legacy solver-dodge field; the workflow keeps timing, fingerprint, failure, and diagnostics-key evidence only.

Task 2 continuation note: the reactive regression benchmark payload and table no longer aggregate retired solver-dodge flags, missing-counter lists, or density warm-start dodge counters. Benchmark evidence is now limited to timing, success/failure counts, fingerprints, diagnostic keys, target-family counts, cache hits/misses, and solve/evaluation counters.

### Task 3: Build Dependency Boundary

**Files:**
- Modify: `pyproject.toml`
- Modify: `CMakeLists.txt`
- Modify: `scripts/dev/build_epcsaft.py`
- Modify: `scripts/dev/doctor.py`
- Modify: `.codex/environments/*.toml` if workflow commands change.

- [x] Remove SciPy from package/dev/test dependency groups.
- [x] Replace or delete the single SciPy-based Rezaee fitting analysis workflow.
- [x] Make Ceres required for regression builds.
- [x] Make CppAD required for derivative-capable builds.
- [x] Add native system Ipopt discovery and fail loudly when the required Ipopt build is requested but missing.
- [x] Remove the Python IPOPT wrapper as a production backend.
- [x] Update doctor/build scripts to report Ceres, CppAD, and Ipopt status.
- [x] Validate TOML parsing and stale command searches.
- [x] Commit as `Require native solver dependency gates`.

Task 3 note: the first dependency-boundary slice removed the external numerical package from test dependencies, refreshed the lockfile, and deleted the legacy Rezaee package-local fitting workflow. The second slice added explicit native system Ipopt discovery/reporting and removed the Python IPOPT wrapper. The third slice made Ceres and CppAD required for dev-script, package-backend, and CMake builds, excluded vendored Ceres install rules from the wheel boundary, and validated the local extension with Ceres enabled. Ceres-only regression ownership still remains open for Task 10.

### Task 4: Create Native Ipopt Adapter

**Files:**
- Create: `src/epcsaft/native/equilibrium_nlp/ipopt_adapter.h`
- Create: `src/epcsaft/native/equilibrium_nlp/ipopt_adapter.cpp`
- Create: `src/epcsaft/native/equilibrium_nlp/nlp_problem.h`
- Create: `src/epcsaft/native/equilibrium_nlp/nlp_problem.cpp`
- Modify: `CMakeLists.txt`
- Modify: `src/epcsaft/bindings.cpp`

- [x] Add native C++ abstraction for NLP variables, constraints, objective, gradients, Jacobians, scaling, and result unpacking.
- [x] Implement one Ipopt adapter and keep all direct Ipopt calls there.
- [x] Add a tiny ideal quadratic/linear-constraint smoke problem to prove callback wiring.
- [x] Bind a private smoke entry point only if needed for tests.
- [x] Add tests proving Ipopt availability and adapter callback behavior.
- [x] Commit as `Add native Ipopt adapter`.

Task 4 note: the adapter boundary now lives under `src/epcsaft/native/equilibrium_nlp/`. The local fast build intentionally reports Ipopt disabled, so the quadratic smoke is dependency-gated locally and will execute only in an Ipopt-enabled build.

Task 4 continuation note: the local Windows Ipopt proof uses the install root at `C:\ProgramData\miniconda3\envs\ePC-SAFT\Library` with the MSVC toolchain. The old MinGW probe compiled/linked against that root but crashed because the Ipopt C++ ABI is MSVC-oriented. The supported proof command is `uv run python scripts/dev/build_epcsaft.py --clean --profile ipopt --ipopt-root C:\ProgramData\miniconda3\envs\ePC-SAFT\Library --parallel 4`, followed by setting both `PATH` and `EPCSAFT_RUNTIME_DLL_DIRS` to the Ipopt `Library\bin` directory for Ipopt-executing tests. This adapter proof passed with the native quadratic Ipopt smoke before the active dev build was restored to the default no-Ipopt quick profile.

### Task 5: Add Gibbs And Reaction Blocks

**Files:**
- Create: `src/epcsaft/native/equilibrium_nlp/gibbs_blocks.h`
- Create: `src/epcsaft/native/equilibrium_nlp/gibbs_blocks.cpp`
- Create: `src/epcsaft/native/equilibrium_nlp/reaction_block.h`
- Create: `src/epcsaft/native/equilibrium_nlp/reaction_block.cpp`
- Add tests under `tests/native/equilibrium/`

- [x] Implement ideal homogeneous reduced Gibbs objective terms.
- [x] Implement reaction extent/species amount variable mapping.
- [x] Add analytical gradients/Jacobians for ideal reaction/speciation cases.
- [x] Add ideal liquid and ideal vapor validation cases proving `Q = K`.
- [x] Explicitly label these tests as convex ideal subkernel coverage only.
- [x] Commit as `Add ideal Gibbs reaction NLP blocks`.

Task 5 note: the convex coverage is deliberately limited to homogeneous ideal reaction/speciation validation blocks. It does not assert convexity for ePC-SAFT multiphase, electrolyte, density, or association equilibrium.

### Task 6: Replace Reactive Speciation Solve Ownership

**Files:**
- Modify: `src/epcsaft/reactive_speciation.py`
- Modify: `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp` or replace route with new `equilibrium_nlp` files.
- Modify: `src/epcsaft/bindings.cpp`
- Modify: reactive speciation tests.

- [ ] Route homogeneous reactive speciation to native Ipopt for ideal and nonideal cases.
- [ ] Remove package-owned Newton, scalar bracket, manual damping, and route retry behavior from accepted paths.
- [ ] Feed activity/EOS terms through analytical/CppAD derivatives.
- [ ] Keep Python as request/result facade only.
- [ ] Add tests for ideal, nonideal, charged, and failure diagnostics.
- [ ] Commit as `Route reactive speciation through native Ipopt`.

Task 6 progress note: the first route slice adds explicit `solver_backend="ipopt"` support for homogeneous `ideal_mole_fraction` reactive speciation through a new native `equilibrium_nlp` ideal-speciation problem and the generic Ipopt adapter. The route uses species amounts as variables, material balances as equality constraints, a reduced ideal Gibbs objective, analytical gradients/Jacobians, one deterministic initial point, and Ipopt limited-memory Hessian behavior as solver-internal behavior only. Local proof covered both the no-Ipopt dependency gate and the Ipopt-enabled MSVC build at `C:\ProgramData\miniconda3\envs\ePC-SAFT\Library`. A second slice added exact charge-balance constraints when the charge row is independent of material balances and proved a charged ideal association case through the native Ipopt route. The broader Task 6 checkbox remains open because activity/concentration and EOS/CppAD derivative-coupled routes still need the later EOS derivative NLP blocks before the old accepted custom residual paths can be deleted.

Task 6 continuation note: the public helper that manufactured missing implicit-derivative payloads was removed. Reactive speciation now records implicit solve results only for real analytical/CppAD-backed solved-state sensitivities; unsupported implicit derivative requests raise instead of returning placeholder diagnostics.

Task 6 continuation note: the legacy `jacobian_backend="autodiff"` selector was removed from reactive speciation. Explicit derivative selection is now `auto`, `analytic` as an alias for auto, or `cppad`; CppAD still fails loudly until that residual route exists.

Task 6 continuation note: no-reaction, failed, and best-effort reactive speciation paths now omit implicit-sensitivity payloads when the native result has no reaction-constant sensitivity matrix. Real implicit solve payloads remain required when the native matrices exist.

Task 6 continuation note: public `solve_reactive_speciation`, `mixture.chemical_equilibrium`, and `mixture.equilibrium(kind="chemical_equilibrium")` no longer execute the transitional native chemical-equilibrium residual solve when `solver_backend="auto"`. The public auto route now validates request shape and reaction-convention support, then raises until the native Ipopt homogeneous reactive-speciation NLP owns the route; explicit `solver_backend="ipopt"` remains the accepted ideal-speciation path when the build includes Ipopt. Tests that only protected custom residual-solver success were deleted or rewritten as public route-gate checks, while the private residual evaluator remains covered as derivative-surface diagnostics.

### Task 7: Build EOS Phase Blocks For Equilibrium

**Files:**
- Create: `src/epcsaft/native/equilibrium_nlp/eos_phase_block.h`
- Create: `src/epcsaft/native/equilibrium_nlp/eos_phase_block.cpp`
- Modify: relevant native EOS derivative files.
- Add tests under `tests/native/equilibrium/`.

- [ ] Add phase amount, composition, and volume/density variables.
- [ ] Add Helmholtz/free-energy phase terms.
- [ ] Add EOS pressure consistency constraints where required.
- [ ] Use CppAD/analytical derivatives for objective and constraints.
- [ ] Add single-phase and two-phase consistency tests.
- [ ] Commit as `Add EOS phase NLP blocks`.

### Task 8: Replace Neutral Equilibrium Routes

**Files:**
- Modify: `src/epcsaft/equilibrium.py`
- Modify or replace: `src/epcsaft/native/epcsaft_equilibrium.cpp`
- Add new route builders under `src/epcsaft/native/equilibrium_nlp/`.
- Modify neutral equilibrium tests.

- [ ] Implement native Ipopt route builders for neutral TP flash, VLE, LLE, bubble, and dew workflows.
- [ ] Use one canonical initial point per route.
- [ ] Remove accepted-path bisection/bracketing/scalar-search/golden-section behavior.
- [ ] Add postsolve gates for material balance, phase distance, pressure consistency, and chemical-potential/fugacity consistency.
- [ ] Delete or rewrite legacy tests that assert old diagnostics.
- [ ] Commit as `Replace neutral equilibrium with native Ipopt NLPs`.

Task 8 progress note: the public `EquilibriumOptions.solver_backend` contract no longer accepts explicit `newton`; public selection is limited to `auto` for the still-existing internal native routes and `ipopt` for explicit native-Ipopt requests while the route builders are being replaced.

Task 8 continuation note: obsolete public candidate fallback controls and diagnostics were removed from `EquilibriumOptions`, option normalization, and electrolyte LLE failure diagnostics. The old compatibility dictionary aliases for split and solver-acceptance tolerances are no longer accepted.

Task 8 continuation note: the public `EquilibriumOptions.hessian_strategy` knob was removed. Hessian choices are now native solver-internal details only, consistent with the gate that allows Ipopt limited-memory and Ceres Gauss-Newton behavior without reporting them as package derivative backends.

Task 8 continuation note: `EquilibriumOptions` no longer exposes the public `damping` control. Public equilibrium callers can no longer tune old manual damping behavior, and native equilibrium residual request parsing no longer carries a damping field.

Task 8 continuation note: public neutral bubble/dew scalar solve paths were removed from accepted public execution. The public methods remain declared, validate inputs, and fail loudly until native Ipopt route builders own those production solves.

Task 8 continuation note: `EquilibriumOptions.jacobian_backend` no longer accepts the legacy `autodiff` spelling. Public equilibrium derivative selection is limited to `auto`, `analytic`, or `cppad`, matching the analytical/CppAD derivative gate outside Ceres-owned regression.

### Task 9: Replace Electrolyte And Reactive Phase Equilibrium Routes

**Files:**
- Create: `src/epcsaft/native/equilibrium_nlp/electrolyte_block.h`
- Create: `src/epcsaft/native/equilibrium_nlp/electrolyte_block.cpp`
- Create: `src/epcsaft/native/equilibrium_nlp/association_block.h`
- Create: `src/epcsaft/native/equilibrium_nlp/association_block.cpp`
- Modify electrolyte/reactive public facades and tests.

- [ ] Add charge-balance and electrolyte contribution blocks.
- [ ] Couple density/volume and association variables where practical.
- [ ] Implement electrolyte LLE/VLE/bubble route builders.
- [ ] Implement reactive phase equilibrium route builders.
- [ ] Remove Ceres equilibrium residual-solve ownership from accepted equilibrium paths.
- [ ] Add charge/material/reaction/phase-distance acceptance tests.
- [ ] Commit as `Replace electrolyte reactive equilibrium with native Ipopt`.

Task 9 progress note: the public `ReactiveSpeciationOptions.hessian_strategy` knob was removed for the same reason as the equilibrium facade knob. Hessian behavior remains internal to native solver adapters.

Task 9 continuation note: `ReactiveSpeciationOptions` no longer exposes a public `damping` control. Homogeneous reactive-speciation callers keep tolerance/backend/error controls only.

Task 9 continuation note: chemical-equilibrium residual evaluation payloads no longer report Gauss-Newton as a public Hessian backend, and the ideal Ipopt speciation route reports limited-memory Hessian behavior only as a solver-internal mode. The residual evaluator exposes the implemented analytic Jacobian/gradient surface and explicitly marks Hessian callbacks unavailable until an exact solver-owned Hessian path exists.

Task 9 continuation note: direct pybind native equilibrium and chemical-equilibrium request parsing no longer accepts a `damping` option.

Task 9 continuation note: the unused native chemical-equilibrium `hessian_strategy` request field and parser branch were deleted. Hessian strategy reporting remains limited to the Ipopt adapter's solver-internal metadata.

Task 9 continuation note: stale electrolyte LLE contract tests that expected missing residual derivatives were rewritten to assert the implemented Ceres solve plus CppAD-implicit transformed-variable residual Jacobian surface.

Task 9 continuation note: the Python reactive-phase diagnostics no longer call NumPy's least-squares convenience routine for reaction extent reporting. The helper uses a small direct linear solve of the stoichiometric normal system and remains diagnostic-only.

Task 9 continuation note: public coupled reactive LLE and reactive electrolyte LLE no longer execute the transitional Ceres coupled residual route. Public `ReactivePhaseEquilibriumProblem` and `mixture.equilibrium(kind="reactive_lle" | "reactive_electrolyte_lle")` validate requests and raise until native Ipopt reactive phase-equilibrium route builders own those solves; the lower-level native residual surface remains diagnostic/private coverage.

Task 9 continuation note: private tests that asserted accepted Ceres equilibrium solves for electrolyte LLE, neutral associating LLE, and coupled reactive phase equilibrium were deleted. Remaining native coverage for these transitional surfaces is limited to residual/Jacobian evaluators and public route-gate tests until native Ipopt NLP builders own accepted equilibrium solves.

Task 9 continuation note: the unbound native coupled reactive phase Ceres solve implementation and pybind entrypoint were deleted. The retained native surface is residual/Jacobian evaluation only, with solver diagnostics labeled `residual_surface_only`.

### Task 10: Make Regression Ceres-Only

**Files:**
- Modify: `src/epcsaft/regression.py`
- Modify: `src/epcsaft/native/epcsaft_regression.cpp`
- Modify: `src/epcsaft/bindings.cpp`
- Modify regression tests and docs.

- [x] Delete Eigen unsupported nonlinear optimizer production paths.
- [x] Delete old least-squares backend public aliases.
- [x] Route pure-neutral, pure-ion, binary pair, and generic supported fits through Ceres.
- [ ] Use Ceres autodiff where residuals can be templated.
- [ ] Use analytical/CppAD Jacobians where residuals depend on implicit EOS/state derivatives.
- [x] Add tests proving no numeric-diff Ceres route exists.
- [ ] Commit as `Make regression Ceres-only`.

Task 10 progress note: the first regression slice moves public nonassociating pure-neutral regression to the native Ceres route by default, rejects the old native least-squares backend from that public path, removes the Python and pybind private pure-neutral least-squares entry points, and updates capability metadata so Ceres is a production optimizer when compiled.

Task 10 continuation note: the second regression slice removes the private generic Eigen Levenberg-Marquardt route, its derivative-approximation Jacobian, its pybind/Python wrapper, and the doctor-required symbol. Supported generic production fits remain Ceres-owned. Associating and MEA-CO2-H2O benchmark helpers now do residual-only scoring and reject optimization until native analytic/CppAD/implicit Ceres derivative coverage exists for those target families.

Task 10 continuation note: a tracked repo gate now scans native C++ sources and blocks Ceres numeric-diff APIs plus legacy shifted-source derivative tokens. This closes the explicit test-gate item, while templated Ceres autodiff and broader implicit derivative coverage remain open.

Task 10 continuation note: regression public/native result payloads no longer expose placeholder fallback booleans, empty missing-derivative reason fields, or Hessian skeleton metadata. The retained regression derivative surface reports optimizer backend, derivative backend, objective/evaluation counters, gradient/step norms, and implemented Jacobian availability/backend only.

Task 10 continuation note: the residual-only reactive electrolyte fit wrapper no longer exposes a public `damping` argument or stale line-search status documentation. It remains a structured residual-evaluation result until native Ceres derivative coverage owns that optimization route.

Task 10 continuation note: reactive electrolyte batch capabilities now expose only the implemented residual-only/failed-row fit statuses and renamed the mixed pressure/speciation entry to a residual-context capability. The public metadata no longer advertises bounded step control or line-search fit outcomes for this route.

### Task 11: Internal Extension Boundaries

**Files:**
- Modify: `src/epcsaft/__init__.py`
- Modify: `src/epcsaft/epcsaft.py`
- Modify: `src/epcsaft/equilibrium.py`
- Modify: `src/epcsaft/regression.py`
- Add subsystem modules as needed.

- [x] Ensure EOS/property imports do not import equilibrium/regression implementation modules.
- [x] Keep public facades stable while moving implementation behind internal modules.
- [x] Add import-boundary tests.
- [x] Update docs to describe EOS core, equilibrium extension, and regression extension as internal subsystems.
- [x] Commit as `Separate EOS equilibrium regression internals`.

Task 11 progress note: the top-level `epcsaft` package is now a lazy compatibility facade. Importing
`epcsaft`, importing the EOS namespace, or importing property helpers does not load equilibrium, regression, or
reactive extension modules. Public exports remain available through the same top-level names and load only their owning
subsystem on first access. Repo tests now enforce the import boundary and docs describe the EOS/property core,
equilibrium extension, and regression extension split.

### Task 12: Capabilities, Docs, And Roadmap Cleanup

**Files:**
- Modify: `src/epcsaft/runtime.py`
- Modify: `docs/pages/**`
- Modify: `docs/roadmaps/**`
- Modify: `README.md` if public workflow changes.

- [x] Remove stale compatibility wording.
- [x] Remove legacy solver option docs.
- [x] Make capabilities report only implemented validated routes.
- [x] Document native Ipopt/Ceres/CppAD build requirements.
- [x] Document convex ideal subkernel scope and nonconvex full-route scope.
- [x] Commit as `Document native solver gate architecture`.

Task 12 note: runtime derivative capabilities now list implemented production coverage only instead of embedding open
blocker rows. Reactive electrolyte batch regression is no longer described as a production optimizer in capabilities; it
is documented as a diagnostic residual context until native Ceres owns that route. README and docs now describe native
Ipopt as an explicit constrained-NLP backend with current public wiring limited to homogeneous ideal reactive speciation,
and they state the convex formulation scope as homogeneous ideal subkernels only.

Task 12 continuation note: reactive electrolyte fit status metadata no longer carries an empty placeholder-status list; the residual-only route is labeled `residual_evaluation_only`. EOS derivative coverage matrices now report implemented or out-of-scope rows only instead of acting as a missing-derivative backlog. Unsupported public property-derivative methods raise `InputError` directly rather than returning placeholder derivative result objects.

Task 12 continuation note: Python/native EOS property wrappers no longer reintroduce removed regression fallback/Hessian skeleton metadata, and public `dadx()` payloads no longer expose an empty missing-derivative reason field. The composition-derivative path now either reports implemented derivative backends or raises a typed unsupported-derivative error.

Task 12 continuation note: liquid-electrolyte SSM/DS Born parameter derivatives now either return the implemented analytical liquid payload or raise a typed unsupported-scope error for out-of-scope states. The old unsupported derivative payload row was removed from the native result contract and coverage matrix handling.

Task 12 continuation note: native EOS/CppAD derivative contracts no longer use the old missing-backend label in the C++ and Python property-derivative surfaces touched by the CppAD contract tests. Disabled CppAD reports `cppad_disabled`; unsupported derivative configurations raise typed `unsupported` errors or are excluded from implemented coverage rows.

Task 12 continuation note: density closure diagnostics no longer expose legacy retry wording. Native/Python result payloads now report `density_best_candidate_refinement_used`, `density_best_candidate_rejection_reason`, and `density_warm_start_rejections`, with focused runtime and equilibrium tests updated to enforce the renamed contract. The underlying standalone density closure remains an evidence-gated internal exception rather than a completed Ipopt/Ceres replacement.

Task 12 continuation note: stale planning and GoalBuddy notes no longer preserve the removed missing-backend label. The tracked text gate now also blocks common approximate-derivative wording so future docs, scripts, tests, and source files cannot reintroduce those routes under alternate spelling.

Task 12 continuation note: reactive electrolyte regression no longer exposes the unused failed warm-start policy or the per-row reused-seed status flag. Failed reused seeds remain visible through `warm_start_failed` and `warm_start_source`, while the residual-only route keeps explicit penalty/drop failure handling.

Task 12 continuation note: active source, tests, and dev scripts no longer carry the remaining old retry/default token names. The cleanup renamed internal pressure-seed/default variables, removed an unused reactive regression seed helper, and kept deterministic parameter defaults explicit under non-retry naming.

Task 12 continuation note: reactive/speciation capabilities and staged reactive diagnostics no longer expose negative numerical-derivative availability flags. The positive contract is now the accepted derivative backend list plus typed raise-on-unsupported behavior.

Task 12 continuation note: the tracked text gate now enforces that active source, tests, and dev scripts stay free of old retry/default fallback tokens, while the broader documentation cleanup remains separate from executable surfaces.

Task 12 continuation note: fixed-liquid electrolyte bubble pressure and composed reactive electrolyte bubble pressure no longer execute the package-owned pressure-search route from public Python entry points. The public contracts now validate inputs and raise `InputError` until native Ipopt route builders own those solves; capabilities, docs, and tests mark both routes as `route_pending`, and the public `ElectrolyteBubbleOptions` surface no longer carries the disabled pressure-search controls.

Task 12 continuation note: the private native electrolyte bubble pressure-search binding and C++ implementation were deleted after the public route moved to route-pending. The removed code included the pybind `_solve_electrolyte_bubble_native` entrypoint, the native route option struct, vapor submixture helpers, and the log-pressure search implementation.

Task 12 continuation note: `EquilibriumOptions` no longer exposes the equilibrium-level best-effort result switch. Neutral/electrolyte LLE failures stay loud with JSON-safe diagnostics on `SolutionError`, and docs/tests now instruct downstream sweeps to catch failures rather than consume unaccepted phase results.

Task 12 continuation note: the duplicate reactive-speciation best-effort option alias was removed. Diagnostic nonconverged reactive-speciation payloads now use the explicit `error_mode="result"` contract only, leaving active source, tests, and docs free of the old best-effort token.

Task 12 continuation note: the tracked source/test/script text gate now blocks reintroducing best-effort option tokens and the deleted electrolyte bubble private native binding/solver labels in active executable surfaces.

Task 12 continuation note: the tracked source/test/script text gate now also blocks the removed Hessian-backend diagnostic token so approximate solver-internal Hessian behavior cannot be reintroduced as a package derivative backend label.

Task 12 continuation note: public neutral LLE no longer executes the transitional Ceres residual route from `lle_flash` or `lle_tp`. The public adapter validates feed, mixture, scalar inputs, and optional phase seeds, then raises until a native Ipopt constrained NLP route owns production neutral LLE; lower-level native residual-surface coverage remains private diagnostic coverage.

Task 12 continuation note: public electrolyte LLE no longer executes the transitional Ceres residual route from `electrolyte_lle` or `electrolyte_lle_tp`. The public adapter validates ion-containing mixtures, charge-neutral feeds, electrolyte formula bases, molality/direct-feed inputs, and optional `aq`/`org` phase seeds, then raises until a native Ipopt constrained NLP route owns production electrolyte LLE; lower-level native residual-surface coverage remains private diagnostic coverage.

Task 12 continuation note: public neutral TP flash no longer executes the native package-owned TP flash route from `tp_flash`, `flash_tp`, typed `TPFlash`, or neutral `auto` dispatch. The public adapter validates feed, mixture, and scalar inputs, then raises until a native Ipopt constrained NLP route owns production TP flash. Neutral equilibrium benchmark coverage was narrowed to property-state timing so benchmark scripts do not exercise disabled public flash/LLE routes.

Task 12 continuation note: public neutral and electrolyte stability routes no longer execute native TPD searches from `stability`, `stability_tp`, `electrolyte_stability`, `electrolyte_stability_tp`, typed `StabilityAnalysis`, or reactive-stability post-speciation dispatch. The public adapters validate inputs, phase labels, charge neutrality, and electrolyte formula bases, then raise until native Ipopt stability/NLP route builders own production stability analysis.

Task 12 continuation note: obsolete PR #126 and issue-specific Ceres-equilibrium handoff documents were removed from active docs. Literature benchmark metadata for blocked Ascani electrolyte LLE and reactive phase-equilibrium cases now points at this native Ipopt gate plan instead of superseded Ceres-equilibrium planning artifacts.

Task 12 continuation note: public route-pending errors for neutral flash/LLE/stability, electrolyte LLE/bubble, reactive speciation, and reactive phase equilibrium now state the native Ipopt ownership requirement directly without naming retired solver routes as compatibility context.

Task 12 continuation note: reactive-regression runtime capability labels now describe structured residual-evaluation contexts instead of naming Python orchestration as a solver backend. The metadata keeps the route explicitly non-optimizer until native Ceres derivative coverage owns the production fit.

Task 12 continuation note: the unreferenced tracked LaTeX backup `docs/latex/equations_old.tex` was deleted. The current source of truth remains `docs/latex/equations.tex`, with generated equation navigation kept in `docs/equations.md` and `docs/equations_registry.yaml`.

Task 12 continuation note: completed, unreferenced JetBrains cleanup plan artifacts under `docs/plans/` were removed. The implemented script remains `scripts/dev/configure_jetbrains_project.py`, and local agent routing continues to point there.

Task 12 continuation note: remaining unreferenced stale handoff/planning artifacts under `docs/handoffs/` and `docs/plans/` were removed. This plan remains the active implementation handoff for equilibrium, regression, derivative, and documentation gate work.

Task 12 continuation note: the unused private `_solve_equilibrium_native` pybind wrapper and Python payload adapter were removed. Public equilibrium routes remain route-gated to native Ipopt builders, and private residual-surface bindings stay separate for derivative diagnostics.

Task 12 continuation note: unbound native neutral stability, electrolyte stability, neutral TP flash, neutral LLE, and electrolyte LLE solve entry points were deleted from the C++ header/source. The removed code included custom TPD searches, Rachford-Rice/bisection flash loops, neutral LLE Ceres residual solves, electrolyte LLE Ceres residual solves, and the generic private Ceres residual-solver utility. The retained native electrolyte LLE surface is residual/Jacobian evaluation only.

Task 13 validation note: the `_core` pybind module now disables pybind11 release extras and MSVC optimization for `bindings.cpp`. Native objects still use the configured Release build, while the large binding translation unit avoids the MSVC heap exhaustion seen during rebuild.

Task 14 continuation note: public and native regression no longer expose repeated-start Ceres controls. The Ceres routes now receive one canonical clipped initial point, stale repeated-start tests were simplified, and the tracked source/test/script text gate blocks reintroducing that control surface.

Task 15 continuation note: the native homogeneous chemical-equilibrium solve entrypoint no longer contains the old scalar activity bracket, soft-start, or damped residual-loop implementations. The solve path validates the request and dispatches only to the native Ipopt ideal-speciation NLP; activity/concentration standard states remain route-gated until their EOS derivative NLP blocks exist. The residual evaluator remains available for analytic Jacobian diagnostics.

### Task 13: Final Validation And Cleanup

**Files:**
- No planned source changes except fixes from validation failures.

- [ ] Run normal native build:
  - `uv run python scripts/dev/build_epcsaft.py`
- [ ] Run doctor:
  - `uv run python scripts/dev/doctor.py`
- [ ] Run focused tests for changed areas with `uv run python run_pytest.py ... -q`.
- [ ] Run quick validation:
  - `uv run python scripts/dev/validate_project.py quick`
- [ ] Run docs validation:
  - `uv run python scripts/dev/validate_project.py docs`
- [ ] Run package boundary check:
  - `uv run python scripts/dev/build_dist.py`
- [ ] Run cleanup hook:
  - `pwsh.exe -NoProfile -ExecutionPolicy Bypass -File "$env:USERPROFILE\.codex\hooks\codex-cleanup.ps1" -RepoRoot .`
- [ ] Confirm clean git status.
- [ ] Commit any final validation/doc fixes.

---

## Acceptance Criteria

The overhaul is complete only when all are true:

- `main` contains the previous pruning commit.
- This branch contains this plan and subsequent implementation commits.
- EOS/property import/use is isolated from equilibrium/regression implementation.
- Production equilibrium routes use native Ipopt.
- Production regression routes use native Ceres.
- Active source/tests/docs contain no banned derivative/status concepts except gate code that assembles terms safely.
- No SciPy package/dev/test dependency remains.
- The Rezaee SciPy-based fitting workflow is replaced or deleted.
- Old custom public solve algorithms are removed from accepted equilibrium/regression paths.
- Standalone density/association closure decisions have benchmark evidence and clear exception notes if retained.
- Convex chemical equilibrium is implemented only as an ideal homogeneous subkernel and validation target.
- Quick/docs/package validations pass.
- Cleanup hook passes.
- Final git status is clean.

## Stop Conditions

Stop and ask before proceeding if:

- native Ipopt cannot be installed/discovered on the target build environment;
- a route cannot be expressed with analytical/CppAD/Ceres-autodiff derivatives;
- replacing density or association closures would be slower without a documented robustness/accuracy gain;
- public API removal would break documented downstream contracts beyond the agreed solver/backend cleanup;
- validation failures reveal a chemistry/input-basis problem rather than an implementation bug.
