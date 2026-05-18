# Ipopt Improvement Plan

> **For future implementation agents:** Treat this as an implementation roadmap, not evidence that any item is complete. Before editing solver behavior, re-check the active branch, current Ipopt build, focused route tests, and any newly added diagnostics.

**Goal:** Improve package-wide native Ipopt robustness, diagnostics, and performance across equilibrium problem families.

**Scope:** Native Ipopt equilibrium NLPs owned by the upstream `epcsaft` package: electrolyte LLE, neutral and reactive LLE, reactive electrolyte LLE, bubble/dew routes, stability routes, and homogeneous/reactive speciation where Ipopt is used.

**Source Basis:** Ipopt project documentation for options, interfaces, output interpretation, and special features; local package implementation in `src/epcsaft/native/equilibrium_nlp/`, `src/epcsaft/native/epcsaft_equilibrium.cpp`, and `src/epcsaft/native/equilibrium/reactive_phase_equilibrium_problem.cpp`.

---

## Scoring Rubric

Scores estimate package-wide impact, not one paper-validation lane.

- `10`: likely improves many equilibrium routes and changes solver evidence quality.
- `8-9`: likely improves multiple hard routes or makes failures much easier to fix.
- `6-7`: useful but route-specific, secondary, or dependent on higher-priority work.
- `5`: mostly observability or configuration hygiene unless paired with other items.

The first implementation tranche should cover only items with score `>= 9.0`.

---

## Ranked Improvement List

| Rank | Score | Improvement | Package-Wide Impact |
|---:|---:|---|---|
| 1 | 10.0 | Real problem scaling, not all-ones `user-scaling` | Better conditioning for all Ipopt routes, especially mixed residual/balance systems. |
| 2 | 9.5 | Exact Hessian support for small/medium NLPs | Better Newton steps and fewer limited-memory failures on tightly coupled equilibrium problems. |
| 3 | 9.2 | Per-iteration diagnostics through Ipopt callbacks | Turns solver failures into actionable traces instead of final-status-only payloads. |
| 4 | 9.0 | Warm starts with primal and multiplier state | Enables continuation across tie-lines, temperatures, pressure sweeps, and staged reactive solves. |
| 5 | 8.8 | Better initial-point strategy: multistart, continuation, phase-seeded starts | Reduces local basin failures across LLE/reactive/stability routes. |
| 6 | 8.5 | Linear solver selection and reporting | Lets users compare sparse linear algebra backends and diagnose platform-specific behavior. |
| 7 | 8.2 | Sparse Jacobian structures instead of dense-by-default structures | Improves performance and memory use as route dimensions grow. |
| 8 | 8.0 | Separate feasibility, dual, complementarity, and acceptable tolerances | Allows stricter or looser acceptance contracts without overloading one `tolerance`. |
| 9 | 7.8 | Barrier strategy controls | Helps hard constrained problems where adaptive versus monotone barrier behavior matters. |
| 10 | 7.6 | Restoration and infeasibility controls plus diagnostics | Makes local infeasibility exits easier to interpret and sometimes avoid. |
| 11 | 7.4 | Ipopt derivative checker exposed through package diagnostics | Provides a native solver-level verification path for exact gradient/Jacobian claims. |
| 12 | 7.2 | Bound handling controls | Gives routes a controlled way to tune interiorization and bound relaxation. |
| 13 | 7.0 | Line-search and watchdog controls | Helps hard nonconvex cases where default globalization stalls. |
| 14 | 6.8 | Dependency detection for linearly dependent constraints | Useful for reaction, charge, and material-balance systems with redundant rows. |
| 15 | 6.6 | Better limited-memory Hessian configuration | Improves the required limited-memory mode while exact Hessians are unavailable or disabled. |
| 16 | 6.4 | Fixed-variable treatment controls | Avoids poor behavior when specs collapse variables or route bounds pin values. |
| 17 | 6.2 | Linear-system scaling and pivot options | Gives more control over sparse factorization robustness. |
| 18 | 5.8 | Option-file or structured option passthrough with retained fingerprint | Makes reproducible solver experiments possible without hard-coding every knob. |
| 19 | 5.5 | Better final-status classification using Ipopt termination causes | Improves user-facing error boundaries and prevents misleading accepted/blocked labels. |
| 20 | 5.0 | Timing and statistics capture | Separates expensive model evaluations, Jacobian work, factorization, and callback overhead. |

---

## Tranche 1: Score >= 9.0

These are the first implementation candidates. They should be planned and implemented in small commits with focused route tests before any broad validation claim.

### 1. Real Problem Scaling

**Score:** 10.0

**Current Gap:** The Ipopt adapter forces `nlp_scaling_method = "user-scaling"`, but important routes return identity scaling. Electrolyte LLE and reactive phase equilibrium currently scale objective, variables, and constraints as `1.0`.

**Primary Owners:**

- `src/epcsaft/native/equilibrium_nlp/ipopt_adapter.cpp`
- `src/epcsaft/native/equilibrium_nlp/nlp_problem.*`
- `src/epcsaft/native/epcsaft_equilibrium.cpp`
- `src/epcsaft/native/equilibrium/reactive_phase_equilibrium_problem.cpp`
- `src/epcsaft/native/equilibrium_nlp/route_builders.cpp`

**Implementation Direction:**

- Add a shared scaling helper for objective, variable, and constraint scales.
- Derive constraint scaling from initial residual/Jacobian row norms where available.
- Keep physically meaningful rows, such as phase-separation and feasibility rows, on deliberate scales.
- Record scaling strategy and scale ranges in route diagnostics.
- Preserve exact-gradient and exact-Jacobian requirements.

**Validation Targets:**

- Contract tests that `NlpScaling` is nontrivial for electrolyte LLE and reactive LLE.
- Native route-builder tests showing scaling vector sizes and positive finite values.
- A/B diagnostics on representative hard cases: one electrolyte LLE, one reactive LLE, one bubble/dew, one speciation or chemical equilibrium route.

**Risk:** Scaling can change solver trajectories and accepted points. Validation must compare route status, residual norms, postsolve acceptance, and thermodynamic consistency, not just pass/fail.

### 2. Exact Hessian Support

**Score:** 9.5

**Current Gap:** The adapter currently requires limited-memory Hessian mode and rejects any request that disables it. Ipopt supports exact Hessians, and many package NLPs are small enough that exact second derivatives may be practical.

**Primary Owners:**

- `src/epcsaft/native/equilibrium_nlp/ipopt_adapter.cpp`
- `src/epcsaft/native/equilibrium_nlp/nlp_problem.*`
- CppAD-backed NLP implementations under `src/epcsaft/native/equilibrium_nlp/`
- Route-specific native files that own objective and constraint derivatives.

**Implementation Direction:**

- Extend `NlpProblem` with optional Hessian structure and Hessian value callbacks.
- Keep limited-memory Hessian as a reportable route option, not the only supported mode.
- Start with one small, well-bounded problem family before generalizing.
- Report `hessian_approximation` and `hessian_backend` in diagnostics for every Ipopt solve.

**Validation Targets:**

- Adapter smoke test for exact Hessian mode on a small analytic problem.
- Focused exact-Hessian comparison against CppAD or analytic reference values.
- Route test showing exact Hessian can be selected and reported without weakening exact gradient/Jacobian gates.

**Risk:** Hessian implementation mistakes are high impact. This work should follow derivative-checker and diagnostic improvements closely if exact Hessian is not implemented first.

### 3. Per-Iteration Diagnostics

**Score:** 9.2

**Current Gap:** The adapter keeps final status, final variables, final constraints, and status strings, but it does not retain Ipopt iteration history. Hard failures such as local infeasibility, restoration failure, or max-iteration exits are therefore hard to diagnose.

**Primary Owners:**

- `src/epcsaft/native/equilibrium_nlp/ipopt_adapter.cpp`
- Python/native result formatting in `src/epcsaft/bindings.cpp`
- Python result diagnostics in `src/epcsaft/equilibrium_core/native_results.py` and `src/epcsaft/equilibrium.py`

**Implementation Direction:**

- Implement an Ipopt intermediate callback.
- Retain bounded iteration history with objective, primal infeasibility, dual infeasibility, complementarity, barrier parameter, step size, line-search count, and restoration-phase marker when available.
- Capture final callback exception context without hiding route failures.
- Add an option to cap history length so default result payloads remain small.

**Validation Targets:**

- Adapter smoke test that a successful solve records at least one iteration row.
- A hard-route test that a rejected solve records the last finite iterate and final infeasibility measures.
- Result payload tests proving diagnostics are present but bounded.

**Risk:** Large diagnostics can bloat result payloads. The default should keep compact history and allow fuller traces only through explicit options.

### 4. Warm Starts With Primal And Multiplier State

**Score:** 9.0

**Current Gap:** The adapter supplies primal initial values only. It rejects bound and constraint multiplier initialization, and route APIs do not expose a solver-state object suitable for continuation across related equilibrium problems.

**Primary Owners:**

- `src/epcsaft/native/equilibrium_nlp/ipopt_adapter.cpp`
- `src/epcsaft/native/equilibrium_nlp/nlp_problem.*`
- Public option/result wrappers in `src/epcsaft/equilibrium.py`
- Route builders for LLE, reactive LLE, bubble/dew, and stability continuation.

**Implementation Direction:**

- Extend the adapter to accept optional starting multipliers.
- Add a compact `continuation_state` payload with primal variables and multipliers when Ipopt returns them.
- Add route-level opt-in warm-start behavior that validates dimension, species order, route family, and fixed specs before reuse.
- Support continuation across ordered sweeps without silently accepting incompatible states.

**Validation Targets:**

- Adapter test proving multiplier starts are accepted when dimensions match.
- Public route test that a second solve can consume a continuation state from a compatible first solve.
- Negative tests for incompatible species order, route kind, variable count, or fixed specifications.

**Risk:** Bad warm starts can make convergence worse. The route must reject incompatible states loudly and report whether a warm start was used.

---

## Tranche 2: High-Value Solver Surface

These items are not first tranche, but they should be designed with Tranche 1 so APIs do not box them out.

### 5. Initial-Point Strategy

**Score:** 8.8

Add route-owned multistart and phase-seeded starts. This should reuse diagnostics from Tranche 1 so each rejected seed is recorded with final status and infeasibility measures.

### 6. Linear Solver Selection

**Score:** 8.5

Expose `linear_solver` through native options, report the selected backend, and reject unavailable solvers with a clear dependency message. Keep default behavior stable until benchmark evidence supports a change.

### 7. Sparse Jacobian Structures

**Score:** 8.2

Replace dense-by-default Jacobian structures where sparsity is known. This matters most for larger reactive and electrolyte systems with structured material, charge, phase, and reaction rows.

### 8. Separate Tolerance Families

**Score:** 8.0

Expose and record feasibility, dual, complementarity, and acceptable tolerances separately. Public acceptance should still use package-level postsolve gates rather than Ipopt status alone.

---

## Tranche 3: Specialized Robustness Controls

### 9. Barrier Strategy Controls

**Score:** 7.8

Expose a controlled subset of barrier strategy options for hard nonlinear equilibrium routes. Retain selected settings in diagnostics.

### 10. Restoration And Infeasibility Controls

**Score:** 7.6

Expose restoration-related diagnostics first, then options. This helps distinguish infeasible equations from bad starts or poor scaling.

### 11. Ipopt Derivative Checker

**Score:** 7.4

Expose derivative-checker mode as an explicit diagnostic run, not normal production behavior. Use it to verify exact-gradient/Jacobian claims on route fixtures.

### 12. Bound Handling Controls

**Score:** 7.2

Expose controlled bound interiorization settings for route owners. Record final bound activity and original-bound violations.

### 13. Line-Search And Watchdog Controls

**Score:** 7.0

Expose only after iteration diagnostics exist, because these settings are hard to interpret without line-search and step traces.

### 14. Constraint Dependency Detection

**Score:** 6.8

Useful for reactive and electrolyte systems where charge, material balance, phase constraints, and reaction rows can become redundant.

### 15. Better Limited-Memory Hessian Configuration

**Score:** 6.6

Improve the current default mode by configuring the nonlinear variable subspace and recording the selected quasi-Newton behavior.

### 16. Fixed-Variable Treatment Controls

**Score:** 6.4

Needed when specs or bounds pin variables. This can avoid poor derivative behavior for fixed dimensions.

### 17. Linear-System Scaling And Pivot Options

**Score:** 6.2

Relevant once linear solver selection is exposed. Do not expose a broad option surface before selected solver availability is detectable.

### 18. Structured Option Passthrough With Fingerprint

**Score:** 5.8

Add a controlled option surface or option-file bridge with retained solver option fingerprint. This is useful for reproducibility but should not become a dodge flag for required routes.

### 19. Better Final-Status Classification

**Score:** 5.5

Translate Ipopt status into package statuses with more context: accepted, acceptable, infeasible, restoration failure, max iterations, callback failure, dependency/configuration failure, and postsolve rejection.

### 20. Timing And Statistics Capture

**Score:** 5.0

Record callback evaluation counts and timing where practical. This separates thermodynamic state cost, derivative cost, linear solve cost, and Python/native transfer overhead.

---

## Suggested Implementation Order For Tranche 1

The score ranking is not the same as the safest implementation order.

1. **Per-iteration diagnostics.** Needed to measure all later changes and explain regressions.
2. **Real problem scaling.** Highest expected improvement once diagnostics can prove what changed.
3. **Warm starts and continuation state.** Builds naturally on diagnostic and route-state payloads.
4. **Exact Hessian support.** Highest derivative risk; implement after diagnostics or with derivative-checker support.

---

## Cross-Route Validation Matrix

Each tranche should be validated against representative routes instead of one paper case.

| Route Family | Minimum Check |
|---|---|
| Electrolyte LLE | Public `mix.equilibrium(kind="electrolyte_lle", ...)` hard case plus retained Ascani/Khudaida diagnostics. |
| Reactive LLE | Public `mix.equilibrium(kind="reactive_lle", ...)` source-backed fixture. |
| Reactive electrolyte LLE | Public `mix.equilibrium(kind="reactive_electrolyte_lle", ...)` fixture with phase-tagged reactions. |
| Bubble/dew | One neutral and one electrolyte bubble/dew route if available under native Ipopt. |
| Homogeneous speciation | `epcsaft.solve_reactive_speciation` when Ipopt-backed paths are used. |
| Stability | TPDF/stability route with rejected-seed diagnostics retained. |

---

## Required Evidence Before Claiming Completion

For any implemented item:

- Focused native tests must pass.
- Public API route proof must exist when the item affects public equilibrium behavior.
- Diagnostics must report whether the feature was used.
- Limited-memory Hessian use must remain explicitly reported when exact Hessian is not used.
- No fallback solver, synthetic payload, monkeypatch-only proof, or diagnostic-only completion may be counted as route capability.
- A quick project validation command should pass, or the final report must name the exact blocker.

---

## Official Ipopt Documentation Anchors

- `https://coin-or.github.io/Ipopt/index.html`
- `https://coin-or.github.io/Ipopt/INTERFACES.html`
- `https://coin-or.github.io/Ipopt/OPTIONS.html`
- `https://coin-or.github.io/Ipopt/OUTPUT.html`
- `https://coin-or.github.io/Ipopt/SPECIALS.html`

