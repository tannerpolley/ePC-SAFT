# ePC-SAFT Master Context Map and Roadmap

## Purpose of this document

Every Codex agent working on `tannerpolley/ePC-SAFT` must read this document before planning, coding, reviewing, or merging.

This document explains:

1. What the package is supposed to be.
2. Why the package needs the current backend, derivative, equilibrium, regression, and benchmark work.
3. What downstream projects require.
4. What prior agents got wrong.
5. What “complete” means now.

This document is authoritative over older roadmap language that allowed audit-only closure, inventory-only closure, staged-only closure, diagnostic-only closure, or documented limitations as completion.

`docs/roadmaps/unified_equilibrium_core_algorithm.md` is the companion architecture contract for shared equilibrium-core implementation. It replaces completed Ipopt tranche handoffs, GoalBuddy receipts, and native-route audit reports as the current source of truth for equilibrium solver structure.

---

# 1. Package identity

`epcsaft` is a general-purpose thermodynamic package.

It must provide generic, reusable ePC-SAFT functionality for:

- EOS/state/property evaluation
- fugacity coefficients
- activity coefficients
- chemical potentials
- density and pressure closures
- pure-component and binary parameter regression
- electrolyte parameter regression
- phase equilibrium
- chemical/speciation equilibrium
- reactive phase equilibrium
- electrolyte LLE
- reactive electrolyte LLE
- literature benchmark validation
- downstream integration through generic APIs

It must not become application-specific.

Forbidden public package API shapes:

```python
fit_mea_absorption(...)
fit_lithium_extraction_parameters(...)
screen_lithium_extractants(...)
fit_absorption_column(...)
calculate_extraction_efficiency(...)
calculate_distribution_coefficient(...)
calculate_selectivity(...)
```

Allowed package API concepts:

```python
ParameterSet.from_records(...)
ParameterSet.from_dataset(...)
ParameterSet.to_runtime_dict()
State(...)
TPFlash(...)
LLEProblem(...)
ElectrolyteLLEProblem(...)
ReactiveSpeciationProblem(...)
ReactivePhaseEquilibriumProblem(...)
TargetRow(...)
TargetDataset(...)
TargetDataset.target_family_summaries()
mixture.solve_equilibrium(problem)
mixture.equilibrium(kind="tp_flash", ...)
fit_pure_parameters(...)
fit_binary_parameters(...)
fit_liquid_electrolyte_parameters(...)
build_reactive_regression_objective(...)
evaluate_reactive_regression_objective(...)
epcsaft.capabilities()
```

Downstream projects compute application metrics from generic outputs.

Current API orientation for agents:

- Prefer typed equilibrium problems plus `mixture.solve_equilibrium(problem)` for non-reactive route ownership and result contracts.
- Treat `mixture.equilibrium(kind=...)` as the convenience facade that adapts explicit non-reactive string requests and owns reactive-specialized entrypoints.
- Treat `ParameterSet` as the canonical parameter-family boundary; runtime payload emission belongs to `ParameterSet.to_runtime_dict()`.
- Treat `TargetDataset.target_family_summaries()` as the shared target-family summary shape across generic and reactive regression evidence.
- Treat `epcsaft.capabilities()` and `capability_evidence` as the authoritative public capability contract.

---

# 2. Downstream projects and why this package matters

## MEA-Thermodynamics

Needs a publishable true-species reactive electrolyte ePC-SAFT workflow for MEA/CO2/H2O.

The downstream project needs:

- pure MEA, CO2, H2O parameter handling
- MEA association scheme support
- MEA-water and CO2-related binary parameters
- true-species speciation
- activity-based reaction equilibrium
- vapor/liquid fugacity equality for CO2 and other neutral volatile species
- pressure/speciation validation
- ionic parameter regression for MEAH+ and MEACOO-
- Born/SSM/DS parameter regression where justified
- reproducible manuscript figures and parameter tables

The package must provide generic equilibrium and regression capabilities. The MEA repo owns manuscript data, plots, interpretation, and application-specific wrappers.

## Lithium_Extraction

Needs electrolyte LLE and reactive electrolyte LLE support for brine/organic solvent systems.

The downstream project needs:

- generic electrolyte LLE with distributed ions
- phase electroneutrality
- mixed salts / common ions
- organic and aqueous phase support
- lithium and sodium species in aqueous and possibly organic phases
- reaction-coupled phase transfer when a chemistry model requires it
- generic parameter regression where literature parameters are missing
- output phase compositions, activities, fugacities, and reaction extents

The package must not compute lithium extraction efficiency, selectivity, or distribution coefficient. The downstream repo computes those from generic equilibrium outputs.

## MEA-Absorption-Column

Needs generic reactive vapor/liquid equilibrium and speciation outputs that can be used inside column calculations.

The package provides thermodynamic state/equilibrium results. The column project owns column discretization, transport, kinetics, and absorber metrics.

---

# 3. Scientific model foundation

The core ePC-SAFT model is based on residual Helmholtz energy:

```text
a_res = a_hard_chain + a_dispersion + a_association + a_Debye_Huckel + a_Born
```

The package must therefore treat the EOS as a coupled thermodynamic model, not as isolated property helpers.

Required contribution families:

- hard-chain
- dispersion
- association
- Debye-Hückel / ionic
- Born
- modified Born / SSM / DS liquid-electrolyte terms
- relative-permittivity and concentration-dependent dielectric behavior

The literature basis requires:

- pure neutral parameters: `m`, `sigma`, `epsilon`
- associating parameters: association energy and association volume
- ion parameters: ionic diameter / solvated diameter and ion-solvent dispersion energy
- modified Born parameters: `d_born`, `f_solv`, and related solvent/ion terms
- binary parameters: `k_ij`, `l_ij`, `k_hb_ij`

The package must treat `k_ij`, `l_ij`, and `k_hb_ij` as real required parameter families. They are not optional vocabulary.

---

# 4. Non-negotiable backend policy

Two literal tokens are banned from committed repository text.

Define them only as fragments in scripts/tests to avoid writing the exact strings:

```python
BANNED_MISSING_BACKEND_TOKEN = "backend" + "_" + "unavailable"
BANNED_NONEXACT_DERIVATIVE_TOKEN = "finite" + "_" + "difference"
```

The exact assembled strings must not appear anywhere in committed repository text, including:

- source
- tests
- docs
- prompts
- GoalBuddy files
- issue-scope notes
- capability strings
- diagnostics
- comments
- generated committed artifacts

The package must not use missing-backend status strings as normal runtime output for required workflows.

The package must not use non-exact derivative substitutes.

Allowed derivative mechanisms:

- exact analytic derivatives
- CppAD for explicit algebraic derivatives
- analytic implicit sensitivities
- CppAD implicit sensitivities

Production solvers and regression loops must be native C++ where the package owns the algorithm.

Python may:

- validate inputs
- construct problem objects
- serialize data
- call native code
- format outputs

Python may not own production optimizer loops, residual packing, or solver iterations for required workflows.

---

# 5. What prior agents got wrong

Prior agents often completed a narrow safe slice and marked the issue done.

That is no longer allowed.

Common failure modes:

- inventory treated as benchmark implementation
- diagnostic route treated as solver implementation
- staged workflow treated as coupled equilibrium
- synthetic payload treated as downstream integration
- missing paths recorded as limitations
- unsupported parameter families called out of scope
- Ceres dependency presence treated as native Ceres production regression
- CppAD availability treated as complete derivative coverage
- association sensitivity deferred without implementing the actual required path
- GoalBuddy final audit marked `full_outcome_complete: true` despite known missing required workflows

A PR may not close an issue by being honest about incompleteness. Honesty is required, but incompleteness means the issue remains open.

---

# 6. Completion standard

An issue is complete only when every in-scope workflow is:

1. Implemented in production code.
2. Exercised through public generic APIs.
3. Backed by native C++ implementation where the package owns the algorithm.
4. Using analytic, CppAD, analytic-implicit, or CppAD-implicit derivatives.
5. Tested with executable tests.
6. Validated with real or literature-backed data when relevant.
7. Free of banned literal tokens in committed text.
8. Free of synthetic-only or monkeypatched proof.
9. Reflected in capabilities without overclaiming.
10. Usable by downstream projects without private package workarounds.

The final audit question is:

```text
Can a downstream project use this generic package workflow now without adding missing package code?
```

If the answer is no, the issue is not complete.

---

# 7. Required architecture

## EOS harness

Must support:

- residual Helmholtz energy and contribution maps
- property evaluation
- fugacity coefficients
- activity coefficients
- chemical potentials
- relative permittivity
- density and pressure closure
- parameter derivatives
- composition derivatives
- temperature/density derivatives where required
- contribution-level traceability

## Derivatives

Must support:

- CppAD scalar substrate across EOS and residual layers
- association solved-state implicit sensitivities
- density-root implicit sensitivities
- bubble/dew solved-state sensitivities
- LLE residual Jacobians
- electrolyte LLE residual Jacobians
- reactive/speciation residual Jacobians
- reactive LLE coupled residual Jacobians
- regression objective Jacobians

No direct differentiation through iterative loops is accepted as production derivative support.

## Association

Association site fractions are solved internal states.

Required:

- stable association value evaluation
- residual equations for association site fractions
- Jacobian of association residuals with respect to site fractions
- parameter/composition/density sensitivities through implicit differentiation
- support for associating pure components and associating mixtures
- support for `k_hb_ij` when association interactions require it
- support for `k_ij` regression even when one or both components associate

The package must not claim `k_ij` or `k_hb_ij` cannot be supported merely because association is involved. The required solution is association implicit sensitivity plus native regression integration.

## Regression

The package must provide native Ceres or equivalent native optimizer loops for:

- pure neutral regression
- pure ion regression
- binary `k_ij`
- binary `l_ij`
- binary `k_hb_ij`
- association parameters
- `d_born`
- `f_solv`
- relative-permittivity parameters
- density rows
- vapor-pressure rows
- fugacity rows
- activity rows
- MIAC rows
- osmotic rows
- VLE partial-pressure rows
- LLE phase-composition rows
- reactive speciation rows
- reactive pressure/speciation rows

Required regression result fields:

- optimizer backend
- derivative backend
- parameter map
- initial parameters
- final parameters
- parameter movement
- active bounds
- objective initial
- objective final
- row diagnostics
- source summaries
- target-family summaries
- residual block norms
- convergence status

Current claim boundary:

- Reactive regression is currently exposed as structured residual evaluation through `build_reactive_regression_objective(...)` and `evaluate_reactive_regression_objective(...)`. Do not describe that surface as a production optimizer until native optimizer ownership and capability evidence say so.

## Equilibrium

Required production solvers:

- neutral VLE / bubble / dew
- direct fugacity-based volatile partial pressure where valid
- neutral LLE
- electrolyte LLE with distributed ions
- activity-based speciation
- reactive VLE / bubble
- reactive LLE
- reactive electrolyte LLE

Staged workflows may exist only as initialization or diagnostics. They do not satisfy coupled equilibrium requirements.

## Benchmarks

A literature benchmark is complete only when it is executable and tolerance-checked.

A manifest, inventory, or data list does not count.

Required benchmark families:

- Gross/Sadowski pure PC-SAFT parameters
- Gross/Sadowski associating systems
- Baygi MEA PC-SAFT reproduction
- Held/Cameretti aqueous electrolyte density/MIAC
- Held alcohol/salt MIAC and mixed-solvent prediction
- Held 2014 revised ePC-SAFT electrolyte phase behavior
- Bülow/Ascani/Held advanced Born and dielectric behavior
- Figiel 2025 modified Born / SSM / DS
- Ascani 2022 electrolyte LLE with distributed ions
- Ascani 2023 reactive phase equilibrium
- Khudaida 2026 electrolyte salting-out LLE
- Rezaee lithium extraction thermodynamic modeling
- MEA-CO2 pressure/speciation data validation with copied fitted parameters
- downstream-facing lithium/column contract workflows with in-worktree inputs

---

# 8. Required milestones

## Milestone 0 — Process hard gates

Complete before more feature work.

- purge banned literal tokens from the repo
- add lexical guard
- add GoalBuddy lifecycle guard
- remove classification-as-completion language
- update future prompts
- ensure no merged goal remains active

## Milestone 1 — Native derivative substrate

- package-wide CppAD scalar support
- scalar-templated EOS residuals
- association implicit sensitivities
- density-root implicit sensitivities
- bubble/dew sensitivities
- property derivative result APIs
- derivative coverage tests

## Milestone 2 — Native regression backend

- Ceres owns production optimizer loops
- CppAD/implicit derivatives feed Jacobians
- Python only validates/serializes
- all required target families supported
- all binary parameter families supported
- real parameter movement tests

## Milestone 3 — Production equilibrium backend

- native coupled activity speciation
- native derivative-backed neutral LLE
- native electrolyte LLE with distributed ions
- native coupled reactive LLE
- native coupled reactive electrolyte LLE
- real residual/Jacobian diagnostics

## Milestone 4 — Literature benchmarks

- every named benchmark executable
- every benchmark has fixture, command, expected result, tolerance
- benchmark suite proves parameter/equilibrium behavior

## Milestone 5 — Real downstream integration

- install local package into each downstream repository
- run one real workflow per downstream repo
- prove generic APIs are sufficient
- no copied EOS implementation
- no private workaround for required package behavior

---

# 9. Agent working rules

Every agent must:

1. Read this file.
2. Read the replacement master issue.
3. Inspect current `origin/main`.
4. State whether its task is derivative, regression, equilibrium, benchmark, process, or downstream.
5. Identify the exact production workflow it must make work.
6. Identify tests that prove it.
7. Refuse to close by inventory, diagnostics, staging, synthetic fixtures, or limitation text.
8. Avoid banned exact literals in committed text.
9. Keep APIs generic.
10. Open a PR only when the production workflow works.

Every final PR must answer:

```text
What downstream workflow can now run that could not run before?
What production native code path now exists?
What derivative path is used?
What real data or benchmark proves it?
What public API exercises it?
What tests would fail if this regressed?
```

If those answers are weak, do not merge.
