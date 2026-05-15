# Modified Born, Electrolyte Regression, And Figiel-Held Benchmark Goal

Source authority: GitHub issue #118, "Complete modified Born / SSM / DS electrolyte regression and Figiel-Held benchmark proof"

Issue URL: https://github.com/tannerpolley/ePC-SAFT/issues/118

## Outcome

Execute issue #118 as a vertical implementation goal. The package must support the production chain from modified Born / SSM / DS liquid-electrolyte property derivatives, through ion/solvent parameter derivative paths, through native Ceres electrolyte regression, to Figiel/Held literature-backed benchmark proof.

This board setup does not complete issue #118. It creates the execution board so a later `/goal` run can proceed through the issue phases in order.

## Hard Constraints

- Follow issue #118 exactly. If this file and the issue disagree, the GitHub issue is authoritative.
- Do not start source edits until `docs/goals/modified-born-electrolyte-regression/notes/intake.md` is fully populated with the Stage 0 evidence required by issue #118.
- Do not close, mark ready, or call issue #118 complete from inventories, manifests, schema-only support, diagnostics-only routes, staged-only helpers, synthetic-only payloads, mocked solver output, documented limitations, capability labels, or assertions without exact passing tests.
- Use only production derivative mechanisms allowed by issue #118: analytic, CppAD, analytic implicit sensitivity, or CppAD implicit sensitivity.
- For solved internal states, use implicit sensitivity through the local residual system, including density roots, association site fractions, reaction/speciation variables, bubble/dew roots, LLE phase-split variables, and reactive phase-equilibrium variables.
- Ceres must own the production optimizer loop for electrolyte regression. Python may validate inputs, call native code, and format outputs, but must not own the production optimizer or nonlinear-solver loops for required workflows.
- Keep `epcsaft` a general-purpose package. Do not add public APIs named after downstream applications or downstream metrics.
- Do not edit phase-equilibrium solvers unless the change is only an import or path compatibility fix.
- The production claim is liquid-electrolyte only unless vapor electrolyte support is implemented and tested inside this issue.
- Do not commit the old forbidden backend token, old forbidden numeric derivative token, or old forbidden numeric derivative phrase as contiguous text. Construct guard searches from fragments when needed.
- The actual modified Born / SSM / DS formula must be taken from the repository equation registry and Figiel/Bulow/Held source equations, not from memory or from the simplified reminder in the issue.

## Required Phase Order

1. Stage 0 - Intake gate: create and fill `notes/intake.md` with current Born / SSM / DS contribution files, current derivative tests, current Ceres electrolyte regression tests, current Figiel/Held fixture paths, and the chosen final benchmark set.
2. Stage 1 - Modified Born / SSM / DS derivatives: implement and test liquid-electrolyte derivatives for `d_born`, `f_solv`, relative-permittivity parameters, ion diameter or solvated diameter, ion-solvent dispersion energy, composition, density, and temperature where required.
3. Stage 2 - Native Ceres electrolyte regression: make native Ceres fit the required electrolyte parameter families against density, osmotic coefficient, mean ionic activity coefficient, relative permittivity, and solvation or transfer Gibbs-energy rows where fixtures exist.
4. Stage 3 - Benchmark proof: run Figiel 2025 modified Born / SSM / DS, Held/Cameretti aqueous electrolyte density and MIAC, and Held alcohol/salt mixed-solvent density/osmotic/MIAC benchmark families with numeric tolerances.
5. Required validation - run the exact build, pytest, docs, quick-validation, and whitespace checks named in issue #118.
6. Final audit - verify every issue #118 completion line is true before calling the issue complete; otherwise leave the issue open or the PR draft with a stopped-state report.

## Benchmark Report Requirements

The final benchmark report must include:

- source
- target rows
- initial parameters
- final parameters
- objective initial and final
- parameter movement
- active bounds
- derivative backend
- tolerances

## Completion Proof

The goal is complete only when every issue #118 completion line is true: liquid-electrolyte modified Born / SSM / DS derivatives are validated, native Ceres fits electrolyte parameters with production derivatives, the Figiel/Held benchmark suite runs with numeric tolerances, parameter movement and objective decrease are reported, docs state liquid-electrolyte scope exactly, and no vapor electrolyte Born claim is made without implementation and tests.

The final audit must answer:

- What production native code path now exists?
- What derivative path is used?
- What public generic API exercises it?
- What real data or literature-backed benchmark proves it?
- What tests would fail if this regressed?
- Can downstream projects use the generic package workflow without private package workarounds?

## Starter Command

`/goal Follow docs/goals/modified-born-electrolyte-regression/goal.md.`
