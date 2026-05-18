# ePC-SAFT Package Context

This is the domain glossary for the `epcsaft` upstream package. It names the package-level thermodynamics, solver, regression, benchmark, and downstream-integration concepts that architecture, issue, and debugging skills should use consistently.

## Language

### Package Identity

**ePC-SAFT Package**:
The upstream thermodynamic package that owns generic ePC-SAFT model interfaces, native implementation paths, tests, documentation, and release behavior.
_Avoid_: downstream study, application wrapper, case-study package

**Downstream Consumer**:
A separate project that uses package outputs for application-specific studies, process models, figures, reports, or manuscript claims.
_Avoid_: package feature, built-in application, private workaround

**Generic Package Workflow**:
A reusable package workflow exposed through public package interfaces and valid beyond one downstream case study.
_Avoid_: custom script, repo-specific workaround, synthetic demonstration

**Capability Contract**:
The package's current, validated statement of what workflows are supported, including the production path and evidence behind each claim.
_Avoid_: aspiration, roadmap item, dependency presence

### Thermodynamic Model

**EOS Harness**:
The package module family that evaluates ePC-SAFT states, residual Helmholtz contributions, closures, properties, activities, fugacities, chemical potentials, and derivatives.
_Avoid_: property helper collection, equation utilities

**Residual Helmholtz Model**:
The coupled thermodynamic model whose residual Helmholtz energy is assembled from hard-chain, dispersion, association, ionic, Born, and related electrolyte contribution families.
_Avoid_: independent property calculators, isolated terms

**Contribution Family**:
A named thermodynamic term family in the residual Helmholtz model with traceable equations, implementation ownership, and derivative behavior.
_Avoid_: miscellaneous correction, helper formula

**Parameter Family**:
A supported class of fitted or supplied model parameters, such as pure-component parameters, binary interaction parameters, association parameters, ionic parameters, or Born-related parameters.
_Avoid_: optional vocabulary, ad hoc coefficient

**State**:
A thermodynamic state specification with resolved temperature, pressure or density closure, composition, phase assumptions, and parameter context.
_Avoid_: input row, loose condition dictionary

### Equilibrium And Regression

**Equilibrium Problem**:
A generic package problem that asks the solver to satisfy phase, chemical, reactive, or electrolyte equilibrium conditions for a declared model system.
_Avoid_: process calculation, downstream metric

**Electrolyte LLE Problem**:
An equilibrium problem for liquid-liquid phase split calculations with distributed ions and phase electroneutrality constraints.
_Avoid_: extraction efficiency calculator, brine-specific script

**Reactive LLE Problem**:
An equilibrium problem that couples liquid-liquid phase equilibrium with chemical or speciation equilibrium.
_Avoid_: staged speciation then split workflow

**Regression Problem**:
A generic package problem that fits supported parameter families to target datasets through a native optimizer loop and derivative-backed residuals.
_Avoid_: curve-fit script, downstream calibration notebook

**Target Dataset**:
A structured set of target rows with source context, quantities, units, uncertainty or tolerance expectations, and target-family classification.
_Avoid_: loose CSV, anonymous fixture

**Production Solver Path**:
The validated native implementation route that owns residual packing, solver iteration, diagnostics, and result contracts for a required workflow.
_Avoid_: diagnostic route, staged-only route, Python-owned optimizer loop

**Derivative Path**:
The validated analytic, CppAD, analytic-implicit, or CppAD-implicit derivative route used by a solver or regression workflow.
_Avoid_: numerical derivative fallback, unspecified derivative support

### Evidence And Completion

**Literature Benchmark**:
An executable, tolerance-checked reproduction of a named literature system or data family using package-owned generic interfaces.
_Avoid_: bibliography entry, manifest, inventory

**Validation Lane**:
An explicit command, fixture set, or workflow that proves a package contract at a defined confidence level.
_Avoid_: hand inspection, partial transcript, stale note

**Completion Evidence**:
The production code, public-interface exercise, native path, derivative path, real data or benchmark, and regression test evidence needed before a workflow can be called complete.
_Avoid_: documented limitation, honest incompleteness, dependency-only proof

## Relationships

- The **ePC-SAFT Package** serves one or more **Downstream Consumers** through **Generic Package Workflows**.
- A **Capability Contract** is valid only when backed by **Completion Evidence**.
- The **EOS Harness** implements the **Residual Helmholtz Model** through traceable **Contribution Families**.
- A **State** combines model parameters, composition, phase assumptions, and closure information for the **EOS Harness**.
- An **Equilibrium Problem** uses a **Production Solver Path** and a **Derivative Path** to produce generic thermodynamic outputs.
- An **Electrolyte LLE Problem** and a **Reactive LLE Problem** are specialized **Equilibrium Problems**.
- A **Regression Problem** fits **Parameter Families** to a **Target Dataset** using a **Production Solver Path** and a **Derivative Path**.
- A **Literature Benchmark** is a high-confidence **Validation Lane** for package behavior.

## Example Dialogue

> **Dev:** "Can the lithium project use the package to calculate selectivity directly?"
> **Domain expert:** "No. The **ePC-SAFT Package** should solve the **Electrolyte LLE Problem** and return generic phase outputs; the **Downstream Consumer** computes selectivity from those outputs."

> **Dev:** "Does a staged speciation script prove the reactive split?"
> **Domain expert:** "No. A **Reactive LLE Problem** needs a coupled **Production Solver Path** plus a validated **Derivative Path** before the **Capability Contract** can claim support."

## Flagged Ambiguities

- "Package workflow" has sometimes meant a source-checkout analysis script; resolved: use **Generic Package Workflow** only for public, reusable package behavior.
- "Benchmark" has sometimes meant an inventory or manifest; resolved: use **Literature Benchmark** only for executable, tolerance-checked evidence.
- "Staged workflow" has sometimes been treated as coupled equilibrium; resolved: staged workflows may initialize or diagnose but do not satisfy a **Reactive LLE Problem** contract.
- "Downstream integration" has sometimes implied adding application metrics to the package; resolved: **Downstream Consumers** own application metrics, while the package owns generic thermodynamic outputs.
