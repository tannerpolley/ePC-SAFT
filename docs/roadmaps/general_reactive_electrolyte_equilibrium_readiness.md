# General Reactive/Electrolyte Equilibrium Readiness Roadmap

## Purpose

This roadmap defines the next staged development layer for `epcsaft`.

The package must remain a generalized, open-source ePC-SAFT package. It should be usable by many downstream engineering projects, including but not limited to MEA thermodynamics, absorption-column modeling, and lithium brine extraction studies.

The public package must not be application-specific.

Forbidden public API examples:

```python
fit_lithium_extraction_parameters(...)
screen_lithium_extractants(...)
fit_mea_absorption(...)
fit_co2_capture_column(...)
```

Preferred public API concepts:

```python
model.equilibrium(problem)
epcsaft.equilibrium(model, problem)
epcsaft.regress_parameters(problem)
ReactionSet(...)
EquilibriumProblem(...)
RegressionProblem(...)
TargetDataset(...)
PhaseSpec(...)
ParameterSet(...)
```

Downstream projects compute application-specific metrics from generic package outputs.

## Package responsibilities

The package owns:

```text
EOS evaluation
residual Helmholtz contributions
activity/fugacity/chemical-potential evaluation
density/pressure closures
phase equilibrium
reaction/speciation equilibrium
electrolyte equilibrium
generic target-row regression
native optimizer integration
derivative backends
capability reporting
literature benchmarks
```

Downstream projects own:

```text
process-specific data curation
application-specific metrics
figures and reports
surrogate models
column models
extraction-efficiency calculations
distribution coefficients
selectivity calculations
```

## Derivative policy

```text
No finite difference.
CppAD for explicit algebraic derivatives.
Exact analytic derivatives may remain where validated.
Solved internal states use analytic_implicit or cppad_implicit sensitivities.
backend_unavailable is allowed only for out-of-scope workflows.
backend_unavailable is a blocker for required workflows.
```

Solved internal states include:

```text
association site fractions
density roots
speciation solves
VLE/bubble/dew roots
LLE phase splits
reactive LLE reaction coordinates
stability analysis solves
```

For solved states:

```text
R(u, theta) = 0
u_theta = - R_u^{-1} R_theta
```

Do not tape iterative solver loops as the production derivative.

## Phase 1 — General EOS, property, and simple regression readiness

Focus:

```text
pure parameter fitting
binary parameter fitting
direct property evaluation
activity/fugacity/chemical-potential outputs
ideal/apparent speciation workflow
MEA Smith-Missen-style benchmark
basic non-electrolyte LLE benchmark
```

## Phase 2 — Activity-based equilibrium readiness

Focus:

```text
reaction-constant convention layer
activity-based speciation
VLE fugacity equilibrium
electrolyte LLE with distributed ions
MDEA ePC-SAFT benchmark
Ascani 2022 electrolyte LLE benchmark
```

## Phase 3 — Chemical + phase equilibrium readiness

Focus:

```text
reactive LLE
chemical + phase equilibrium
implicit sensitivities for solved states
Ascani 2023 reactive LLE benchmark
MEA activity-based speciation/VLE benchmark
```

## Phase 4 — General regression and downstream integration readiness

Focus:

```text
generic target-row regression
native optimizer with real Jacobians
optional reaction-constant fitting
downstream smoke tests
MEA-Thermodynamics integration
Lithium_Extraction integration
MEA-Absorption-Column integration
```

## LLE progression

The solver progression should be:

```text
plain non-electrolyte LLE
electrolyte LLE with distributed ions
reactive LLE / chemical phase equilibrium
```

This keeps the phase split logic modular while acknowledging that electrolyte accounting adds enough complexity to deserve a separate issue.

## Literature anchors

Use paper-specific systems as benchmarks, not as public API names.

Relevant benchmark groups:

```text
PC-SAFT pure/associating systems
MEA simple Smith-Missen style workflow
MDEA activity-based ePC-SAFT workflow
CO2 electrolyte solubility / carbonate benchmark
Figiel 2025 SSM+DS Born benchmark
Held 2014 ePC-SAFT revised benchmark
Bülow/Ascani dielectric/Born benchmarks
non-electrolyte LLE benchmark
Ascani 2022 electrolyte LLE benchmark
Ascani 2023 reactive LLE benchmark
Khudaida electrolyte salting-out LLE benchmark
Hubach/Yu lithium-related electrolyte LLE benchmark
Rezaee DES parameterization benchmark
```

## Issue sequence

```text
A. Derivative backend completion audit and coverage matrix hard gate
B. Explicit CppAD parameter derivatives for EOS/property APIs
C. Generic implicit sensitivity framework for solved states
D. General reaction/equilibrium-constant convention layer
E. Generic target-row and dataset schema
F. Generic speciation solver using ePC-SAFT activities
G. Generic VLE/fugacity-equilibrium solver for volatile neutral species
H. Generic non-electrolyte LLE benchmark and solver hardening
I. Generic electrolyte LLE with distributed ions
J. Generic reactive LLE and chemical phase equilibrium
K. Generic regression row schema and native optimizer backend
L. Literature benchmark suite
M. Downstream integration smoke tests
```

## API preference

Concise but typed:

```python
result = model.equilibrium(problem)
fit = epcsaft.regress_parameters(problem)
```

The `problem` object determines route and diagnostics.

Examples:

```python
LLE(...)
ElectrolyteLLE(...)
ReactiveEquilibrium(...)
ReactiveLLE(...)
BubblePoint(...)
TPFlash(...)
RegressionProblem(...)
```

## Final acceptance standard

The package is ready for the next downstream stage only when:

```text
No finite difference remains.
Derivative coverage is explicit and truthful.
Required paths are not backend_unavailable.
Generic equilibrium APIs are stable.
Generic regression APIs are stable.
Literature benchmarks prove capabilities.
Downstream projects can use generic APIs without copied EOS code or application-specific package methods.
```
