# Algorithm Index

This file is generated from `docs/latex/algorithms.tex` by `scripts/docs/sync_algorithm_registry.py`.
The LaTeX document remains the current source of truth; this Markdown view and `docs/algorithms_registry.yaml` stay aligned with it.

## Core Optimizer and Adapter Algorithms

### `native_nlp_problem_contract`
- Family: Core optimizer adapter
- Status: Implemented
- Public API: Private native contract used by package equilibrium routes
- Backend: Native C++
- Dependency: None
- Derivative backend: Route-provided exact gradient and Jacobian callbacks
- Solver role: Shared interface for Ipopt-owned constrained NLP routes
- Implementation owner: src/epcsaft/native/equilibrium_nlp/nlp_problem.h
- Validation: tests/native/equilibrium/test_ipopt_adapter_contract.py
- Capability key: native_nlp_problem_contract
- Description: Defines the common native NLP shape consumed by the Ipopt adapter.
- Change note: Initial algorithm-registry entry for the shared native NLP contract.
- LaTeX: `docs/latex/algorithms.tex:17`
- Code owners: `src/epcsaft/native/equilibrium_nlp/nlp_problem.h:33` (class NlpProblem {)

The `NlpProblem` contract is the common native interface for objective,
constraint, Jacobian, bounds, scaling, and diagnostics payloads used by
Ipopt-backed equilibrium and speciation routes.

**LaTeX source**

```tex
The \texttt{NlpProblem} contract is the common native interface for objective,
constraint, Jacobian, bounds, scaling, and diagnostics payloads used by
Ipopt-backed equilibrium and speciation routes.
```


### `ipopt_tnlp_adapter`
- Family: Core optimizer adapter
- Status: Implemented
- Public API: Private native adapter reached through package equilibrium routes
- Backend: Native C++ Ipopt TNLP adapter
- Dependency: Ipopt
- Derivative backend: Route-provided exact gradient and Jacobian callbacks
- Solver role: Maps NlpProblem into Ipopt::TNLP and returns IpoptSolveResult
- Implementation owner: src/epcsaft/native/equilibrium_nlp/ipopt_adapter.cpp
- Validation: tests/native/equilibrium/test_ipopt_adapter_contract.py
- Capability key: native_ipopt_equilibrium_nlp
- Description: Bridges the package NLP contract to Ipopt without owning route equations.
- Change note: Initial algorithm-registry entry for native Ipopt dispatch.
- LaTeX: `docs/latex/algorithms.tex:35`
- Code owners: `src/epcsaft/native/equilibrium_nlp/ipopt_adapter.cpp:325` (class IpoptTnlpAdapter final : public Ipopt::TNLP {)

The Ipopt adapter owns solver invocation, option transfer, status mapping, and
diagnostic capture for constrained NLP routes. Route-specific residuals and
Jacobians remain in their route problem owners.

**LaTeX source**

```tex
The Ipopt adapter owns solver invocation, option transfer, status mapping, and
diagnostic capture for constrained NLP routes. Route-specific residuals and
Jacobians remain in their route problem owners.
```


### `native_ceres_regression_adapter`
- Family: Core optimizer adapter
- Status: Implemented
- Public API: Private native adapter reached through regression helpers
- Backend: Native C++ Ceres least-squares adapter
- Dependency: Ceres
- Derivative backend: Native analytic or CppAD-backed residual Jacobians, route dependent
- Solver role: Maps regression residual/Jacobian evaluators into ceres::CostFunction and ceres::Problem
- Implementation owner: src/epcsaft/native/epcsaft_regression.cpp
- Validation: tests/regression/literature/test_literature_binary_kij_regression.py
- Capability key: native_ceres_regression
- Description: Shared adapter pattern for currently supported native Ceres regression routes.
- Change note: Initial algorithm-registry entry for native Ceres regression dispatch.
- LaTeX: `docs/latex/algorithms.tex:53`
- Code owners: `src/epcsaft/bindings.cpp:3601` (m.def("_fit_generic_native_ceres", &fit_generic_native_ceres_binding);), `src/epcsaft/native/epcsaft_regression.cpp:565` (ceres::Solver::Options ceres_regression_options_cpp(int max_num_iterations) {)

The Ceres adapter family owns native least-squares solve setup for supported
regression routes. It is separate from Ipopt and does not call Ipopt routes.

**LaTeX source**

```tex
The Ceres adapter family owns native least-squares solve setup for supported
regression routes. It is separate from Ipopt and does not call Ipopt routes.
```


## Equilibrium Algorithms

### `neutral_tp_flash_ipopt`
- Family: Phase equilibrium
- Status: Implemented
- Public API: mix.equilibrium(kind="tp_flash", ...); mix.solve_equilibrium(TPFlash(...))
- Backend: Native C++ Ipopt equilibrium NLP
- Dependency: Ipopt
- Derivative backend: Exact native objective gradient and constraint Jacobian
- Solver role: Native neutral TP flash two-phase EOS constrained NLP
- Implementation owner: src/epcsaft/equilibrium.py; src/epcsaft/native/equilibrium_nlp/route_builders.cpp; src/epcsaft/bindings.cpp
- Validation: tests/equilibrium/core/test_tp_flash.py
- Capability key: neutral_tp_flash
- Description: Solves neutral two-phase TP flash through the native Ipopt EOS route.
- Change note: Initial algorithm-registry entry for neutral TP flash.
- LaTeX: `docs/latex/algorithms.tex:72`
- Code owners: `src/epcsaft/bindings.cpp:2208` (m.def("_native_neutral_tp_flash_eos_nlp_contract", [](), `src/epcsaft/bindings.cpp:2879` (m.def("_native_neutral_tp_flash_eos_route_result", [](), `src/epcsaft/equilibrium.py:297` (@dataclass(frozen=True, slots=True)), `src/epcsaft/native/equilibrium_nlp/route_builders.cpp:1831` (NeutralTwoPhaseEosNlpContract evaluate_neutral_two_phase_eos_tp_flash_nlp_contract(), `src/epcsaft/native/equilibrium_nlp/route_builders.cpp:2720` (NeutralTwoPhaseEosRouteResult solve_neutral_two_phase_eos_tp_flash_route()

This entry covers the generic neutral TP flash route, including the typed
`TPFlash` problem object and the string facade that adapts to it.

**LaTeX source**

```tex
This entry covers the generic neutral TP flash route, including the typed
\texttt{TPFlash} problem object and the string facade that adapts to it.
```


### `neutral_lle_ipopt`
- Family: Phase equilibrium
- Status: Implemented
- Public API: mix.equilibrium(kind="lle_flash", ...); mix.solve_equilibrium(LLEProblem(...))
- Backend: Native C++ Ipopt equilibrium NLP
- Dependency: Ipopt
- Derivative backend: Exact native objective gradient and constraint Jacobian
- Solver role: Native neutral liquid-liquid two-phase EOS constrained NLP
- Implementation owner: src/epcsaft/equilibrium.py; src/epcsaft/native/equilibrium_nlp/route_builders.cpp; src/epcsaft/bindings.cpp
- Validation: tests/equilibrium/core/test_lle.py
- Capability key: neutral_lle_flash
- Description: Solves neutral LLE through the native Ipopt EOS route.
- Change note: Initial algorithm-registry entry for neutral LLE.
- LaTeX: `docs/latex/algorithms.tex:89`
- Code owners: `src/epcsaft/bindings.cpp:2227` (m.def("_native_neutral_lle_eos_nlp_contract", [](), `src/epcsaft/bindings.cpp:2924` (m.def("_native_neutral_lle_eos_route_result", [](), `src/epcsaft/equilibrium.py:410` (@dataclass(frozen=True, slots=True)), `src/epcsaft/native/equilibrium_nlp/route_builders.cpp:1860` (NeutralTwoPhaseEosNlpContract evaluate_neutral_two_phase_eos_lle_nlp_contract(), `src/epcsaft/native/equilibrium_nlp/route_builders.cpp:2749` (NeutralTwoPhaseEosRouteResult solve_neutral_two_phase_eos_lle_route()

This entry covers the generic neutral liquid-liquid split route and its typed
`LLEProblem` public problem object.

**LaTeX source**

```tex
This entry covers the generic neutral liquid-liquid split route and its typed
\texttt{LLEProblem} public problem object.
```


### `electrolyte_lle_ipopt`
- Family: Electrolyte phase equilibrium
- Status: Implemented
- Public API: mix.equilibrium(kind="electrolyte_lle", ...); mix.solve_equilibrium(ElectrolyteLLEProblem(...))
- Backend: Native C++ Ipopt equilibrium NLP
- Dependency: Ipopt
- Derivative backend: Exact native objective gradient and constraint Jacobian
- Solver role: Native charge-constrained electrolyte LLE EOS constrained NLP
- Implementation owner: src/epcsaft/equilibrium.py; src/epcsaft/native/equilibrium_nlp/route_builders.cpp; src/epcsaft/bindings.cpp
- Validation: tests/equilibrium/electrolyte/test_electrolyte_lle_smokes.py
- Capability key: electrolyte_lle
- Description: Solves electrolyte liquid-liquid split routes with charge-constrained phase amounts.
- Change note: Initial algorithm-registry entry for electrolyte LLE.
- LaTeX: `docs/latex/algorithms.tex:106`
- Code owners: `src/epcsaft/bindings.cpp:2246` (m.def("_native_electrolyte_lle_eos_nlp_contract", [](), `src/epcsaft/bindings.cpp:2969` (m.def("_native_electrolyte_lle_eos_route_result", [](), `src/epcsaft/equilibrium.py:433` (@dataclass(frozen=True, slots=True)), `src/epcsaft/native/equilibrium_nlp/route_builders.cpp:1889` (NeutralTwoPhaseEosNlpContract evaluate_electrolyte_lle_eos_nlp_contract(), `src/epcsaft/native/equilibrium_nlp/route_builders.cpp:2778` (NeutralTwoPhaseEosRouteResult solve_electrolyte_lle_eos_route()

This entry covers the generic electrolyte LLE route. Downstream extraction or
selectivity metrics remain outside the package algorithm.

**LaTeX source**

```tex
This entry covers the generic electrolyte LLE route. Downstream extraction or
selectivity metrics remain outside the package algorithm.
```


### `bubble_dew_ipopt`
- Family: Phase equilibrium
- Status: Implemented
- Public API: mix.equilibrium(kind="bubble_p"|"bubble_t"|"dew_p"|"dew_t", ...); mix.solve_equilibrium(BubblePoint(...) or DewPoint(...))
- Backend: Native C++ Ipopt equilibrium NLP
- Dependency: Ipopt
- Derivative backend: Exact native objective gradient and constraint Jacobian
- Solver role: Native fixed-composition bubble/dew constrained NLP routes
- Implementation owner: src/epcsaft/equilibrium.py; src/epcsaft/native/equilibrium_nlp/bubble_dew_route_builders.cpp; src/epcsaft/bindings.cpp
- Validation: tests/equilibrium/core/test_bubble_dew.py
- Capability key: neutral_bubble_dew
- Description: Solves neutral bubble and dew route variants through native Ipopt.
- Change note: Initial algorithm-registry entry for neutral bubble/dew routes.
- LaTeX: `docs/latex/algorithms.tex:123`
- Code owners: `src/epcsaft/bindings.cpp:2457` (m.def("_native_neutral_bubble_p_eos_nlp_contract", [](), `src/epcsaft/bindings.cpp:3019` (m.def("_native_neutral_bubble_p_eos_route_result", [](), `src/epcsaft/equilibrium.py:358` (@dataclass(frozen=True, slots=True)), `src/epcsaft/equilibrium.py:384` (@dataclass(frozen=True, slots=True)), `src/epcsaft/native/equilibrium_nlp/bubble_dew_route_builders.cpp:1863` (NeutralTwoPhaseEosRouteResult solve_pressure_route(), `src/epcsaft/native/equilibrium_nlp/bubble_dew_route_builders.cpp:1999` (NeutralTwoPhaseEosRouteResult solve_temperature_route()

This entry groups the neutral bubble-pressure, bubble-temperature,
dew-pressure, and dew-temperature routes because they share the same native
route-builder family.

**LaTeX source**

```tex
This entry groups the neutral bubble-pressure, bubble-temperature,
dew-pressure, and dew-temperature routes because they share the same native
route-builder family.
```


### `stability_tpd_ipopt`
- Family: Stability analysis
- Status: Implemented
- Public API: mix.equilibrium(kind="stability", ...); mix.equilibrium(kind="electrolyte_stability", ...)
- Backend: Native C++ Ipopt stability NLP
- Dependency: Ipopt
- Derivative backend: Exact native objective gradient and constraint Jacobian
- Solver role: Native tangent-plane-distance stability route
- Implementation owner: src/epcsaft/equilibrium.py; src/epcsaft/native/equilibrium_nlp/stability_route_builders.cpp; src/epcsaft/bindings.cpp
- Validation: tests/equilibrium/core/test_stability.py
- Capability key: neutral_stability; electrolyte_stability
- Description: Solves neutral and electrolyte TPD trial routes through native Ipopt.
- Change note: Initial algorithm-registry entry for stability routes.
- LaTeX: `docs/latex/algorithms.tex:141`
- Code owners: `src/epcsaft/bindings.cpp:2538` (m.def("_native_neutral_stability_tpd_nlp_contract", [](), `src/epcsaft/bindings.cpp:3225` (m.def("_native_neutral_stability_tpd_route_result", [](), `src/epcsaft/equilibrium.py:315` (@dataclass(frozen=True, slots=True)), `src/epcsaft/native/equilibrium_nlp/stability_route_builders.cpp:1016` (StabilityRouteResult solve_neutral_stability_tpd_route(), `src/epcsaft/native/equilibrium_nlp/stability_route_builders.cpp:1090` (StabilityRouteResult solve_electrolyte_stability_tpd_route()

This entry covers neutral and electrolyte tangent-plane-distance trial routes.
It does not describe downstream phase-selection policy beyond the returned
generic stability payload.

**LaTeX source**

```tex
This entry covers neutral and electrolyte tangent-plane-distance trial routes.
It does not describe downstream phase-selection policy beyond the returned
generic stability payload.
```


## Reactive and Speciation Algorithms

### `ideal_speciation_ipopt`
- Family: Chemical equilibrium
- Status: Implemented
- Public API: mix.chemical_equilibrium(...); ReactiveSpeciationProblem
- Backend: Native C++ Ipopt ideal Gibbs speciation NLP
- Dependency: Ipopt
- Derivative backend: Analytic or CppAD ideal log-amount residual Jacobian
- Solver role: Native homogeneous ideal reactive speciation route
- Implementation owner: src/epcsaft/reactive_speciation.py; src/epcsaft/native/equilibrium_nlp/ideal_speciation_problem.cpp; src/epcsaft/native/epcsaft_chemical_equilibrium.cpp; src/epcsaft/bindings.cpp
- Validation: tests/api/reactive/test_reactive_speciation.py
- Capability key: reactive_speciation
- Description: Solves ideal homogeneous reactive speciation through native Ipopt.
- Change note: Initial algorithm-registry entry for ideal speciation.
- LaTeX: `docs/latex/algorithms.tex:161`
- Code owners: `src/epcsaft/bindings.cpp:3607` (m.def("_solve_chemical_equilibrium_native", &solve_chemical_equilibrium_native_binding);), `src/epcsaft/equilibrium.py:489` (@dataclass(frozen=True, slots=True)), `src/epcsaft/native/equilibrium_nlp/ideal_speciation_problem.cpp:699` (IdealSpeciationIpoptResult solve_ideal_speciation_ipopt(), `src/epcsaft/native/equilibrium_nlp/ideal_speciation_problem.cpp:735` (ChemicalEquilibriumResultNative solve_ideal_speciation_chemical_equilibrium_ipopt(), `src/epcsaft/reactive_speciation.py:625` (def _solve_reactive_speciation_native()

This entry covers the ideal-standard-state speciation route. It is a chemical
equilibrium algorithm, not a phase-equilibrium split.

**LaTeX source**

```tex
This entry covers the ideal-standard-state speciation route. It is a chemical
equilibrium algorithm, not a phase-equilibrium split.
```


### `nonideal_speciation_ipopt`
- Family: Chemical equilibrium
- Status: Implemented
- Public API: mix.chemical_equilibrium(...); ReactiveSpeciationProblem
- Backend: Native C++ Ipopt nonideal residual speciation NLP
- Dependency: Ipopt
- Derivative backend: CppAD implicit phase-state sensitivity path where required
- Solver role: Native homogeneous nonideal reactive speciation route
- Implementation owner: src/epcsaft/reactive_speciation.py; src/epcsaft/native/epcsaft_chemical_equilibrium.cpp; src/epcsaft/bindings.cpp
- Validation: tests/api/reactive/test_reactive_speciation.py
- Capability key: reactive_speciation
- Description: Solves nonideal homogeneous reactive speciation through native Ipopt.
- Change note: Initial algorithm-registry entry for nonideal speciation.
- LaTeX: `docs/latex/algorithms.tex:178`
- Code owners: `src/epcsaft/bindings.cpp:3607` (m.def("_solve_chemical_equilibrium_native", &solve_chemical_equilibrium_native_binding);), `src/epcsaft/equilibrium.py:489` (@dataclass(frozen=True, slots=True)), `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp:2419` (ChemicalEquilibriumResultNative solve_nonideal_speciation_chemical_equilibrium_ipopt(), `src/epcsaft/reactive_speciation.py:625` (def _solve_reactive_speciation_native()

This entry covers nonideal standard-state routes whose residuals require
package thermodynamic activity/fugacity calculations.

**LaTeX source**

```tex
This entry covers nonideal standard-state routes whose residuals require
package thermodynamic activity/fugacity calculations.
```


### `reactive_lle_liquid_root_ipopt`
- Family: Reactive phase equilibrium
- Status: Implemented
- Public API: mix.equilibrium(kind="reactive_lle", ...); ReactivePhaseEquilibriumProblem
- Backend: Native C++ Ipopt reactive liquid-root phase-equilibrium NLP
- Dependency: Ipopt
- Derivative backend: Native route residual Jacobian
- Solver role: Coupled reactive LLE liquid-root constrained route
- Implementation owner: src/epcsaft/equilibrium.py; src/epcsaft/native/equilibrium/reactive_phase_equilibrium_problem.cpp; src/epcsaft/bindings.cpp
- Validation: tests/api/reactive/test_reactive_phase_equilibrium.py
- Capability key: reactive_lle_liquid_root
- Description: Solves coupled neutral reactive LLE through the native liquid-root Ipopt route.
- Change note: Initial algorithm-registry entry for reactive LLE.
- LaTeX: `docs/latex/algorithms.tex:195`
- Code owners: `src/epcsaft/bindings.cpp:2301` (m.def("_native_reactive_lle_eos_nlp_contract", [](), `src/epcsaft/equilibrium.py:625` (@dataclass(frozen=True, slots=True)), `src/epcsaft/native/equilibrium/reactive_phase_equilibrium_problem.cpp:2960` (epcsaft::native::equilibrium_nlp::ReactiveTwoPhaseEosRouteResult solve_reactive_phase_liquid_root_route_native(), `src/epcsaft/native/equilibrium_nlp/route_builders.cpp:2230` (ReactiveTwoPhaseEosRouteResult solve_reactive_lle_eos_route()

This entry is a coupled reactive phase-equilibrium route, not a staged
speciation-then-split workflow.

**LaTeX source**

```tex
This entry is a coupled reactive phase-equilibrium route, not a staged
speciation-then-split workflow.
```


### `reactive_electrolyte_lle_liquid_root_ipopt`
- Family: Reactive electrolyte phase equilibrium
- Status: Implemented
- Public API: mix.equilibrium(kind="reactive_electrolyte_lle", ...); ReactivePhaseEquilibriumProblem
- Backend: Native C++ Ipopt reactive electrolyte liquid-root phase-equilibrium NLP
- Dependency: Ipopt
- Derivative backend: Native route residual Jacobian
- Solver role: Coupled reactive electrolyte LLE liquid-root constrained route
- Implementation owner: src/epcsaft/equilibrium.py; src/epcsaft/native/equilibrium/reactive_phase_equilibrium_problem.cpp; src/epcsaft/bindings.cpp
- Validation: tests/api/reactive/test_reactive_phase_equilibrium.py
- Capability key: reactive_electrolyte_lle_liquid_root
- Description: Solves coupled reactive electrolyte LLE through the native liquid-root Ipopt route.
- Change note: Initial algorithm-registry entry for reactive electrolyte LLE.
- LaTeX: `docs/latex/algorithms.tex:212`
- Code owners: `src/epcsaft/bindings.cpp:2338` (m.def("_native_reactive_electrolyte_lle_eos_nlp_contract", [](), `src/epcsaft/equilibrium.py:625` (@dataclass(frozen=True, slots=True)), `src/epcsaft/native/equilibrium/reactive_phase_equilibrium_problem.cpp:2960` (epcsaft::native::equilibrium_nlp::ReactiveTwoPhaseEosRouteResult solve_reactive_phase_liquid_root_route_native(), `src/epcsaft/native/equilibrium_nlp/route_builders.cpp:2277` (ReactiveTwoPhaseEosRouteResult solve_reactive_electrolyte_lle_eos_route()

This entry preserves the boundary that the package returns generic reactive
phase outputs; downstream projects compute application metrics from those
outputs.

**LaTeX source**

```tex
This entry preserves the boundary that the package returns generic reactive
phase outputs; downstream projects compute application metrics from those
outputs.
```


## Regression Algorithms

### `pure_neutral_ceres_regression`
- Family: Regression
- Status: Implemented
- Public API: fit_pure_parameters(...); fit_pure_neutral(...)
- Backend: Native C++ Ceres regression
- Dependency: Ceres
- Derivative backend: Native residual Jacobian for supported target families
- Solver role: Native least-squares pure-neutral parameter regression
- Implementation owner: src/epcsaft/regression.py; src/epcsaft/native/epcsaft_regression.cpp; src/epcsaft/bindings.cpp
- Validation: tests/regression/test_pure_neutral_regression.py
- Capability key: regression:pure_neutral
- Description: Fits currently supported pure-neutral parameter targets through native Ceres.
- Change note: Initial algorithm-registry entry for pure-neutral regression.
- LaTeX: `docs/latex/algorithms.tex:232`
- Code owners: `src/epcsaft/bindings.cpp:3596` (m.def("_fit_pure_neutral_native_ceres", &fit_pure_neutral_native_ceres_binding);), `src/epcsaft/native/epcsaft_regression.cpp:578` (class PureNeutralCeresCostFunction final : public ceres::CostFunction {), `src/epcsaft/regression.py:2386` (def fit_pure_neutral(), `src/epcsaft/regression.py:2723` (def fit_pure_parameters()

This entry covers the implemented pure-neutral native Ceres route and does not
broaden support to every possible pure-component parameter family.

**LaTeX source**

```tex
This entry covers the implemented pure-neutral native Ceres route and does not
broaden support to every possible pure-component parameter family.
```


### `pure_ion_ceres_regression`
- Family: Regression
- Status: Implemented
- Public API: fit_liquid_electrolyte_parameters(...); fit_pure_ion(...)
- Backend: Native C++ Ceres regression
- Dependency: Ceres
- Derivative backend: CppAD implicit residual Jacobian for supported pure-ion targets
- Solver role: Native least-squares pure-ion parameter regression
- Implementation owner: src/epcsaft/regression.py; src/epcsaft/native/epcsaft_regression.cpp; src/epcsaft/bindings.cpp
- Validation: tests/regression/test_pure_ion_regression.py
- Capability key: regression:pure_ion
- Description: Fits currently supported pure-ion and Born-related target sets through native Ceres.
- Change note: Initial algorithm-registry entry for pure-ion regression.
- LaTeX: `docs/latex/algorithms.tex:249`
- Code owners: `src/epcsaft/bindings.cpp:3601` (m.def("_fit_generic_native_ceres", &fit_generic_native_ceres_binding);), `src/epcsaft/native/epcsaft_regression.cpp:1537` (class PureIonCeresCostFunction final : public ceres::CostFunction {), `src/epcsaft/regression.py:2771` (def fit_pure_ion(), `src/epcsaft/regression.py:2974` (def fit_liquid_electrolyte_parameters()

This entry is limited to the currently implemented pure-ion target surface.

**LaTeX source**

```tex
This entry is limited to the currently implemented pure-ion target surface.
```


### `binary_kij_ceres_regression`
- Family: Regression
- Status: Implemented with caveat
- Public API: fit_binary_parameters(...)
- Backend: Native C++ Ceres regression
- Dependency: Ceres
- Derivative backend: CppAD-backed binary k_ij residual Jacobian
- Solver role: Native least-squares constant-k_ij binary regression
- Implementation owner: src/epcsaft/regression.py; src/epcsaft/native/epcsaft_regression.cpp; src/epcsaft/bindings.cpp
- Validation: tests/regression/literature/test_literature_binary_kij_regression.py
- Capability key: regression:binary_pair
- Description: Fits the currently implemented constant-k_ij binary parameter route through native Ceres.
- Change note: Initial algorithm-registry entry keeps l_ij and k_hb_ij out of the claim until implementation evidence exists.
- LaTeX: `docs/latex/algorithms.tex:265`
- Code owners: `src/epcsaft/bindings.cpp:3601` (m.def("_fit_generic_native_ceres", &fit_generic_native_ceres_binding);), `src/epcsaft/native/epcsaft_regression.cpp:1660` (class BinaryKijCeresCostFunction final : public ceres::CostFunction {), `src/epcsaft/regression.py:2800` (def fit_binary_parameters()

This entry intentionally does not claim native optimizer support for every
binary parameter family, including `l_{ij}` or `k^{hb}_{ij}`.

**LaTeX source**

```tex
This entry intentionally does not claim native optimizer support for every
binary parameter family, including \(l_{ij}\) or \(k^{hb}_{ij}\).
```

**Rendered formulae**

$$
l_{ij}
$$

$$
k^{hb}_{ij}
$$


### `reactive_electrolyte_batch_residual_context`
- Family: Reactive regression context
- Status: Implemented residual evaluator
- Public API: build_reactive_regression_objective(...); evaluate_reactive_regression_objective(...)
- Backend: Python orchestration over package thermodynamic calls
- Dependency: None for the context itself
- Derivative backend: Not a registered native optimizer derivative path
- Solver role: Structured residual/context evaluator, not a production Ceres optimizer
- Implementation owner: src/epcsaft/reactive_regression.py; src/epcsaft/runtime.py
- Validation: tests/api/reactive/test_reactive_regression.py
- Capability key: reactive_electrolyte_batch_context
- Description: Evaluates mixed reactive-electrolyte residual contexts without claiming production optimizer support.
- Change note: Initial algorithm-registry entry preserves the current residual-evaluator claim boundary.
- LaTeX: `docs/latex/algorithms.tex:282`
- Code owners: `src/epcsaft/reactive_regression.py:756` (def build_reactive_regression_objective(), `src/epcsaft/reactive_regression.py:787` (def evaluate_reactive_regression_objective(), `src/epcsaft/runtime.py:542` ("reactive_electrolyte_batch_context": {)

This entry records the implemented residual/context evaluator and explicitly
does not claim native Ceres production optimization for reactive electrolyte
batch fitting.

**LaTeX source**

```tex
This entry records the implemented residual/context evaluator and explicitly
does not claim native Ceres production optimization for reactive electrolyte
batch fitting.
```
