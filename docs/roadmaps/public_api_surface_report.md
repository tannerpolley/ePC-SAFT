# Public API Surface Report

This report is the Phase 9 deliverable for issue #120. It audits the public
package surface without removing documented imports or changing runtime
behavior.

## Summary

- `epcsaft.__all__` currently exposes 104 names.
- `pyproject.toml` defines no `project.scripts` and no entry-point groups.
- `python -m epcsaft` is supported through `src/epcsaft/__main__.py`.
- The package intentionally provides both direct top-level exports and organized
  public namespace modules such as `epcsaft.eos`, `epcsaft.equilibrium`,
  `epcsaft.electrolyte`, `epcsaft.reactive`, `epcsaft.regression`, and
  `epcsaft.diagnostics`.
- No export was removed in this cleanup phase. Removals require a deprecation
  or compatibility-shim decision because tests, docs, and downstream smoke
  checks currently exercise the broad surface.

## Top-Level Export Groups

### Stable Core And Runtime API

- Classes and exceptions: `ePCSAFTMixture`, `ePCSAFTState`,
  `ActivityCoefficientResult`, `InputError`, `SolutionError`.
- Runtime metadata: `__version__`, `__git_commit__`, `runtime_build_info`,
  `capabilities`.
- Fast property helpers: `evaluate_fugacity_coefficients`,
  `evaluate_fugacity_coefficients_batch`.

### Parameter And Dataset API

- Dataset access: `DATASET_ROOT`, `available_datasets`, `get_prop_dict`,
  `validate_dataset_bundle`.
- User parameter folders: `create_parameter_template`.
- Schema objects: `AssociationSite`, `BinaryRecord`, `ComponentIdentifier`,
  `ParameterSet`, `PermittivityRecord`, `PureRecord`.
- Electrolyte composition helpers: `molality_to_molefraction`,
  `molefraction_to_molality`.

### Phase Equilibrium API

- Problem/result objects: `EquilibriumProblem`, `EquilibriumOptions`,
  `EquilibriumPhase`, `EquilibriumResult`, `TPFlash`, `LLEProblem`,
  `BubblePoint`, `DewPoint`, `StabilityAnalysis`, `StabilityResult`,
  `StabilityTrial`, `ElectrolyteLLEProblem`, `ElectrolyteBubblePoint`,
  `ReactiveSpeciationProblem`, `ReactivePhaseEquilibriumProblem`,
  `ReactiveElectrolyteBubbleProblem`.
- Route helpers: `bubble_p`, `bubble_t`, `dew_p`, `dew_t`,
  `electrolyte_feed_from_molality`, `initial_phases_from_result`.

### Reactive Workflow API

- Reactive speciation: `ReactionConstantConvention`, `ReactionDefinition`,
  `ReactiveSpeciationOptions`, `ReactiveSpeciationResult`,
  `solve_reactive_speciation`, `solve_reactive_speciation_sweep`.
- Reactive electrolyte bubble: `ReactiveElectrolyteBubbleOptions`,
  `ReactiveElectrolyteBubbleResult`, `solve_reactive_electrolyte_bubble`,
  `solve_reactive_electrolyte_bubble_sweep`.
- Sequential reactive workflow: `ReactiveStagedEquilibriumResult`,
  `solve_reactive_staged_equilibrium`. These remain public but must be described
  as sequential coupling, not as a full coupled reactive flash.

### Regression API

- Generic fit objects: `FitParameter`, `FitBounds`, `FitTerm`, `FitProblem`,
  `FitResult`, `BinaryInteraction`, `RelativePermittivityResidual`,
  `TargetDataset`, `TargetRow`.
- Fit helpers: `fit_pure_parameters`, `fit_binary_parameters`,
  `fit_liquid_electrolyte_parameters`, `fit_pure_neutral`, `fit_pure_ion`,
  `fit_binary_pair`, `load_regression_records`, `validate_regression_provenance`,
  `write_fit_result`.
- Derivative/result helpers: `evaluate_pure_neutral_derivatives`,
  `evaluate_generic_regression_derivatives`,
  `evaluate_reactive_electrolyte_bubble_residuals`.
- Reactive electrolyte regression: `ReactiveElectrolyteRegressionResult`,
  `ReactiveElectrolyteBatch`, `ReactiveElectrolyteBatchOptions`,
  `ReactiveElectrolyteBatchResult`, `ReactiveElectrolyteRow`,
  `ReactiveElectrolyteRowResult`, `ReactiveElectrolyteRegressionContext`,
  `ReactiveRegressionObjective`, `ReactiveRegressionObjectiveResult`,
  `ReactiveRegressionJacobianResult`, `ReactiveRegressionFitResult`,
  `build_reactive_regression_objective`,
  `evaluate_reactive_regression_objective`,
  `fit_reactive_electrolyte_parameters`, `summarize_regression_result`,
  `write_regression_parameter_table`, `write_regression_residual_table`,
  `write_regression_row_table`, `write_regression_summary`.

### Sensitivity API

- Public implicit-sensitivity helpers: `ImplicitSolveResult`,
  `implicit_backend_for_residual_backend`,
  `implicit_sensitivity_from_jacobians`.
- The previous helper that created missing implicit-derivative payloads is no
  longer exported; unsupported implicit derivative requests now raise instead
  of returning a placeholder result.

## Top-Level Module Classification

| Module | Classification | Decision |
| --- | --- | --- |
| `epcsaft` | Stable public package root | Keep. |
| `epcsaft.__main__` | Public CLI module for `python -m epcsaft` | Keep; no project script entrypoint is declared. |
| `epcsaft.epcsaft` | Historical core class module | Keep as public/legacy import path. |
| `epcsaft.eos` | Organized public import namespace | Keep. |
| `epcsaft.equilibrium` | Stable public equilibrium API | Keep. |
| `epcsaft.electrolyte` | Organized public electrolyte namespace | Keep. |
| `epcsaft.electrolyte_bubble` | Stable public fixed-liquid bubble contracts | Keep. |
| `epcsaft.reactive` | Organized public reactive namespace | Keep. |
| `epcsaft.reactive_speciation` | Stable public reactive speciation API | Keep. |
| `epcsaft.reactive_electrolyte` | Stable public sequential reactive electrolyte bubble API | Keep. |
| `epcsaft.reactive_staged` | Public sequential workflow API | Keep, but wording must remain explicit about sequential coupling. |
| `epcsaft.reactive_regression` | Public reactive regression API | Keep. |
| `epcsaft.regression` | Public generic regression API | Keep. |
| `epcsaft.parameters` | Public parameter dataset API | Keep. |
| `epcsaft.parameter_schema` | Public schema API | Keep. |
| `epcsaft.parameter_templates` | Public user-parameter-folder API | Keep. |
| `epcsaft.properties` | Public fast property-helper API | Keep. |
| `epcsaft.dataset_validation` | Public validation helper API | Keep. |
| `epcsaft.diagnostics` | Organized public diagnostics namespace | Keep. |
| `epcsaft.runtime` | Public runtime metadata implementation module | Keep. |
| `epcsaft.implicit_sensitivity` | Public sensitivity helper module | Keep, but status-token helper naming should be revisited only with API migration. |
| `epcsaft._types` | Private implementation module | Do not promote. |
| `epcsaft.equilibrium_core` | Internal equilibrium support package | Keep internal; do not add top-level exports. |

## Entry Points

`pyproject.toml` currently declares:

- `project.scripts`: none.
- `project.entry-points`: none.
- Supported command-style surface: `python -m epcsaft`.

No entrypoint changes are needed for issue #120.

## Tightening Decisions

- Public namespaces are now documented in `docs/pages/package_architecture.rst`.
- Benchmark runtime logic has already moved out of the installable runtime
  package; benchmark helpers are imported from `scripts/benchmarks/helpers`.
- Native IPOPT implementation work belongs under native equilibrium internals;
  no Python optional backend module remains in the runtime package.
- Sequential reactive workflow names remain public for compatibility, but docs
  now state the workflow boundary explicitly.
- No application-specific downstream helper was found in `epcsaft.__all__` that
  can be safely hidden without a deprecation decision.

## Future API Work

- Decide whether the legacy implicit-sensitivity status helper should be renamed
  or retained permanently.
- Decide whether `DATASET_ROOT` should stay public or move behind a documented
  dataset accessor in a future major cleanup.
- Do not use issue #120 to remove public names that tests or docs currently
  exercise.
