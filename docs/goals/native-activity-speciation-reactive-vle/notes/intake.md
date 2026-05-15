# Stage 0 Intake Gate

This note is the required issue #115 Stage 0 intake artifact. It must be completed by active task `T001` before any source edits happen.

## GitHub Source

- Issue: https://github.com/tannerpolley/ePC-SAFT/issues/115
- Title: Implement native coupled activity speciation and reactive VLE pressure benchmark
- Relevant downstream comment: https://github.com/tannerpolley/ePC-SAFT/issues/115#issuecomment-4457708894

## Current Reactive / Speciation Solver Files

- `src/epcsaft/reactive_speciation.py`
  - Public `solve_reactive_speciation(...)` entry point starts at line 318.
  - Activity-coupled nonideal standard states currently divert through `_solve_reactive_speciation_activity_fixed_point(...)` from line 667.
  - Native route `_solve_reactive_speciation_native(...)` starts at line 1132 and calls `_core._solve_chemical_equilibrium_native(...)`.
- `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp`
  - Current native homogeneous chemical-equilibrium/speciation implementation owner.
  - ePC-SAFT activity/fugacity coefficient evaluation helpers are present at lines 168 and 204.
  - Residual evaluation evaluates activity coefficients inside the native residual path around line 636.
  - Current derivative selection allows analytic ideal mole-fraction derivatives and rejects activity/concentration-coupled derivative support around lines 270-289.
  - Native residual and solve entry points start at lines 791 and 963.
- `src/epcsaft/native/epcsaft_chemical_equilibrium.h`
  - Native option/result declarations for chemical equilibrium, residual evaluation, activity coefficients, and Jacobian payloads.
- `src/epcsaft/bindings.cpp`
  - Pybind bridge exposes `_solve_chemical_equilibrium_native` and `_evaluate_chemical_equilibrium_residual_native` at lines 1591-1592.
- `src/epcsaft/reactive_electrolyte.py`
  - Current reactive pressure path is staged: `solve_reactive_speciation(...)` followed by `electrolyte_bubble_pressure(...)`.
- `src/epcsaft/reactive_staged.py`
  - Current generic staged reactive-equilibrium helper.
- `src/epcsaft/epcsaft.py`
  - Mixture methods route public calls to reactive/speciation and staged equilibrium helpers.
- No current `src/epcsaft/native/reactive_speciation.*` files were found; the active native speciation owner is `epcsaft_chemical_equilibrium.*`.

## Current VLE / Bubble / Dew Files

- `src/epcsaft/equilibrium.py`
  - `BubblePointProblem`, `DewPointProblem`, `ReactiveSpeciationProblem`, `ReactiveElectrolyteBubbleProblem`, and `ReactivePhaseEquilibriumProblem`.
  - Neutral bubble/dew outer solve starts at `_neutral_bubble_dew_outer(...)` around line 1502 and result shaping at line 1747.
- `src/epcsaft/electrolyte_bubble.py`
  - Python wrapper/result layer for fixed-liquid electrolyte bubble pressure.
- `src/epcsaft/native/epcsaft_equilibrium.cpp`
  - Native equilibrium owner for electrolyte bubble pressure, LLE, reactive LLE, reactive electrolyte LLE, fugacity residuals, and partial-pressure diagnostics.
  - Electrolyte bubble evaluation starts around line 818; native electrolyte bubble entry point starts around line 1919.
- `src/epcsaft/native/epcsaft_equilibrium.h`
  - Native equilibrium option/result declarations.
- `src/epcsaft/native/equilibrium/equilibrium_helpers.cpp`
- `src/epcsaft/native/equilibrium/equilibrium_helpers.h`
- `src/epcsaft/reactive_electrolyte.py`
  - Current reactive bubble route is staged speciation plus fixed-liquid electrolyte bubble pressure, not simultaneous reactive VLE.

## Current Python Wrapper Files

- `src/epcsaft/equilibrium.py`
  - Generic problem dataclasses and neutral/electrolyte equilibrium functions.
- `src/epcsaft/reactive.py`
  - Public reactive import surface.
- `src/epcsaft/reactive_speciation.py`
  - Public speciation API, native-call normalization, activity fixed-point route, and result shaping.
- `src/epcsaft/reactive_electrolyte.py`
  - Reactive speciation followed by electrolyte bubble pressure.
- `src/epcsaft/reactive_staged.py`
  - Staged chemical equilibrium followed by explicit phase route.
- `src/epcsaft/electrolyte_bubble.py`
  - Electrolyte bubble wrapper/result layer.
- `src/epcsaft/runtime.py`
  - Capability metadata currently describes fixed-point activity handling for reactive speciation and staged reactive routes.
- `src/epcsaft/epcsaft.py`
  - `ePCSAFTMixture` public methods and `equilibrium(...)` dispatcher.
- `src/epcsaft/bindings.cpp`
  - Pybind bridge to native chemical/equilibrium routines.

## Current Tests That Use Staged Or Fixed-Point Behavior

- `tests/api/reactive/test_reactive_speciation_results.py`
  - Asserts activity fixed-point diagnostics and derivative policy that is not coupled into the outer activity iteration.
- `tests/api/reactive/test_reactive_staged_equilibrium.py`
  - Exercises staged chemical equilibrium followed by a phase route.
- `tests/api/reactive/test_reactive_staged_workflow_contract.py`
  - Asserts staged workflow diagnostics and a public generic staged contract.
- `tests/equilibrium/reactive/test_reactive_lle.py`
  - Current reactive LLE test asserts staged coupling diagnostics.
- `tests/native/equilibrium/test_chemical_equilibrium_native_residuals.py`
  - Current native residual coverage includes activity fixed-point expectations.
- `tests/api/reactive/test_reactive_electrolyte_bubble_setup.py`
- `tests/api/reactive/test_reactive_electrolyte_bubble_results.py`
  - Current reactive electrolyte bubble coverage follows staged speciation into fixed-liquid electrolyte bubble pressure.
- `tests/api/reactive/test_reactive_speciation_options.py`
  - Covers solver and Jacobian option boundaries for reactive speciation.

## Chosen Pressure / Speciation Benchmark

- Issue #115 requires a generic CO2 + amine + water pressure/speciation case with repo-contained fixture data.
- Current repo-contained starting points:
  - `tests/api/package/test_downstream_integration_smokes.py` has generic `CO2`, `H2O`, `Amine`, `AmineH+`, and `HCO3-` public API smoke rows including speciation and partial-pressure families.
  - `scripts/benchmarks/helpers/reactive_regression.py` has pressure/speciation reactive-regression surrogate cases, including a 35-row pressure/speciation surrogate shape.
  - `tests/fixtures/literature/pure_neutral/mea_co2_h2o_benchmark.json` and `tests/regression/literature/test_mea_co2_h2o_pure_parameter_benchmark.py` provide repo-contained MEA/CO2/H2O parameter-regression fixture material, but not yet the full issue #115 reactive VLE pressure/speciation proof.
- Chosen Stage 5 benchmark target: extend or add a generic CO2 + amine + water native reactive pressure/speciation benchmark based on the public downstream-smoke species family, with repo-contained fixture values and required numeric checks for reaction residual norm, charge residual norm, material residual norm, CO2 partial pressure, liquid speciation, activity convention, derivative backend, and native solver iteration diagnostics.
- The preferred MEA-style nine-species benchmark remains desirable, but Stage 0 did not find an existing repo-contained nine-species pressure/speciation fixture that already satisfies the issue #115 proof contract.

## Intake Completion Checklist

- Current implementation files listed with exact paths.
- Current wrapper files listed with exact paths.
- Current staged or fixed-point tests listed with exact paths.
- Benchmark fixture path or fixture gap identified.
- Native-regression no-edit boundary confirmed.
- Banned-token hygiene checked for this goal directory.

## Stage 0 Evidence Commands

- GitHub connector loaded issue #115 body and the downstream blocker comment.
- JetBrains MCP `ide_index_status` reported indexing ready for `C:/Users/Tanner/Documents/git/ePC-SAFT`.
- JetBrains MCP `ide_find_references` on `ReactiveSpeciationProblem` found public imports and test usage, confirming the public generic API surface.
- `rg --files src/epcsaft/native src/epcsaft tests/native/equilibrium tests/api/reactive tests/equilibrium`
- `rg -n` targeted searches for reactive/speciation, VLE/bubble/dew, staged/fixed-point, benchmark, and pressure/speciation terms.
- Focused file reads of `src/epcsaft/reactive_speciation.py`, `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp`, `src/epcsaft/equilibrium.py`, and `src/epcsaft/reactive_electrolyte.py`.
