T028 audit: issue #53 completion classification
================================================

Classification
--------------

Partial production slice. The branch has real native-regression production
progress, but it does not fully complete issue #53.

Completed production pieces
---------------------------

- Native thermodynamic regression rows can be evaluated in C++ with fixed-shape
  residuals.
- The supported Ceres thermodynamic fit slice is wired through
  `fit_native_thermo_regression(...)` and now reachable from the high-level
  `fit_reactive_electrolyte_parameters(..., backend="native")` wrapper for
  reactive-speciation/logK/linear-speciation-target contexts.
- Python compatibility fitting remains explicit behind `backend="python_compat"`
  and diagnostics label it non-production.
- Production Backend unavailables remain gated/rejected rather than silently used.
- Generic binary parameters are documented as a separate `fit_binary_pair(...)`
  path; reactive-regression tests for the Born SSM+DS gap now focus on
  `d_born` and `f_solv`.

Remaining blockers
------------------

- Born-SSM+DS `d_born` and `f_solv` parameters are applied to native mixtures,
  but production Ceres sensitivities through the activity/fugacity path are not
  implemented. The branch correctly reports `backend_unavailable`.
- Reactive electrolyte bubble-pressure rows are native-evaluated but not
  production-differentiated for Ceres. The coupled log-pressure,
  vapor-composition, fugacity-equality residual system still needs analytic,
  CppAD, or implicit sensitivities.
- The local build reports Ceres and CppAD disabled, so benchmark evidence for
  the Ceres loop in this checkout is fixed-shape/native-boundary evidence, not
  a successful local Ceres solve.

PR recommendation
-----------------

Commit this branch as an honest partial native-regression production slice.
Do not close issue #53 or claim full mixed pressure/speciation MEA-style
production regression until the remaining sensitivity gaps are implemented and
validated. A PR title should explicitly say "partial native thermodynamic
regression slice" or similar.

Validation evidence
-------------------

- `uv run python run_pytest.py tests/native/test_native_ceres_thermodynamic_regression.py tests/native/test_cppad_eos_derivatives.py tests/api/test_runtime.py tests/native/test_native_ceres_reactive_pressure_speciation.py tests/native/test_cppad_bubble_derivatives.py tests/api/test_reactive_regression.py tests/workflows/test_benchmark_reactive_regression.py tests/native/test_native_reactive_regression.py tests/native/test_native_regression_types.py tests/workflows/test_benchmark_native_regression.py -q`
  - pass: 80 passed in 110.89s
- `uv run python scripts/benchmark_native_ceres_thermo_regression.py --warmup 1 --repeat 3`
  - pass: native thermodynamic benchmark ran and reported `backend_unavailable`
    honestly for the current Ceres-disabled local build
- `uv run python scripts/validate_project.py quick`
  - pass: 23 passed in 12.37s after doctor
- `uv run python scripts/validate_project.py docs`
  - pass
- `uv run ruff check ...`
  - pass on touched Python files
- `uv run black --check ...`
  - pass on touched Python files

