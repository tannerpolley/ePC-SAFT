T047 audit: post-activity derivative disposition
===============================================

Status
------

T047 is complete.

Audit conclusion
----------------

The user-added blocker statement

- "The required native pressure-composition derivative substrate is still missing."

is no longer true.

That substrate is already present and validated:

- `src/epcsaft/native/epcsaft_Z.cpp`
  - `pressure_composition_derivative_result_cpp(...)`
- `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp`
  - concentration and activity standard-state autodiff paths consume native
    pressure-closure derivatives
- `tests/native/test_runtime_contracts.py`
  - runtime contracts already validate the derivative payloads

What is now proven
------------------

The supported native Ceres thermodynamic regression hot loop now has three real
CppAD-backed implicit derivative benchmark cases:

1. ideal mole-fraction `logK`
2. concentration-standard-state `logK`
3. activity-standard-state `logK`

Validation evidence:

- `uv run python scripts/benchmark_native_ceres_thermo_regression.py --case reactive_speciation_logk_implicit --warmup 1 --repeat 1`
  - `converged ceres cppad_implicit`
- `uv run python scripts/benchmark_native_ceres_thermo_regression.py --case reactive_speciation_concentration_logk_implicit --warmup 1 --repeat 1`
  - `converged ceres cppad_implicit`
- `uv run python scripts/benchmark_native_ceres_thermo_regression.py --case reactive_speciation_activity_logk_implicit --warmup 1 --repeat 1`
  - `converged ceres cppad_implicit`

What is still not done
----------------------

Two production derivative surfaces remain explicitly unsupported:

1. Reactive-electrolyte bubble regression differentiation

   Evidence:
   - `tests/native/test_native_ceres_reactive_pressure_speciation.py`
     - row diagnostics still report `derivative_backend == "not_differentiated"`
     - fitting still reports `status == "backend_unavailable"`
   - `src/epcsaft/native/regression/thermo_regression.cpp`
     - derivative support gate still only accepts `row_mode == "reactive_speciation"`

2. Non-`logK` thermodynamic parameter sensitivities on the supported
   `reactive_speciation` slice

   Evidence:
   - `tests/native/test_native_ceres_thermodynamic_regression.py::test_native_thermo_regression_reports_ssmds_born_derivatives_unavailable`
   - `src/epcsaft/native/regression/thermo_regression.cpp`
     - parameter support gate still rejects everything except
       `reaction_equilibrium_constant` / `log_equilibrium_constant`

Decision
--------

The next concrete tranche should be:

- supported non-`logK` thermodynamic parameter sensitivities for native
  `reactive_speciation` rows first

Reasoning:

- It stays on the already-working `reactive_speciation` regression surface.
- It directly addresses the remaining parameter kinds called out by the issue
  and current tests (`d_born`, `f_solv`, `k_ij`).
- Bubble differentiation still requires a broader phase-equilibrium implicit
  derivative system and remains the larger architecture jump.

Rejected next steps
-------------------

- Bubble first:
  rejected because the row mode is still structurally undifferentiated and
  needs a larger solver-coupling tranche than the remaining
  `reactive_speciation` parameter work.

- Stop and claim completion:
  rejected because production derivatives still do not cover bubble rows or
  non-`logK` thermodynamic parameters.
