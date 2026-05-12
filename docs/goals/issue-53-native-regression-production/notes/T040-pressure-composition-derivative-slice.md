T040 receipt: native pressure-composition derivative slice
=========================================================

Result
------

Completed the first broader EOS derivative architecture slice after the T037/T039
checkpoint.

What changed
------------

Added a backend-unavailable-free native pressure-composition derivative result
surface for the currently supported nonassociating EOS slice.

Implemented surfaces:

- `PressureCompositionDerivativeResult`
- `pressure_composition_derivative_result_cpp(...)`
- `ePCSAFTStateNative::pressure_composition_derivative_result()`
- Python state wrapper: `state.pressure_composition_derivative()`

Scope
-----

Supported:

- neutral nonassociating states
- ionic nonassociating states
- ionic pressure derivatives through HC + DISP + DH terms
- Born terms remain pressure-neutral because `dadrho_born_cpp()` is zero in the
  current model

Explicitly still gated:

- associating states
- any tranche that still needs backend-unavailable-free `dp/drho`
- concentration-standard-state reactive-speciation derivatives in native thermo
  regression

Why this matters
----------------

The T037 blocker said there was no native pressure-composition derivative
substrate at all. That is no longer true for the supported nonassociating slice.

The remaining blocker is now narrower and more precise:

- native concentration-standard-state regression still needs a real
  backend-unavailable-free `dp/drho` path through density closure

Validation evidence
-------------------

Focused runtime/native validation passed after rebuilding `_core`:

- `uv run python scripts/build_epcsaft.py --build-only --parallel 10`
- `uv run python run_pytest.py tests/native/test_runtime_contracts.py tests/native/test_native_ceres_thermodynamic_regression.py -q`

Observed result:

- `29 passed in 2.01s`

Key test coverage added
-----------------------

- supported-state `pressure_composition_derivative()` matches constrained
  backend-unavailable pressure changes
- associating ionic runtime reports the explicit nonassociating gate instead of
  pretending support


