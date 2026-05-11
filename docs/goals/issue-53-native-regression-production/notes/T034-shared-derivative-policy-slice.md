T034 receipt: shared derivative-policy and Born scalar-helper slice
==================================================================

Result
------

Done.

What changed
------------

This slice intentionally stayed structural:

- added a shared native finite-difference debug-gate helper under
  `src/epcsaft/native/autodiff/debug_gate.h`;
- rewired native chemical-equilibrium and native electrolyte-LLE derivative
  gating to use that shared helper instead of local duplicated environment
  checks;
- removed duplicated double versus autodiff solvent-reference dielectric logic
  in the Born contribution path by introducing one scalar-templated local
  helper in `epcsaft_contrib_born.cpp`.

What did not change
-------------------

- no production capability claims were widened;
- no `d_born`, `f_solv`, or `k_ij` parameter sensitivities were claimed as
  implemented for native thermodynamic regression;
- reactive-electrolyte bubble rows remain a later implicit-derivative task.

Why this slice was worth doing
------------------------------

The broader tranche needs one shared derivative-policy surface and less scalar
duplication in the Born path before the next worker attempts real `R_theta`
assembly for reactive-speciation parameter sensitivities.

Validation
----------

- `uv run python scripts/build_epcsaft.py --build-only --parallel 10`
- `uv run python run_pytest.py tests/native/test_cppad_eos_derivatives.py tests/native/test_cppad_reactive_speciation_derivatives.py tests/native/test_cppad_bubble_derivatives.py tests/native/test_chemical_equilibrium_native.py tests/native/test_equilibrium_native_contracts.py -q`

Observed result
---------------

- build-only native rebuild passed in about 16.5s;
- focused native derivative and contract slice passed: `32 passed in 3.31s`.
