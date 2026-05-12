## T050 result

### Outcome

This tranche completed the safe default-policy move for **ideal chemical-equilibrium** Jacobians, but it did **not** complete the broader public runtime `dadx()` auto-policy shift.

### What landed

- `chemical_equilibrium` with `jacobian_backend="auto"` now prefers `autodiff` on the supported ideal mole-fraction log-amount slice when CppAD is compiled.
- The explicit finite-difference debug-only gate remains unchanged.

Files:

- `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp`
- `tests/native/test_chemical_equilibrium_native.py`

### What was tried and reverted

I also tested pushing the public runtime `dadx()` auto/default contribution modes toward autodiff for HC/dispersion/ion/Born/association.

That broader switch was **not** retained because it was not a clean production move:

1. It perturbed established public runtime reference values on the focused test slice.
2. It exposed unsupported/default-sensitive behavior on the SSM+DS Born activity path.
3. The supported/public `dadx()` auto story is therefore still narrower than the chemical-equilibrium auto story.

So the runtime `dadx()` auto/default route remains analytic by default for now.

### Validation

- `uv run python scripts/build_epcsaft.py --build-only --parallel 10`
- `uv run python run_pytest.py tests/api/test_runtime.py tests/api/test_reactive_speciation.py tests/native/test_chemical_equilibrium_native.py tests/native/test_cppad_reactive_speciation_derivatives.py -q`

Result:

- `80 passed, 1 skipped`

### Board truth after T050

- Safe and validated:
  - ideal chemical-equilibrium `auto -> autodiff` on the supported log-amount slice
- Still not safe to claim:
  - public runtime `dadx()` auto/default prefers autodiff in production
- Still later blocker:
  - reactive-electrolyte bubble differentiation
