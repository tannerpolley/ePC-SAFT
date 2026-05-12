## T053 audit

### Question answered

After `T052`, which derivative surfaces still do not use CppAD in production, which ones are intentionally not on CppAD by default, and what is the next honest implementation tranche?

### Confirmed CppAD-backed production surfaces

The supported native thermodynamic regression slices are already on the real CppAD production path:

- `reactive_speciation_logk_implicit`
- `reactive_speciation_concentration_logk_implicit`
- `reactive_speciation_activity_logk_implicit`
- `reactive_speciation_activity_ssmds_born_radius_implicit`

These benchmark as `backend=ceres`, `derivative=cppad_implicit`, and `native_hot_loop=True`.

### Deliberate analytic-by-design surfaces

The public runtime `state.dadx()` auto/default policy is still intentionally analytic for issue #53.

That is a policy choice, not a hidden missing-production-derivative bug. `T051` already rejected a broad public auto-autodiff flip because it perturbed established runtime reference outputs and exposed unsupported/default-sensitive behavior on the SSM+DS Born activity path.

### Remaining non-CppAD derivative gaps

1. Public/runtime EOS derivative helpers still implemented on the older `AutoDual` substrate rather than the CppAD substrate.
   - `src/epcsaft/native/epcsaft_Z.cpp`
   - `src/epcsaft/native/epcsaft_parameter_setup.cpp`
   - `src/epcsaft/native/epcsaft_regression.cpp`
   - `src/epcsaft/native/contributions/epcsaft_contrib_internal.h`

2. Reactive-electrolyte bubble differentiation is still explicitly unsupported in native thermo regression.
   - bubble rows still report `derivative_backend = "not_differentiated"`
   - bubble fitting still returns `backend_unavailable`

### Decision

The next honest tranche is **not** to reopen runtime auto/default policy, and it is **not** to jump straight into bubble differentiation first.

The next implementation tranche should be:

1. migrate the remaining supported public/runtime derivative helpers from `AutoDual` to the native CppAD substrate for supported nonassociating states
2. keep public runtime `dadx()` auto/default analytic unless explicit runtime-contract validation later justifies changing policy
3. leave reactive-electrolyte bubble differentiation as the next larger downstream blocker after the substrate migration

### Why this is the next tranche

- It matches the updated goal direction to push real differential surfaces onto CppAD where differentials are used.
- It is narrower and safer than the full bubble-pressure implicit residual/Jacobian architecture jump.
- It removes a real remaining substrate inconsistency without pretending the bubble path is solved.

### Evidence

- `src/epcsaft/runtime.py`
- `src/epcsaft/native/epcsaft_Z.cpp`
- `src/epcsaft/native/epcsaft_parameter_setup.cpp`
- `src/epcsaft/native/epcsaft_regression.cpp`
- `src/epcsaft/native/regression/thermo_regression.cpp`
- `tests/native/test_native_ceres_reactive_pressure_speciation.py`
- `docs/goals/issue-53-native-regression-production/notes/T051-runtime-auto-scope-audit.md`
- `docs/goals/issue-53-native-regression-production/notes/T052-runtime-reporting-alignment.md`
