## T051 audit

### Question answered

After `T050`, should issue #53 keep pushing the public runtime `dadx()` auto policy toward CppAD, or should that stay analytic by design while the regression-production work continues elsewhere?

### Decision

For issue #53, the public runtime `dadx()` auto/default policy should stay **analytic by design** for now.

The next worker should **not** re-open a broad runtime auto-autodiff implementation pass.

Instead, the next worker should update the public capability/reporting layer so it truthfully reflects what already landed:

- ideal chemical-equilibrium `auto` now prefers `autodiff` when CppAD is compiled
- concentration and activity auto paths are CppAD-backed on the supported native slice
- public runtime `dadx()` auto/default still remains analytic

### Why

The attempted public runtime auto-autodiff shift was not a clean production move:

1. it perturbed established public runtime reference outputs on the focused slice
2. it exposed unsupported/default-sensitive behavior on the SSM+DS Born activity path
3. it did not advance the core issue #53 regression-production objective as directly as the chemical-equilibrium Jacobian work did

That means the runtime `dadx()` auto story is still narrower than the native regression-production story.

### Evidence

- `docs/goals/issue-53-native-regression-production/notes/T050-default-derivative-slice.md`
- `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp`
- `tests/native/test_chemical_equilibrium_native.py`
- `tests/api/test_runtime.py`
- `src/epcsaft/runtime.py`

### Next worker scope

Update the reporting and docs surfaces so they match the implemented truth:

1. `epcsaft.capabilities()["equilibrium"]["reactive_speciation"]`
2. runtime-facing tests for that metadata
3. docs that still describe the old ideal-auto-analytic behavior

### Still later blocker

- reactive-electrolyte bubble differentiation remains the larger unresolved derivative architecture gap after the reporting cleanup
