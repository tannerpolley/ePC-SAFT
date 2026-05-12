## T049 audit

### Question answered

After `T048`, what is the next real derivative scope if the user wants:

1. the stale pressure-composition blocker removed from the board truth, and
2. the codebase pushed further toward CppAD-backed differential surfaces?

### Confirmed facts

The stale blocker

- "The required native pressure-composition derivative substrate is still missing."

is still false.

Evidence:

- `src/epcsaft/native/epcsaft_Z.cpp`
  - `pressure_composition_derivative_result_cpp(...)`
  - `pressure_density_derivative_result_cpp(...)`
- `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp`
  - concentration/activity autodiff paths consume those pressure derivatives
- `tests/native/test_runtime_contracts.py`
  - runtime derivative payloads are already exercised

### What T048 actually closed

`T048` closed the first honest non-`logK` thermo derivative slice for the supported reactive-speciation path:

- SSM+DS Born activity-standard-state row support
- `d_born`
- `f_solv`

It did **not** close the broader "every differential uses CppAD" objective.

### What is still not CppAD-backed by default

1. Public runtime `dadx()` auto/default policy still reports analytic backends

Evidence:

- `tests/api/test_runtime.py::test_default_dadx_reports_auto_derivative_policy`
  - `hc == analytic`
  - `disp == analytic`
  - `ion == analytic`
  - `born == analytic`
  - `assoc == analytic`
- `src/epcsaft/native/epcsaft_ares.cpp`
  - backend map still resolves auto/default contribution modes to `"analytic"`
- `src/epcsaft/native/contributions/epcsaft_contrib_hc.cpp`
- `src/epcsaft/native/contributions/epcsaft_contrib_disp.cpp`
- `src/epcsaft/native/contributions/epcsaft_contrib_ion.cpp`
- `src/epcsaft/native/contributions/epcsaft_contrib_born.cpp`
  - explicit autodiff exists, but auto/default still keeps analytic when available

2. Ideal chemical-equilibrium Jacobian path still defaults to analytic, not CppAD

Evidence:

- `tests/native/test_chemical_equilibrium_native.py::test_native_chemical_equilibrium_residual_evaluator_uses_analytic_jacobian_by_default`
- `src/epcsaft/native/epcsaft_chemical_equilibrium.cpp`
  - auto selection still routes the ideal log-amount case through `analytic_ideal_log_amount_jacobian(...)`

3. Reactive-electrolyte bubble derivatives remain structurally unsupported

Evidence:

- `tests/native/test_native_ceres_reactive_pressure_speciation.py`
- `src/epcsaft/native/regression/thermo_regression.cpp`
  - derivative support gate is still `row_mode == "reactive_speciation"` only

### Scope decision

The next worker tranche should **not** be thermo `k_ij`.

Reason:

- The latest user direction narrowed the native thermo regression focus to the Born SSM+DS regression parameters (`d_born`, `f_solv`).
- Generic binary parameters such as `k_ij` can remain handled by the generic/binary regression path rather than forcing them into the reactive-speciation thermo tranche right now.

### Next tranche

The next concrete worker should target:

1. default/auto runtime composition-derivative policy:
   - prefer CppAD for supported `dadx()` contribution surfaces when CppAD is compiled
2. default/auto ideal chemical-equilibrium Jacobian policy:
   - prefer CppAD on the supported ideal log-amount path
3. keep explicit finite-difference debug gating intact

Bubble differentiation stays as the later blocker after that.

### Completion truth after T049

The goal is still not complete.

Why:

- default runtime derivative surfaces are still not uniformly CppAD-backed
- ideal chemical-equilibrium auto/default is still analytic
- bubble regression differentiation is still missing
