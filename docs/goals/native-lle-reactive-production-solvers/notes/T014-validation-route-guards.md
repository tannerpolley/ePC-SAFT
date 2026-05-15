# T014 Validation And Route Guard Receipt

## Result

Done. The issue #117 Stage 10 validation ladder passed on the current branch after the T013 public route commit.

## Validation

- `uv run python scripts/dev/build_epcsaft.py --clean --enable-ceres --enable-cppad`: pass.
- `uv run python run_pytest.py tests/native/equilibrium -q`: pass, 26 tests.
- `uv run python run_pytest.py tests/equilibrium/reactive -q`: pass, 4 tests.
- `uv run python run_pytest.py tests/api/reactive -q`: pass, 65 passed and 1 skipped for optional `cyipopt`.
- `uv run python scripts/dev/validate_project.py quick`: pass, 31 tests.
- `uv run python scripts/dev/validate_project.py docs`: pass.
- `git diff --check`: pass.

## Combined Issue #116 Checks

The issue #116 validation slices that are outside the issue #117 ladder also passed on the same clean Ceres+CppAD build:

- `uv run python run_pytest.py tests/equilibrium/electrolyte -q`: pass, 37 passed and 7 skipped opt-in/optional cases.
- `uv run python run_pytest.py tests/equilibrium/core -q`: pass, 74 tests.
- `uv run python run_pytest.py tests/api/equilibrium -q`: pass, 1 test.

## Route Guard Classification

The issue #117 guard search:

```powershell
rg "ReactivePhaseEquilibriumProblem|reactive_staged_equilibrium|staged chemical|native_derivative_free_nelder_mead|not_required.*phase" src tests docs
```

returned matches, but none are accepted reactive production paths:

- `src/epcsaft/equilibrium.py` defines `ReactivePhaseEquilibriumProblem` and routes production LLE/electrolyte LLE to `_solve_reactive_phase_equilibrium_native(...)`.
- `src/epcsaft/epcsaft.py`, `src/epcsaft/reactive_staged.py`, `src/epcsaft/reactive.py`, and `src/epcsaft/__init__.py` keep the explicitly named staged compatibility route.
- `tests/api/reactive/test_reactive_phase_equilibrium_problem_routes_native.py` and `tests/api/reactive/test_staged_reactive_route_not_production.py` intentionally monkeypatch the staged helper to fail if production calls it.
- `tests/equilibrium/reactive/test_reactive_lle.py` now exercises only explicit staged compatibility behavior.
- `native_derivative_free_nelder_mead` remains in neutral LLE diagnostics and neutral tests, not in accepted electrolyte or reactive production diagnostics.
- Older docs and goal notes are historical intake, planning, or route-guard text.

The issue #116 guard search:

```powershell
rg 'native_derivative_free_nelder_mead|not_required.*phase split solve|newton_step\(.*not_available|reactive_staged_equilibrium' src tests docs
```

returned the same classes of acceptable hits plus seed-only/optional references:

- Old electrolyte LLE helper names remain as seed-generation support, while accepted electrolyte production diagnostics are asserted as Ceres with `cppad_implicit` Jacobian.
- Tests assert the old labels are absent from accepted electrolyte production diagnostics.
- The `native_transformed_newton_seed` string is in the optional IPOPT backend seed path, not the accepted native Ceres production route.

## Stop Conditions

No T014 stop condition was triggered. The remaining guard matches are compatibility-only, historical-plan text, neutral-only diagnostics, or tests proving the old route is not production.
